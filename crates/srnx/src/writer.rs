//! Writing observation files as SRNX: the inverse of [`crate::reader::read`].
//!
//! ```no_run
//! let bytes = std::fs::read("abmf1190.26o")?;
//! let file = srnx::reader::read(&bytes)?;
//! let out = std::io::BufWriter::new(std::fs::File::create("abmf1190.srnx")?);
//! srnx::writer::write(&file, out)?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! The writer follows the reference encoder, `rnx2srnx`: it writes the
//! `SRNX`, `RHDR` and epoch chunks, then the `EVTF` chunks, then each
//! satellite's observation chunks followed by its `SATE` chunk, then the
//! `SDIR` chunk and the file digest.  It chooses delta orders and blocks to
//! make each observation chunk small, and with `SOC2` chunks, uses
//! cross-signal predictors and SSI references where they help.  Satellites
//! are encoded in parallel; the output does not depend on the number of
//! threads.

use std::collections::{btree_map, BTreeMap};
use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use thiserror::Error;

use crate::digest::DigestId;
use crate::epoc::{advance_time, Epoch, EpochChunkKind};
use crate::leb128::{sleb128_len, uleb128_len, write_sleb128, write_uleb128};
use crate::model::{EpochRecord, EventRecord, Obs, ObsCode, ObsFile, SatId};
use crate::pack::{self, Coding, Entry, Scratch};
use crate::plan::{self, Plan};
use crate::rinex2_obs::byte_at;
use crate::rinex_header::{self, Header};
use crate::soc2::{prediction, rounded, Predictor};
use crate::socd::ABSENT;
use crate::srnx_reader::{continuous_times, obs_codes_for_sat};

/// The SRNX revision 1 chunks of this writer, as the `SRNX` chunk's major
/// version.
const MAJOR_VERSION: u64 = 1;
/// Length of the `SRNX` chunk's payload: room for the version, digest
/// identifiers and any `SDIR` offset, padded with zeros.
const SRNX_PAYLOAD_LEN: usize = 16;
/// Explicit scales above this many milli-units are not used.
const MAX_SCALE: u64 = 1_000_000;
/// A `SOC2` coding is not coded in full if the other's estimate beats it
/// by more than this fraction.
const PREDICTOR_MARGIN: f64 = 0.05;
/// The LLI and SSI alphabet.
const INDICATORS: &[u8; 11] = b" 0123456789";

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum WriteError {
    #[error(transparent)]
    Io(#[from] std::io::Error),
    /// The model cannot be represented, for example an observation code
    /// that is not in the header for its system.
    #[error("{0}")]
    Invalid(String),
}

fn invalid(message: impl Into<String>) -> WriteError {
    WriteError::Invalid(message.into())
}

/// Which chunks hold the observations.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum ObsChunk {
    /// `SOC2` (SRNX revision 1.1): cross-signal prediction and split-Rice
    /// blocks, for smaller files.
    #[default]
    Soc2,
    /// `SOCD` (SRNX revision 1.0), for readers of the original revision.
    Socd,
}

/// A per-chunk or file-level digest.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Digest {
    None,
    Crc32c,
    /// BLAKE3, 16-byte output.
    Blake3_16,
    /// BLAKE3, 32-byte output.
    Blake3_32,
}

impl Digest {
    fn id(self) -> DigestId {
        match self {
            Self::None => DigestId::Null,
            Self::Crc32c => DigestId::Crc32c,
            Self::Blake3_16 => DigestId::Blake3_16,
            Self::Blake3_32 => DigestId::Blake3_32,
        }
    }
}

/// How to write an SRNX file.  The default matches `rnx2srnx`: `SOC2`
/// chunks, CRC32C chunk digests and a 32-byte BLAKE3 file digest.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct WriteOptions {
    pub obs_chunk: ObsChunk,
    pub chunk_digest: Digest,
    pub file_digest: Digest,
    /// Threads to encode satellites with; 0 (the default) uses the
    /// available parallelism.
    pub threads: usize,
}

impl Default for WriteOptions {
    fn default() -> Self {
        Self {
            obs_chunk: ObsChunk::Soc2,
            chunk_digest: Digest::Crc32c,
            file_digest: Digest::Blake3_32,
            threads: 0,
        }
    }
}

/// Writes `file` as SRNX to `out`, with the default [`WriteOptions`].
pub fn write<W: Write>(file: &ObsFile, out: W) -> Result<(), WriteError> {
    write_with(file, &WriteOptions::default(), out)
}

/// Writes `file` as SRNX to `out`.
///
/// The file is encoded in memory and then written with one `write_all`,
/// so nothing is written if the model cannot be represented.  `out` is not
/// flushed.
pub fn write_with<W: Write>(
    file: &ObsFile,
    options: &WriteOptions,
    mut out: W,
) -> Result<(), WriteError> {
    out.write_all(&encode(file, options)?)?;
    Ok(())
}

/// Append a chunk: FOURCC, ULEB128 payload length, payload and digest.
fn chunk(out: &mut Vec<u8>, fourcc: [u8; 4], payload: &[u8], digest: DigestId) {
    let start = out.len();
    out.extend_from_slice(&fourcc);
    write_uleb128(out, payload.len() as u64);
    out.extend_from_slice(payload);
    append_digest(out, start, digest);
}

/// Append the digest of `out[start..]`.
fn append_digest(out: &mut Vec<u8>, start: usize, digest: DigestId) {
    let mut value = Vec::with_capacity(digest.byte_len());
    digest.append(&out[start..], &mut value);
    out.extend_from_slice(&value);
}

/// Everything a satellite's encoding needs from the rest of the file.
struct Context<'a> {
    epoch_count: usize,
    /// Continuous epoch times, in units of `1 / units_per_second` s.
    times: &'a [i128],
    units_per_second: u64,
    obs_chunk: ObsChunk,
    digest: DigestId,
    glonass: [Option<i8>; 100],
}

/// One satellite's observations, by index into its system's code list.
struct Satellite<'a> {
    id: SatId,
    name: [u8; 3],
    codes: &'a [ObsCode],
    series: Vec<(usize, &'a [(usize, Obs)])>,
    /// Epochs that list the satellite, from `ObsFile::satellite_epochs`.
    listed: &'a [usize],
}

/// A satellite's encoded chunks, ending with its `SATE` chunk.
struct Encoded {
    bytes: Vec<u8>,
    sate_offset: usize,
    any_soc2: bool,
}

fn encode(file: &ObsFile, options: &WriteOptions) -> Result<Vec<u8>, WriteError> {
    let digest = options.chunk_digest.id();
    let rhdr = rhdr_text(file.header.as_str())?;
    let header =
        rinex_header::parse(rhdr.as_bytes()).map_err(|e| invalid(format!("RINEX header: {e}")))?;
    if !matches!(header.version.major, 2..=4) {
        return Err(invalid("RINEX header is not version 2, 3 or 4"));
    }

    let kind = if file.epochs.iter().any(|e| e.extra_ps != 0) {
        EpochChunkKind::Epo2
    } else {
        EpochChunkKind::Epoc
    };
    let units_per_second = kind.units_per_second();
    let epochs = file
        .epochs
        .iter()
        .enumerate()
        .map(|(i, e)| to_epoch(i, e))
        .collect::<Result<Vec<_>, _>>()?;
    let times = continuous_times(&epochs, units_per_second);
    if kind == EpochChunkKind::Epo2 {
        if let (Some(min), Some(max)) = (times.iter().min(), times.iter().max()) {
            if max - min >= 1 << 63 {
                return Err(invalid("epochs span too long for an EPO2 chunk"));
            }
        }
    }

    let mut out = vec![0; 5 + SRNX_PAYLOAD_LEN + digest.byte_len()];
    chunk(&mut out, *b"RHDR", rhdr.as_bytes(), digest);
    let epoch_offset = out.len();
    chunk(
        &mut out,
        kind.fourcc(),
        &epoch_payload(&epochs, units_per_second),
        digest,
    );
    let evtf_offset = write_events(&mut out, &file.events, &header, epochs.len(), digest)?;

    let satellites = satellites(file, &header)?;
    let context = Context {
        epoch_count: epochs.len(),
        times: &times,
        units_per_second,
        obs_chunk: options.obs_chunk,
        digest,
        glonass: plan::glonass_channels(&rhdr),
    };
    let mut sdir = Vec::new();
    write_uleb128(&mut sdir, epoch_offset as u64);
    write_uleb128(&mut sdir, evtf_offset as u64);
    let mut any_soc2 = false;
    for (satellite, encoded) in
        satellites
            .iter()
            .zip(encode_satellites(&satellites, &context, options.threads))
    {
        let encoded = encoded?;
        sdir.extend_from_slice(&satellite.name);
        write_uleb128(&mut sdir, (out.len() + encoded.sate_offset) as u64);
        out.extend_from_slice(&encoded.bytes);
        any_soc2 |= encoded.any_soc2;
    }
    let sdir_offset = out.len();
    chunk(&mut out, *b"SDIR", &sdir, digest);

    let minor = match (kind, any_soc2) {
        (EpochChunkKind::Epo2, _) => 2,
        (_, true) => 1,
        _ => 0,
    };
    let mut payload = Vec::with_capacity(SRNX_PAYLOAD_LEN);
    for value in [
        MAJOR_VERSION,
        minor,
        digest.raw(),
        options.file_digest.id().raw(),
        sdir_offset as u64,
    ] {
        write_uleb128(&mut payload, value);
    }
    payload.resize(SRNX_PAYLOAD_LEN, 0);
    let mut head = Vec::new();
    chunk(&mut head, *b"SRNX", &payload, digest);
    out[..head.len()].copy_from_slice(&head);

    append_digest(&mut out, 0, options.file_digest.id());
    Ok(out)
}

/// The header text from its `RINEX VERSION / TYPE` line (after any
/// Compact RINEX lines), ending with a newline.
fn rhdr_text(text: &str) -> Result<String, WriteError> {
    let mut start = 0;
    for line in text.split_inclusive('\n') {
        if line.get(60..80) == Some("RINEX VERSION / TYPE") {
            let mut rhdr = text[start..].replace("\r\n", "\n");
            if !rhdr.ends_with('\n') {
                rhdr.push('\n');
            }
            return Ok(rhdr);
        }
        start += line.len();
    }
    Err(invalid("RINEX header has no RINEX VERSION / TYPE record"))
}

fn to_epoch(i: usize, e: &EpochRecord) -> Result<Epoch, WriteError> {
    let bad = |what: &str| invalid(format!("epoch {i}: {what}"));
    if e.flag != 0 {
        return Err(bad(&format!(
            "flag {} (only flag-0 epochs are observation epochs; others are events)",
            e.flag
        )));
    }
    if e.year < 100 {
        return Err(bad("two-digit year"));
    }
    if !(1..=12).contains(&e.month) || !(1..=31).contains(&e.day) {
        return Err(bad("invalid date"));
    }
    if e.hour > 23 || e.minute > 59 || e.seconds_e7 >= 1_000_000_000 {
        return Err(bad("invalid time of day"));
    }
    if e.extra_ps >= 100_000 {
        return Err(bad("picoseconds out of range"));
    }
    Ok(Epoch {
        year: e.year,
        month: e.month,
        day: e.day,
        hour: e.hour,
        minute: e.minute,
        seconds_e7: e.seconds_e7,
        extra_ps: e.extra_ps,
        clock_offset_e12: e.clock_offset_e12,
    })
}

/// The `EPOC` or `EPO2` payload: epoch count, spans, clock offsets.
fn epoch_payload(epochs: &[Epoch], units: u64) -> Vec<u8> {
    let date = |e: &Epoch| u64::from(e.year) * 10_000 + u64::from(e.month) * 100 + u64::from(e.day);
    let time_of_day = |e: &Epoch| {
        (u64::from(e.hour) * 3600 + u64::from(e.minute) * 60) * units + e.seconds_in(units)
    };
    let mut payload = Vec::new();
    write_uleb128(&mut payload, epochs.len() as u64);
    let mut start = 0;
    while start < epochs.len() {
        let e0 = &epochs[start];
        // The interval from the first two epochs, if on the same day and
        // forward in time; the span extends while the decoder's time
        // arithmetic gives the next epoch.
        let mut interval = 0;
        let mut len = 1;
        if let Some(e1) = epochs.get(start + 1).filter(|e1| date(e1) == date(e0)) {
            interval = i64::try_from(time_of_day(e1))
                .unwrap()
                .saturating_sub(i64::try_from(time_of_day(e0)).unwrap())
                .max(0);
        }
        if interval > 0 {
            let (mut hour, mut minute, mut seconds) = (e0.hour, e0.minute, e0.seconds_in(units));
            for e in &epochs[start + 1..] {
                advance_time(&mut hour, &mut minute, &mut seconds, interval, units);
                if date(e) != date(e0)
                    || (hour, minute, seconds) != (e.hour, e.minute, e.seconds_in(units))
                {
                    break;
                }
                len += 1;
            }
        }
        if len == 1 {
            interval = 0;
        }
        let units_i = i64::try_from(units).unwrap();
        let encoded = if interval != 0 && interval % units_i == 0 {
            -(interval / units_i)
        } else {
            interval
        };
        write_sleb128(&mut payload, encoded);
        write_uleb128(&mut payload, len as u64 - 1);
        write_uleb128(&mut payload, date(e0));
        let (h, m) = (u64::from(e0.hour), u64::from(e0.minute));
        write_uleb128(
            &mut payload,
            (h * 100 + m) * 100 * units + e0.seconds_in(units),
        );
        start += len;
    }

    // Receiver clock offsets, run-length coded; trailing zeros are implied.
    let offsets: Vec<i64> = epochs.iter().map(|e| e.clock_offset_e12).collect();
    let end = offsets.iter().rposition(|&o| o != 0).map_or(0, |i| i + 1);
    let mut i = 0;
    while i < end {
        let run = offsets[i..end]
            .iter()
            .take_while(|&&o| o == offsets[i])
            .count();
        write_sleb128(&mut payload, offsets[i]);
        write_uleb128(&mut payload, run as u64 - 1);
        i += run;
    }
    payload
}

/// Append one `EVTF` chunk per event; returns the offset of the first, or
/// 0 if there are none.
fn write_events(
    out: &mut Vec<u8>,
    events: &[(usize, EventRecord)],
    header: &Header,
    epoch_count: usize,
    digest: DigestId,
) -> Result<usize, WriteError> {
    let flag_col = if header.version.is_v3_or_later() {
        31
    } else {
        28
    };
    let mut first = 0;
    for (k, (index, event)) in events.iter().enumerate() {
        let bad = |what: &str| invalid(format!("event {k}: {what}"));
        if *index > epoch_count {
            return Err(bad("epoch index is past the last epoch"));
        }
        if !(1..=6).contains(&event.flag) {
            return Err(bad("flag is not 1 to 6"));
        }
        if byte_at(&event.epoch_line, flag_col) != b'0' + event.flag {
            return Err(bad("epoch line does not have the event's flag"));
        }
        let mut payload = Vec::new();
        write_uleb128(&mut payload, *index as u64);
        for line in std::iter::once(&event.epoch_line).chain(&event.lines) {
            let line = line.strip_suffix('\r').unwrap_or(line);
            if line.contains(['\n', '\r']) {
                return Err(bad("line contains a line break"));
            }
            payload.extend_from_slice(line.as_bytes());
            payload.push(b'\n');
        }
        if k == 0 {
            first = out.len();
        }
        chunk(out, *b"EVTF", &payload, digest);
    }
    Ok(first)
}

/// Group the observations by satellite, checking each code against the
/// header's list for the satellite's system.
fn satellites<'a>(file: &'a ObsFile, header: &'a Header) -> Result<Vec<Satellite<'a>>, WriteError> {
    let mut satellites: BTreeMap<SatId, Satellite<'a>> = BTreeMap::new();
    let codes = |sat: SatId| obs_codes_for_sat(header, sat.system);
    for (&sat, epochs) in &file.satellite_epochs {
        if !epochs.is_empty() {
            satellite(&mut satellites, sat, codes(sat))?.listed = epochs;
        }
    }
    for (&(sat, code), series) in &file.observations {
        if series.is_empty() {
            continue;
        }
        let satellite = satellite(&mut satellites, sat, codes(sat))?;
        let index = satellite
            .codes
            .iter()
            .position(|&c| c == code)
            .ok_or_else(|| {
                invalid(format!(
                    "{sat}: observation code {} is not in the header for system {}",
                    code.as_str(),
                    sat.system.as_char()
                ))
            })?;
        satellite.series.push((index, series));
    }
    Ok(satellites.into_values().collect())
}

/// The entry for `sat` in `satellites`, added if new.
fn satellite<'s, 'a>(
    satellites: &'s mut BTreeMap<SatId, Satellite<'a>>,
    sat: SatId,
    codes: &'a [ObsCode],
) -> Result<&'s mut Satellite<'a>, WriteError> {
    match satellites.entry(sat) {
        btree_map::Entry::Occupied(entry) => Ok(entry.into_mut()),
        btree_map::Entry::Vacant(entry) => {
            let name = format!("{sat}");
            let name: [u8; 3] = name
                .as_bytes()
                .try_into()
                .ok()
                .filter(|n: &[u8; 3]| n[0].is_ascii_uppercase())
                .ok_or_else(|| invalid(format!("satellite {name:?} has no SRNX name")))?;
            Ok(entry.insert(Satellite {
                id: sat,
                name,
                codes,
                series: Vec::new(),
                listed: &[],
            }))
        }
    }
}

/// Encode `satellites` on up to `threads` threads (0 for the available
/// parallelism), returning the results in order.
fn encode_satellites(
    satellites: &[Satellite<'_>],
    context: &Context<'_>,
    threads: usize,
) -> Vec<Result<Encoded, WriteError>> {
    let threads = if threads == 0 {
        std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get)
    } else {
        threads
    }
    .min(satellites.len());
    if threads <= 1 {
        let mut scratch = Scratch::default();
        return satellites
            .iter()
            .map(|s| encode_satellite(s, context, &mut scratch))
            .collect();
    }
    let results: Vec<Mutex<Option<Result<Encoded, WriteError>>>> =
        satellites.iter().map(|_| Mutex::new(None)).collect();
    let next = AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| {
                let mut scratch = Scratch::default();
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    let Some(satellite) = satellites.get(i) else {
                        break;
                    };
                    let result = encode_satellite(satellite, context, &mut scratch);
                    *results[i].lock().unwrap() = Some(result);
                }
            });
        }
    });
    results
        .into_iter()
        .map(|r| r.into_inner().unwrap().unwrap())
        .collect()
}

/// One observation code's values over a satellite's epochs.
struct Column {
    /// Milli-units, `ABSENT` where missing.
    values: Vec<i64>,
    lli: Vec<u8>,
    ssi: Vec<u8>,
}

/// The columns of `satellite` over its satellite epochs `epochs`, by
/// index into its code list.
fn columns(satellite: &Satellite<'_>, epochs: &[usize]) -> Result<Vec<Option<Column>>, WriteError> {
    let sat = satellite.id;
    let mut columns: Vec<Option<Column>> = satellite.codes.iter().map(|_| None).collect();
    for &(index, series) in &satellite.series {
        let code = satellite.codes[index].as_str();
        let n = epochs.len();
        let mut column = Column {
            values: vec![ABSENT; n],
            lli: vec![b' '; n],
            ssi: vec![b' '; n],
        };
        let mut k = 0;
        let mut previous = None;
        for &(epoch, obs) in series {
            if previous.is_some_and(|p| p >= epoch) {
                return Err(invalid(format!(
                    "{sat}/{code}: epoch indices are not strictly increasing"
                )));
            }
            previous = Some(epoch);
            while epochs[k] < epoch {
                k += 1;
            }
            if obs.value == ABSENT {
                return Err(invalid(format!("{sat}/{code}: value {ABSENT} is reserved")));
            }
            for (what, ch) in [("LLI", obs.lli), ("SSI", obs.ssi)] {
                if !INDICATORS.contains(&ch) {
                    return Err(invalid(format!(
                        "{sat}/{code}: {what} {:?} is not a blank or digit",
                        char::from(ch)
                    )));
                }
            }
            column.values[k] = obs.value;
            column.lli[k] = obs.lli;
            column.ssi[k] = obs.ssi;
        }
        columns[index] = Some(column);
    }
    Ok(columns)
}

fn encode_satellite(
    satellite: &Satellite<'_>,
    context: &Context<'_>,
    scratch: &mut Scratch,
) -> Result<Encoded, WriteError> {
    let sat = satellite.id;
    // The satellite epochs: every epoch that lists the satellite or has
    // one of its observations.
    let mut listed = vec![false; context.epoch_count];
    let all = satellite
        .series
        .iter()
        .flat_map(|(_, series)| series.iter().map(|&(e, _)| e));
    for e in all.chain(satellite.listed.iter().copied()) {
        *listed
            .get_mut(e)
            .ok_or_else(|| invalid(format!("{sat}: epoch index past the last epoch")))? = true;
    }
    let epochs: Vec<usize> = (0..listed.len()).filter(|&e| listed[e]).collect();
    let columns = columns(satellite, &epochs)?;

    let times: Vec<i128> = epochs.iter().map(|&e| context.times[e]).collect();
    let mut bytes = Vec::new();
    let mut offsets: Vec<Option<usize>> = vec![None; columns.len()];
    let mut any_soc2 = false;
    let name = |code: ObsCode| {
        let mut name = [0u8; 8];
        name[..3].copy_from_slice(&satellite.name);
        let code = code.as_str().as_bytes();
        name[4..4 + code.len()].copy_from_slice(code);
        name
    };
    match context.obs_chunk {
        ObsChunk::Socd => {
            for (j, column) in columns.iter().enumerate() {
                if let Some(column) = column {
                    offsets[j] = Some(bytes.len());
                    let payload = socd_payload(name(satellite.codes[j]), column, scratch);
                    chunk(&mut bytes, *b"SOCD", &payload, context.digest);
                }
            }
        }
        ObsChunk::Soc2 => {
            let values: Vec<Option<&[i64]>> = columns
                .iter()
                .map(|c| c.as_ref().map(|c| &c.values[..]))
                .collect();
            let plans = plan::plan_satellite(
                sat,
                satellite.codes,
                &values,
                &context.glonass,
                context.units_per_second,
            );
            let coder = Soc2Coder {
                columns: &columns,
                times: &times,
                units_per_second: context.units_per_second,
            };
            // Referenced chunks first: unpredicted ones (the reference
            // phases), then predicted phases, then the rest.
            for pass in 0..3 {
                for (j, column) in columns.iter().enumerate() {
                    let predicted = plans[j].predictor != Predictor::None;
                    let phase = satellite.codes[j].as_str().starts_with('L');
                    if column.is_none()
                        || (pass == 0) == predicted
                        || (pass == 1 && !phase)
                        || (pass == 2 && phase)
                    {
                        continue;
                    }
                    offsets[j] = Some(bytes.len());
                    let payload = coder.payload(j, name(satellite.codes[j]), plans[j], scratch);
                    chunk(&mut bytes, *b"SOC2", &payload, context.digest);
                    any_soc2 = true;
                }
            }
        }
    }

    let sate_offset = bytes.len();
    let payload = sate_payload(satellite.name, &offsets, sate_offset, &epochs);
    chunk(&mut bytes, *b"SATE", &payload, context.digest);
    Ok(Encoded {
        bytes,
        sate_offset,
        any_soc2,
    })
}

/// The `SATE` payload: the satellite's name, the offsets of its
/// observation chunks relative to the `SATE` chunk at `sate_offset`, and
/// its satellite epochs as runs.
fn sate_payload(
    name: [u8; 3],
    offsets: &[Option<usize>],
    sate_offset: usize,
    epochs: &[usize],
) -> Vec<u8> {
    let mut payload = name.to_vec();
    payload.push(0);
    for offset in offsets {
        let relative = offset.map_or(0, |o| {
            i64::try_from(o).unwrap() - i64::try_from(sate_offset).unwrap()
        });
        write_sleb128(&mut payload, relative);
    }
    // Every satellite has at least one epoch (it is in the model only
    // through a series or a presence list), so the run count-minus-1 exists.
    debug_assert!(!epochs.is_empty(), "a satellite has at least one epoch");
    // Epoch presence: runs of consecutive epochs.
    let mut runs = Vec::new();
    for &e in epochs {
        match runs.last_mut() {
            Some((_, end)) if *end == e => *end += 1,
            _ => runs.push((e, e + 1)),
        }
    }
    write_uleb128(&mut payload, runs.len() as u64 - 1);
    let mut previous_end = 0;
    for &(start, end) in &runs {
        write_uleb128(&mut payload, (start - previous_end) as u64);
        write_uleb128(&mut payload, (end - start - 1) as u64);
        previous_end = end;
    }
    payload
}

fn indicator_symbol(ch: u8) -> u8 {
    u8::try_from(INDICATORS.iter().position(|&c| c == ch).unwrap()).unwrap()
}

/// Append a ULEB128-length-prefixed run-length coding of `symbols` over
/// an alphabet of `alphabet` symbols, leaving out trailing blanks (symbol 0).
fn write_rle(out: &mut Vec<u8>, symbols: impl Iterator<Item = u8>, alphabet: u64) {
    let symbols: Vec<u8> = symbols.collect();
    let end = symbols.iter().rposition(|&s| s != 0).map_or(0, |i| i + 1);
    let mut rle = Vec::new();
    let mut i = 0;
    while i < end {
        let run = symbols[i..end]
            .iter()
            .take_while(|&&s| s == symbols[i])
            .count();
        write_uleb128(
            &mut rle,
            (run as u64 - 1) * alphabet + u64::from(symbols[i]),
        );
        i += run;
    }
    write_uleb128(out, rle.len() as u64);
    out.extend_from_slice(&rle);
}

/// The explicit scale for `values`: the greatest common divisor of the
/// present values, or 1 if that is 0 or too large.
fn scale(values: &[i64]) -> u64 {
    let mut gcd = 0u64;
    for &v in values {
        if v == ABSENT {
            continue;
        }
        let mut b = v.unsigned_abs();
        while b != 0 {
            (gcd, b) = (b, gcd % b);
        }
        if gcd == 1 {
            break;
        }
    }
    if gcd == 0 || gcd > MAX_SCALE {
        1
    } else {
        gcd
    }
}

/// Divide a present value by the column's scale (which divides it).
fn scaled(value: i64, scale: u64) -> i64 {
    if scale == 1 {
        value
    } else {
        value / i64::try_from(scale).unwrap()
    }
}

fn socd_payload(name: [u8; 8], column: &Column, scratch: &mut Scratch) -> Vec<u8> {
    // A column exists only for a satellite with at least one epoch.
    debug_assert!(!column.values.is_empty(), "a column has at least one epoch");
    let mut payload = name.to_vec();
    write_uleb128(&mut payload, column.values.len() as u64 - 1);
    write_rle(
        &mut payload,
        column.lli.iter().map(|&c| indicator_symbol(c)),
        11,
    );
    write_rle(
        &mut payload,
        column.ssi.iter().map(|&c| indicator_symbol(c)),
        11,
    );
    let scale = scale(&column.values);
    let entries: Vec<Entry> = column
        .values
        .iter()
        .map(|&v| {
            if v == ABSENT {
                Entry::Absent
            } else {
                Entry::Coded(scaled(v, scale))
            }
        })
        .collect();
    pack::encode(&mut payload, &entries, scale, Coding::Matrix, scratch);
    payload
}

/// [`Soc2Coder::entries`] for a value predictor: `reference × num / den`,
/// rounded as [`rounded`] does, in 64 bits where the product fits.
fn value_entries(
    values: &[i64],
    reference: &[i64],
    scale: u64,
    num: i64,
    den: u64,
) -> Option<Vec<Entry>> {
    let den64 = i64::try_from(den).ok();
    let mut entries = Vec::with_capacity(values.len());
    for (&v, &r) in values.iter().zip(reference) {
        entries.push(if v == ABSENT {
            Entry::Absent
        } else if r == ABSENT {
            Entry::Unpredicted(scaled(v, scale))
        } else {
            let p = match (r.checked_mul(num), den64) {
                (Some(product), Some(den)) => i128::from(rounded64(product, den)),
                _ => rounded(i128::from(r) * i128::from(num), i128::from(den)).ok()?,
            };
            let q = if scale == 1 {
                p
            } else {
                rounded(p, i128::from(scale)).ok()?
            };
            Entry::Coded(i64::try_from(i128::from(scaled(v, scale)) - q).ok()?)
        });
    }
    Some(entries)
}

/// [`rounded`] in 64 bits: `num / den` to the nearest integer, ties away
/// from zero; `den` is positive.
fn rounded64(num: i64, den: i64) -> i64 {
    let (q, rem) = (num / den, num % den);
    if rem.abs() >= den / 2 + den % 2 {
        q + num.signum()
    } else {
        q
    }
}

/// Codes one satellite's `SOC2` chunks.
struct Soc2Coder<'a> {
    columns: &'a [Option<Column>],
    times: &'a [i128],
    units_per_second: u64,
}

/// The SSI digit `SOC2` predicts from a signal strength in milli-units.
fn predicted_ssi(strength: i64) -> u8 {
    b'0' + u8::try_from(strength.div_euclid(6000).clamp(1, 9)).unwrap()
}

/// Bytes of a predictor's fields (its kind is in the mode value).
fn predictor_bytes(predictor: Predictor) -> usize {
    match predictor {
        Predictor::None => 0,
        Predictor::Value {
            reference,
            num,
            den,
        } => uleb128_len(reference as u64) + sleb128_len(num) + uleb128_len(den),
        Predictor::Rate {
            reference,
            num,
            den,
            stencil,
            max_span,
        } => {
            uleb128_len(reference as u64)
                + sleb128_len(num)
                + uleb128_len(den)
                + uleb128_len(u64::from(stencil))
                + uleb128_len(max_span)
        }
    }
}

impl Soc2Coder<'_> {
    fn column(&self, j: usize) -> &Column {
        self.columns[j].as_ref().unwrap()
    }

    /// The entries of column `j` at `scale` with `predictor`, or `None` if
    /// a residual does not fit or the prediction cannot be computed.
    fn entries(&self, j: usize, scale: u64, predictor: Predictor) -> Option<Vec<Entry>> {
        let values = &self.column(j).values;
        let reference = match predictor {
            Predictor::None => None,
            Predictor::Value { reference, .. } | Predictor::Rate { reference, .. } => {
                Some(&self.column(reference).values)
            }
        };
        if let (Predictor::Value { num, den, .. }, Some(reference)) = (predictor, reference) {
            return value_entries(values, reference, scale, num, den);
        }
        let milli = |i: usize| {
            reference
                .and_then(|r| r.get(i))
                .filter(|&&v| v != ABSENT)
                .map(|&v| i128::from(v))
        };
        values
            .iter()
            .enumerate()
            .map(|(k, &v)| {
                if v == ABSENT {
                    return Some(Entry::Absent);
                }
                let x = scaled(v, scale);
                if predictor == Predictor::None {
                    return Some(Entry::Coded(x));
                }
                match prediction(predictor, k, milli, self.times, self.units_per_second).ok()? {
                    None => Some(Entry::Unpredicted(x)),
                    Some(p) => {
                        let q = rounded(p, i128::from(scale)).ok()?;
                        let residual = i64::try_from(i128::from(x) - q).ok()?;
                        Some(Entry::Coded(residual))
                    }
                }
            })
            .collect()
    }

    fn payload(&self, j: usize, name: [u8; 8], plan: Plan, scratch: &mut Scratch) -> Vec<u8> {
        let column = self.column(j);
        let scale = scale(&column.values);

        // SSIs: plain, or against the planned signal-strength code.
        let symbols: Vec<u8> = column.ssi.iter().map(|&c| indicator_symbol(c)).collect();
        let mut ssi = (0, symbols.clone(), 11);
        let mut unpredicted_rle = Vec::new();
        write_rle(&mut unpredicted_rle, symbols.iter().copied(), 11);
        if let Some(s) = plan.ssi_reference {
            let strengths = &self.column(s).values;
            let predicted: Vec<u8> = symbols
                .iter()
                .zip(&column.ssi)
                .zip(strengths)
                .map(|((&sym, &ch), &strength)| {
                    if strength != ABSENT && ch == predicted_ssi(strength) {
                        11
                    } else {
                        sym
                    }
                })
                .collect();
            let mut predicted_rle = Vec::new();
            write_rle(&mut predicted_rle, predicted.iter().copied(), 12);
            let mode_len = uleb128_len((s as u64 + 1) << 2);
            if mode_len + predicted_rle.len() < 1 + unpredicted_rle.len() {
                ssi = (s + 1, predicted, 12);
            }
        }

        // Values: without and (if planned) with the predictor; keep the
        // smaller.  A coding whose estimate is clearly worse is skipped.
        let mut candidates = vec![(
            Predictor::None,
            self.entries(j, scale, Predictor::None).unwrap(),
        )];
        if plan.predictor != Predictor::None {
            if let Some(entries) = self.entries(j, scale, plan.predictor) {
                candidates.push((plan.predictor, entries));
            }
        }
        if candidates.len() == 2 {
            #[allow(clippy::cast_precision_loss)]
            let estimates: Vec<f64> = candidates
                .iter()
                .map(|(p, e)| pack::estimate(e, scratch) + predictor_bytes(*p) as f64)
                .collect();
            if estimates[1] < estimates[0] * (1.0 - PREDICTOR_MARGIN) {
                candidates.remove(0);
            } else if estimates[0] < estimates[1] * (1.0 - PREDICTOR_MARGIN) {
                candidates.pop();
            }
        }
        let (predictor, packed) = candidates
            .iter()
            .map(|(p, entries)| {
                let mut packed = Vec::new();
                pack::encode(&mut packed, entries, scale, Coding::Rice, scratch);
                (*p, packed)
            })
            .min_by_key(|(p, packed)| predictor_bytes(*p) + packed.len())
            .unwrap();

        debug_assert!(!column.values.is_empty(), "a column has at least one epoch");
        let mut payload = name.to_vec();
        write_uleb128(&mut payload, column.values.len() as u64 - 1);
        write_rle(
            &mut payload,
            column.lli.iter().map(|&c| indicator_symbol(c)),
            11,
        );
        let kind = match predictor {
            Predictor::None => 0,
            Predictor::Value { .. } => 1,
            Predictor::Rate { .. } => 2,
        };
        write_uleb128(&mut payload, (ssi.0 as u64) << 2 | kind);
        write_rle(&mut payload, ssi.1.into_iter(), ssi.2);
        match predictor {
            Predictor::None => {}
            Predictor::Value {
                reference,
                num,
                den,
            } => {
                write_uleb128(&mut payload, reference as u64);
                write_sleb128(&mut payload, num);
                write_uleb128(&mut payload, den);
            }
            Predictor::Rate {
                reference,
                num,
                den,
                stencil,
                max_span,
            } => {
                write_uleb128(&mut payload, reference as u64);
                write_sleb128(&mut payload, num);
                write_uleb128(&mut payload, den);
                write_uleb128(&mut payload, u64::from(stencil));
                write_uleb128(&mut payload, max_span);
            }
        }
        payload.extend_from_slice(&packed);
        payload
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::epoc;

    fn epoch(hour: u8, minute: u8, seconds_e7: u64, extra_ps: u32) -> Epoch {
        Epoch {
            year: 2016,
            month: 12,
            day: 31,
            hour,
            minute,
            seconds_e7,
            extra_ps,
            clock_offset_e12: 0,
        }
    }

    /// The epochs decode back, and the payload has `spans` spans.
    fn check(epochs: &[Epoch], kind: EpochChunkKind, spans: usize) {
        let payload = epoch_payload(epochs, kind.units_per_second());
        assert_eq!(epoc::parse(&payload, kind).unwrap().epochs, epochs);
        let mut cursor = &payload[..];
        crate::leb128::read_uleb128(&mut cursor).unwrap();
        let mut found = 0;
        let mut covered = 0;
        while covered < epochs.len() {
            crate::leb128::read_sleb128(&mut cursor).unwrap();
            covered +=
                1 + usize::try_from(crate::leb128::read_uleb128(&mut cursor).unwrap()).unwrap();
            crate::leb128::read_uleb128(&mut cursor).unwrap();
            crate::leb128::read_uleb128(&mut cursor).unwrap();
            found += 1;
        }
        assert_eq!(found, spans);
    }

    #[test]
    fn spans_carry_minutes_and_hours() {
        let every = |step: u64, n: u64| -> Vec<Epoch> {
            (0..n)
                .map(|i| {
                    let s = 23 * 3600 + i * step;
                    let (minute, second) = ((s / 60) % 60, s % 60);
                    epoch(23, u8::try_from(minute).unwrap(), second * 10_000_000, 0)
                })
                .collect()
        };
        // 30 s data: one span, over minute boundaries.
        check(&every(30, 100), EpochChunkKind::Epoc, 1);
        // Steps that carry several minutes, or into the next hour.
        check(&every(90, 4), EpochChunkKind::Epoc, 1);
        check(&every(300, 12), EpochChunkKind::Epoc, 1);
        let hourly: Vec<Epoch> = (0..5).map(|h| epoch(h * 5, 0, 0, 0)).collect();
        check(&hourly, EpochChunkKind::Epoc, 1);
        // A leap second stays in the day's last minute.
        let leap = [
            epoch(23, 59, 590_000_000, 0),
            epoch(23, 59, 600_000_000, 0),
            epoch(23, 59, 605_000_000, 0),
        ];
        check(&leap, EpochChunkKind::Epoc, 2);
        // Picosecond intervals.
        let ps = [
            epoch(13, 10, 360_000_000, 12_345),
            epoch(13, 10, 540_000_000, 0),
            epoch(13, 11, 120_000_000, 0),
            epoch(13, 14, 480_000_000, 0),
        ];
        check(&ps, EpochChunkKind::Epo2, 2);
    }
}
