//! Converts a parsed SRNX file into the shared [`ObsFile`] semantic model.
//!
//! Entry points: [`read`] and [`catalog`].

use std::collections::{BTreeMap, HashSet};

use thiserror::Error;

use crate::digest::DigestId;
use crate::epoc::{Epoch, EpochChunkKind};
use crate::file::{FileError, SrnxFile};
use crate::model::{
    Catalog, EpochRecord, EventRecord, Obs, ObsFile, SatId, SignalFilter, SignalMasks,
};
use crate::rinex2_obs::byte_at;
use crate::rinex_header::{self, Header, HeaderError, HeaderText, ObsCode, SatSystem};
use crate::sdir::Sdir;
use crate::soc2::{Column, SatelliteColumns};
use crate::socd::{self, ABSENT};

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum SrnxReadError {
    #[error("SRNX file error: {0}")]
    File(#[from] FileError),
    #[error("RINEX header error: {0}")]
    Header(#[from] HeaderError),
    #[error("file-level digest mismatch")]
    FileDigestMismatch,
    #[error("SRNX chunk digest mismatch")]
    SrnxChunkDigestMismatch,
    #[error("per-chunk digest mismatch at file offset {0}")]
    ChunkDigestMismatch(u64),
    #[error("no RHDR chunk in SRNX file")]
    NoRhdr,
    #[error("no SDIR chunk in SRNX file")]
    NoSdir,
    #[error("no EPOC or EPO2 chunk in SRNX file")]
    NoEpoc,
    #[error("EPO2 chunk spans 2^63 picoseconds or more")]
    EpochRangeTooLarge,
    #[error("scaled observation value does not fit in 64 bits")]
    ValueOverflow,
    #[error("EVTF chunk {index} (from 0): {reason}")]
    BadEvent { index: usize, reason: &'static str },
    #[error("satellite name {0:?} is not a RINEX satellite identifier")]
    BadSatelliteName(String),
    #[error("SDIR names satellite {0} more than once")]
    DuplicateSatellite(SatId),
}

/// List the satellites in an SRNX file's `SDIR` chunk, each with the
/// observation codes that its `SATE` chunk gives observations for.
pub fn catalog(bytes: &[u8]) -> Result<Catalog, SrnxReadError> {
    let file = open(bytes)?;
    let rhdr_bytes = file.rhdr_payload().ok_or(SrnxReadError::NoRhdr)?;
    let header = rinex_header::parse(rhdr_bytes)?;
    let sdir = file.sdir()?.ok_or(SrnxReadError::NoSdir)?;
    check_satellites(&sdir)?;
    let mut satellites = Vec::with_capacity(sdir.satellites.len());
    for entry in &sdir.satellites {
        let sat_id = satellite_id(entry.sat_id)?;
        let obs_codes = obs_codes_for_sat(&header, sat_id.system);
        let offsets = file.sate_offsets_at(entry.sate_offset, obs_codes.len())?;
        let observed = obs_codes
            .iter()
            .zip(&offsets)
            .filter(|&(_, &offset)| offset != 0)
            .map(|(&code, _)| code)
            .collect();
        satellites.push((sat_id, observed));
    }
    Ok(Catalog {
        codes_by_system: header.obs_codes.by_system,
        satellites: Some(satellites),
    })
}

/// Parse a complete SRNX observation file from bytes and return the
/// normalized [`ObsFile`] semantic model, keeping only the satellite/signal
/// pairs that `filter` selects (all of them for `None`).
///
/// Only the selected observation chunks, and those they reference for
/// prediction, are decoded.
pub fn read(bytes: &[u8], filter: Option<SignalFilter<'_>>) -> Result<ObsFile, SrnxReadError> {
    let file = open(bytes)?;

    // RHDR → HeaderText (public) + internal Header (for obs codes).
    let rhdr_bytes = file.rhdr_payload().ok_or(SrnxReadError::NoRhdr)?;
    let header_text = HeaderText::from_bytes(rhdr_bytes);
    let header = rinex_header::parse(rhdr_bytes)?;

    // EPOC or EPO2 → epoch list.
    let epoc = file.epoc()?.ok_or(SrnxReadError::NoEpoc)?;
    let units_per_second = epoc.kind.units_per_second();
    let total_epochs = epoc.epochs.len();
    let timestamps = continuous_times(&epoc.epochs, units_per_second);
    check_epoch_range(epoc.kind, &timestamps)?;
    let epochs: Vec<EpochRecord> = epoc.epochs.iter().map(epoch_to_record).collect();

    // SDIR.
    let sdir = file.sdir()?.ok_or(SrnxReadError::NoSdir)?;
    check_satellites(&sdir)?;

    let events = events(&file, &header, sdir.evtf_offset, total_epochs)?;

    let mut observations: BTreeMap<(SatId, ObsCode), Vec<(usize, Obs)>> = BTreeMap::new();
    let mut masks = SignalMasks::new(filter);
    let mut satellite_epochs = BTreeMap::new();

    for entry in &sdir.satellites {
        let sat_id = satellite_id(entry.sat_id)?;
        let sate_file_offset = entry.sate_offset;

        // Choose obs-code list for this satellite's system.
        let obs_codes: &[ObsCode] = obs_codes_for_sat(&header, sat_id.system);
        let num_obs_codes = obs_codes.len();
        let keep = masks.get(sat_id, obs_codes);
        if keep.is_some_and(|keep| !keep.contains(&true)) {
            continue;
        }

        let sate = file.sate_at(sate_file_offset, num_obs_codes, total_epochs)?;
        let present = sate.present_epoch_indices;
        let offsets = sate.socd_offsets;

        let times: Vec<i128> = present.iter().map(|&i| timestamps[i]).collect();
        let mut load = |code_idx: usize| -> Result<Option<Column>, FileError> {
            let rel_offset = offsets[code_idx];
            if rel_offset == 0 {
                return Ok(None);
            }
            let offset =
                sate_file_offset
                    .checked_add_signed(rel_offset)
                    .ok_or(FileError::Truncated {
                        what: "observation offset out of range",
                    })?;
            let chunk = file.chunk_at_offset(offset).ok_or(FileError::Truncated {
                what: "observation chunk at advertised offset",
            })?;
            // A chunk holds one observation per satellite epoch.  Check that
            // before decoding it: the declared count is file-supplied and
            // bounds the chunk's allocations.
            let declared = socd::declared_count(chunk.payload)?;
            if usize::try_from(declared).ok() != Some(present.len()) {
                return Err(FileError::ObservationCount {
                    declared,
                    present: present.len(),
                });
            }
            match &chunk.fourcc {
                b"SOCD" => Ok(Some(Column::Socd(file.socd_at(offset)?))),
                b"SOC2" => Ok(Some(Column::Soc2(file.soc2_at(offset)?))),
                _ => Err(FileError::Truncated {
                    what: "observation FOURCC mismatch",
                }),
            }
        };
        let mut columns = SatelliteColumns::new(obs_codes, &times, units_per_second);
        for code_idx in 0..num_obs_codes {
            if keep.map_or(true, |keep| keep[code_idx]) {
                columns.decode(code_idx, &mut load)?;
            }
        }
        for (code_idx, socd) in columns.into_columns().into_iter().enumerate() {
            let Some(socd) = socd else {
                continue;
            };
            let obs_code = obs_codes[code_idx];

            let series = observations.entry((sat_id, obs_code)).or_default();
            // `SatelliteColumns` and the declared-count check above both
            // guarantee one value per satellite epoch.
            for (i, (&ep_idx, &val)) in present.iter().zip(&socd.values).enumerate() {
                if val == ABSENT {
                    continue;
                }
                let value = val
                    .checked_mul(i64::from(socd.scale))
                    .ok_or(SrnxReadError::ValueOverflow)?;
                let lli = socd.lli.get(i).copied().unwrap_or(b' ');
                let ssi = socd.ssi.get(i).copied().unwrap_or(b' ');
                series.push((ep_idx, Obs { value, lli, ssi }));
            }
        }
        // The list is the model's own; it is moved rather than cloned.
        satellite_epochs.insert(sat_id, present);
    }

    Ok(ObsFile {
        header: header_text,
        epochs,
        events,
        observations,
        satellite_epochs,
        // Malformed chunks are errors here, so nothing is skipped.
        skipped_records: 0,
    })
}

/// Parse the file framing and check the file-level digest.
/// Read the `EVTF` chunks as event records.  Each one's text starts with the
/// epoch line, whose flag (1–6) is in column 29 (v2) or 32 (v3/v4); like the
/// RINEX readers, keep only the lines after it.
fn events(
    file: &SrnxFile<'_>,
    header: &Header,
    evtf_offset: u64,
    total_epochs: usize,
) -> Result<Vec<(usize, EventRecord)>, SrnxReadError> {
    let flag_col = if header.version.is_v3_or_later() {
        31
    } else {
        28
    };
    let mut events = Vec::new();
    for (index, evtf) in file.evtf_chunks(evtf_offset)?.into_iter().enumerate() {
        let bad = |reason| SrnxReadError::BadEvent { index, reason };
        let epoch_index = usize::try_from(evtf.epoch_index)
            .ok()
            .filter(|&i| i <= total_epochs)
            .ok_or_else(|| bad("epoch index is past the last epoch"))?;
        let mut lines = evtf.text.iter().map(|l| l.trim_end_matches('\n'));
        let epoch_line = lines.next().ok_or_else(|| bad("no epoch line"))?;
        let flag = match byte_at(epoch_line, flag_col) {
            b @ b'1'..=b'6' => b - b'0',
            _ => return Err(bad("epoch flag is not 1 to 6")),
        };
        events.push((
            epoch_index,
            EventRecord {
                flag,
                epoch_line: epoch_line.to_owned(),
                lines: lines.map(str::to_owned).collect(),
            },
        ));
    }
    Ok(events)
}

/// Parse the file framing and verify whichever digests the file carries.
///
/// A valid file-level digest covers every byte, so it implies valid
/// per-chunk digests and is the only one checked.  Without one, the
/// per-chunk digests (and the SRNX chunk's own) are the only integrity
/// check the file has, so each is enforced.
fn open(bytes: &[u8]) -> Result<SrnxFile<'_>, SrnxReadError> {
    let file = SrnxFile::parse(bytes)?;
    if !file.file_digest_ok {
        return Err(SrnxReadError::FileDigestMismatch);
    }
    if file.header.file_digest == DigestId::Null {
        if !file.srnx_digest_ok {
            return Err(SrnxReadError::SrnxChunkDigestMismatch);
        }
        if let Some(chunk) = file.chunks.iter().find(|chunk| !chunk.digest_ok) {
            return Err(SrnxReadError::ChunkDigestMismatch(chunk.file_offset));
        }
    }
    Ok(file)
}

/// Check that the `SDIR` names each satellite once.  The spec says a name
/// MUST be used in at most one `SATE` chunk, and `sdir::parse` checks the
/// raw names, but the aliases that [`SatId::parse_name`] accepts (`G01`,
/// `G 1`, ` 01`) would otherwise each add the satellite's observations.
fn check_satellites(sdir: &Sdir) -> Result<(), SrnxReadError> {
    let mut seen = HashSet::with_capacity(sdir.satellites.len());
    for entry in &sdir.satellites {
        let sat = satellite_id(entry.sat_id)?;
        if !seen.insert(sat) {
            return Err(SrnxReadError::DuplicateSatellite(sat));
        }
    }
    Ok(())
}

// Gregorian day numbers, plus any leap seconds represented by the epoch
// chunk, in units of `1 / units_per_second` s. A leap second remains in its
// day's final minute and shifts all later days by 1 s.
pub(crate) fn continuous_times(epochs: &[Epoch], units_per_second: u64) -> Vec<i128> {
    let mut leap_days = std::collections::BTreeSet::new();
    let day_number = |e: &Epoch| {
        let year = i128::from(e.year) - 1;
        let mut days = year * 365 + year / 4 - year / 100 + year / 400;
        for month in 1..e.month {
            days += match month {
                4 | 6 | 9 | 11 => 30,
                2 => {
                    if e.year % 4 == 0 && (e.year % 100 != 0 || e.year % 400 == 0) {
                        29
                    } else {
                        28
                    }
                }
                _ => 31,
            };
        }
        days + i128::from(e.day) - 1
    };
    for e in epochs {
        if e.hour == 23 && e.minute == 59 && e.seconds_e7 >= 600_000_000 {
            leap_days.insert(day_number(e));
        }
    }
    epochs
        .iter()
        .map(|e| {
            let day = day_number(e);
            let leaps = i128::try_from(leap_days.range(..day).count()).unwrap();
            (day * 86_400 + leaps + i128::from(e.hour) * 3600 + i128::from(e.minute) * 60)
                * i128::from(units_per_second)
                + i128::from(e.seconds_in(units_per_second))
        })
        .collect()
}

// EPO2 requires its timestamp differences to fit in a signed 64-bit integer.
fn check_epoch_range(kind: EpochChunkKind, timestamps: &[i128]) -> Result<(), SrnxReadError> {
    if kind == EpochChunkKind::Epo2 {
        if let (Some(min), Some(max)) = (timestamps.iter().min(), timestamps.iter().max()) {
            if max - min >= 1 << 63 {
                return Err(SrnxReadError::EpochRangeTooLarge);
            }
        }
    }
    Ok(())
}

fn epoch_to_record(e: &Epoch) -> EpochRecord {
    EpochRecord {
        year: e.year,
        month: e.month,
        day: e.day,
        hour: e.hour,
        minute: e.minute,
        seconds_e7: e.seconds_e7,
        extra_ps: e.extra_ps,
        flag: 0,
        clock_offset_e12: e.clock_offset_e12,
    }
}

/// A satellite name from an `SDIR` entry, which holds the RINEX identifier.
fn satellite_id(name: [u8; 3]) -> Result<SatId, SrnxReadError> {
    SatId::parse_name(name).ok_or_else(|| {
        SrnxReadError::BadSatelliteName(String::from_utf8_lossy(&name).into_owned())
    })
}

/// Return the observation-code slice to use for `system`.
///
/// For v2 MIXED files, all satellites share the single MIXED list.
/// For v3 files, each satellite system has its own list.
pub(crate) fn obs_codes_for_sat(header: &Header, system: SatSystem) -> &[ObsCode] {
    if let Some(codes) = header.obs_codes.for_system(SatSystem::MIXED) {
        return codes;
    }
    if let Some(codes) = header.obs_codes.for_system(system) {
        return codes;
    }
    header
        .obs_codes
        .by_system
        .first()
        .map_or(&[], |(_, v)| v.as_slice())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> ObsFile {
        read(&crate::test_data(name), None).unwrap()
    }

    #[test]
    fn timestamps_cross_month_year_and_leap_second() {
        let e = |year, month, day, hour, minute, seconds_e7| Epoch {
            year,
            month,
            day,
            hour,
            minute,
            seconds_e7,
            extra_ps: 0,
            clock_offset_e12: 0,
        };
        let u7 = EpochChunkKind::Epoc.units_per_second();
        let times = continuous_times(
            &[
                e(2016, 12, 31, 23, 59, 590_000_000),
                e(2016, 12, 31, 23, 59, 600_000_000),
                e(2017, 1, 1, 0, 0, 0),
            ],
            u7,
        );
        assert_eq!(times[1] - times[0], 10_000_000);
        assert_eq!(times[2] - times[1], 10_000_000);
        let times = continuous_times(&[e(2020, 2, 28, 0, 0, 0), e(2020, 3, 1, 0, 0, 0)], u7);
        assert_eq!(times[1] - times[0], 2 * 86400 * 10_000_000);

        // EPO2: picoseconds, including the sub-1e-7 s digits.
        let mut late = e(2016, 12, 31, 23, 59, 609_999_999);
        late.extra_ps = 99_999;
        let times = continuous_times(
            &[late, e(2017, 1, 1, 0, 0, 0)],
            EpochChunkKind::Epo2.units_per_second(),
        );
        assert_eq!(times[1] - times[0], 1);
    }

    /// Read the events of `rinex211_flag1_flag6.srnx` (two epochs, then two
    /// consecutive `EVTF` chunks: flag 1 with 14 lines and flag 6 with one)
    /// after `edit` changes its chunks, which it gets with the indexes of
    /// the two `EVTF` chunks.
    fn flag1_flag6_events(
        edit: impl for<'a> FnOnce(&mut Vec<crate::file::BodyChunk<'a>>, [usize; 2]),
    ) -> Result<Vec<(usize, EventRecord)>, SrnxReadError> {
        let bytes = crate::test_data("rinex211_flag1_flag6.srnx");
        let mut file = SrnxFile::parse(&bytes).unwrap();
        let header = rinex_header::parse(file.rhdr_payload().unwrap()).unwrap();
        let evtf_offset = file.sdir().unwrap().unwrap().evtf_offset;
        let evtf: Vec<usize> = (0..file.chunks.len())
            .filter(|&i| file.chunks[i].fourcc == *b"EVTF")
            .collect();
        assert_eq!(evtf.len(), 2);
        assert_eq!(evtf[1], evtf[0] + 1);
        edit(&mut file.chunks, [evtf[0], evtf[1]]);
        events(&file, &header, evtf_offset, 2)
    }

    /// Replace the payload of chunk `i` with `edit` applied to a copy of it.
    fn edit_payload(chunks: &mut [crate::file::BodyChunk<'_>], i: usize, edit: fn(&mut Vec<u8>)) {
        let mut payload = chunks[i].payload.to_vec();
        edit(&mut payload);
        chunks[i].payload = Box::leak(payload.into_boxed_slice());
    }

    #[test]
    fn events_from_evtf() {
        let events = flag1_flag6_events(|_, _| {}).unwrap();
        let summary: Vec<(usize, u8, usize)> = events
            .iter()
            .map(|(i, e)| (*i, e.flag, e.lines.len()))
            .collect();
        assert_eq!(summary, [(1, 1, 14), (1, 6, 1)]);
        assert_eq!(events[0].1.lines[0].trim(), "G13");
        assert_eq!(events[1].1.lines[0], "                         1.000 1");
    }

    #[test]
    fn events_need_not_be_consecutive() {
        // Move the second EVTF chunk after the chunk that followed it.
        let events = flag1_flag6_events(|chunks, [_, b]| chunks.swap(b, b + 1)).unwrap();
        let flags: Vec<u8> = events.iter().map(|(_, e)| e.flag).collect();
        assert_eq!(flags, [1, 6]);
    }

    #[test]
    fn sdir_must_point_at_the_first_evtf() {
        let err = flag1_flag6_events(|chunks, [a, _]| chunks[a].fourcc = *b"XXXX");
        assert!(matches!(
            err,
            Err(SrnxReadError::File(FileError::EvtfOffset { .. }))
        ));
    }

    #[test]
    fn bad_events_are_errors() {
        let bad = |result: Result<_, SrnxReadError>| match result {
            Err(SrnxReadError::BadEvent { index, reason }) => (index, reason),
            other => panic!("expected BadEvent, got {other:?}"),
        };
        let err = flag1_flag6_events(|chunks, [_, b]| {
            edit_payload(chunks, b, |p| p[0] = 3);
        });
        assert_eq!(bad(err), (1, "epoch index is past the last epoch"));
        let err = flag1_flag6_events(|chunks, [a, _]| edit_payload(chunks, a, |p| p.truncate(1)));
        assert_eq!(bad(err), (0, "no epoch line"));
        // Flag 6 in column 29 of " 24 10  1  0  0  2.0000000  6  1G01".
        let err = flag1_flag6_events(|chunks, [_, b]| {
            edit_payload(chunks, b, |p| {
                assert_eq!(p[1 + 28], b'6');
                p[1 + 28] = b'7';
            });
        });
        assert_eq!(bad(err), (1, "epoch flag is not 1 to 6"));

        let err = flag1_flag6_events(|chunks, [a, _]| edit_payload(chunks, a, |p| p[5] = 0xFF));
        assert!(matches!(
            err,
            Err(SrnxReadError::File(FileError::Evtf {
                source: crate::evtf::EvtfError::NotUtf8,
                ..
            }))
        ));
    }

    #[test]
    fn a_sdir_name_must_be_a_rinex_identifier() {
        // "g01" is not a satellite name: a lower-case system letter would
        // otherwise silently collide with G01.
        let bytes = one_epoch_srnx_named(*b"SOCD", &socd_payload(1, 7), &[*b"g01"]);
        assert!(matches!(
            read(&bytes, None),
            Err(SrnxReadError::BadSatelliteName(name)) if name == "g01"
        ));
    }

    #[test]
    fn a_sdir_name_alias_is_a_duplicate() {
        // " 01" and "G 1" are G01 too: listing them as well would give G01
        // two values for the epoch.
        for alias in [*b" 01", *b"G 1"] {
            let bytes = one_epoch_srnx_named(*b"SOCD", &socd_payload(1, 7), &[*b"G01", alias]);
            let g01 = SatId::new(SatSystem::GPS, 1);
            assert!(matches!(
                read(&bytes, None),
                Err(SrnxReadError::DuplicateSatellite(sat)) if sat == g01
            ));
            assert!(matches!(
                catalog(&bytes),
                Err(SrnxReadError::DuplicateSatellite(sat)) if sat == g01
            ));
        }
    }

    #[test]
    fn fixture_at01_epoch_count() {
        let obs = fixture("at01_v211_1hz.srnx");
        assert_eq!(obs.epochs.len(), 120);
    }

    #[test]
    fn fixture_at01_first_epoch() {
        let obs = fixture("at01_v211_1hz.srnx");
        let e = &obs.epochs[0];
        assert_eq!(e.year, 2020);
        assert_eq!(e.month, 7);
        assert_eq!(e.day, 18);
        assert_eq!(e.hour, 0);
        assert_eq!(e.minute, 0);
        assert_eq!(e.seconds_e7, 0);
    }

    #[test]
    fn fixture_at01_g13_l1_first_epoch() {
        let obs = fixture("at01_v211_1hz.srnx");
        let l1 = ObsCode::from_v2("L1").unwrap();
        let g13 = SatId::new(SatSystem::GPS, 13);
        let series = obs.observations.get(&(g13, l1)).expect("G13/L1 present");
        let &(ep_idx, first) = series.first().expect("at least one G13/L1 observation");
        assert_eq!(ep_idx, 0);
        assert_eq!(
            first.value, 128_250_890_432,
            "G13/L1 milli-cycles at epoch 0"
        );
    }

    #[test]
    fn rejects_tampered_file_digest() {
        let mut tampered = crate::test_data("at01_v211_1hz.srnx");
        let last = tampered.len() - 1;
        tampered[last] ^= 0x01; // flip a bit in the file-level digest
        match read(&tampered, None) {
            Ok(_) => panic!("tampered file should error"),
            Err(SrnxReadError::FileDigestMismatch) => {} // expected
            Err(e) => panic!("expected FileDigestMismatch, got {e:?}"),
        }
    }

    #[test]
    fn fixture_at01_observations_sorted() {
        let obs = fixture("at01_v211_1hz.srnx");
        let num_epochs = obs.epochs.len();
        for (key, series) in &obs.observations {
            assert!(
                series.windows(2).all(|w| w[0].0 < w[1].0),
                "series for {}/{} not sorted",
                key.0,
                key.1.as_str()
            );
            assert!(
                series.last().map_or(0, |&(i, _)| i) < num_epochs,
                "epoch index out of range for {}/{}",
                key.0,
                key.1.as_str()
            );
        }
    }

    /// Append a chunk with no digest.
    fn push_chunk(out: &mut Vec<u8>, fourcc: [u8; 4], payload: &[u8]) {
        out.extend_from_slice(&fourcc);
        crate::leb128::write_uleb128(out, payload.len() as u64);
        out.extend_from_slice(payload);
    }

    /// An SRNX file with null digests: a one-code v2 header, one epoch, one
    /// satellite named `G01`, and `payload` as its only observation chunk.
    fn one_epoch_srnx(fourcc: [u8; 4], payload: &[u8]) -> Vec<u8> {
        one_epoch_srnx_named(fourcc, payload, &[*b"G01"])
    }

    /// [`one_epoch_srnx`] with the satellite named `names[0]`, and an `SDIR`
    /// entry for each of `names` that points at its `SATE` chunk.
    fn one_epoch_srnx_named(fourcc: [u8; 4], payload: &[u8], names: &[[u8; 3]]) -> Vec<u8> {
        let name = names[0];
        let head = 21; // FOURCC + one-byte length + 16-byte padded payload
        let mut body = Vec::new();
        let rhdr = format!(
            "{:<60}{}\n{:<60}{}\n{:<60}{}\n",
            "     2.11           OBSERVATION DATA    G (GPS)",
            "RINEX VERSION / TYPE",
            "     1    L1",
            "# / TYPES OF OBSERV",
            "",
            "END OF HEADER"
        );
        push_chunk(&mut body, *b"RHDR", rhdr.as_bytes());

        let epoc_offset = body.len();
        let mut epoc = Vec::new();
        crate::leb128::write_uleb128(&mut epoc, 1); // one epoch
        crate::leb128::write_sleb128(&mut epoc, 0); // interval
        crate::leb128::write_uleb128(&mut epoc, 0); // epochs in span minus 1
        crate::leb128::write_uleb128(&mut epoc, 20_241_001); // date
        crate::leb128::write_uleb128(&mut epoc, 0); // time of day
        push_chunk(&mut body, *b"EPOC", &epoc);

        // The SATE offset is relative to its own FOURCC, so frame it once
        // with a placeholder to learn its length.
        let sate_offset = body.len();
        let sate = |relative: i64| {
            let mut p = name.to_vec();
            p.push(0);
            crate::leb128::write_sleb128(&mut p, relative);
            crate::leb128::write_uleb128(&mut p, 0); // one presence run
            crate::leb128::write_uleb128(&mut p, 0); // zero absent epochs
            crate::leb128::write_uleb128(&mut p, 0); // one present epoch
            p
        };
        let mut framed = Vec::new();
        push_chunk(&mut framed, *b"SATE", &sate(0));
        let relative = i64::try_from(framed.len()).unwrap();
        push_chunk(&mut body, *b"SATE", &sate(relative));
        push_chunk(&mut body, fourcc, payload);

        let mut sdir = Vec::new();
        crate::leb128::write_uleb128(&mut sdir, (head + epoc_offset) as u64);
        crate::leb128::write_uleb128(&mut sdir, 0); // no events
        for name in names {
            sdir.extend_from_slice(name);
            crate::leb128::write_uleb128(&mut sdir, (head + sate_offset) as u64);
        }
        let sdir_offset = head + body.len();
        push_chunk(&mut body, *b"SDIR", &sdir);

        let mut header = vec![1, 0, 0, 0]; // v1.0, null chunk and file digests
        crate::leb128::write_uleb128(&mut header, sdir_offset as u64);
        header.resize(16, 0);
        let mut out = Vec::new();
        push_chunk(&mut out, *b"SRNX", &header);
        assert_eq!(out.len(), head);
        out.extend(body);
        out
    }

    /// A `SOCD` payload with `count` SLEB128-coded values.
    fn socd_payload(count: u64, value: i64) -> Vec<u8> {
        let mut p = b"G01\0L1\0\0".to_vec();
        crate::leb128::write_uleb128(&mut p, count - 1);
        p.extend([0, 0]); // no LLI or SSI runs
        let mut obs = vec![0, 0xFF]; // schema 0, SLEB128 run
        crate::leb128::write_uleb128(&mut obs, count - 1);
        for _ in 0..count {
            crate::leb128::write_sleb128(&mut obs, value);
        }
        crate::leb128::write_uleb128(&mut p, obs.len() as u64);
        p.extend(obs);
        p
    }

    #[test]
    fn observation_count_must_match_the_satellite_epochs() {
        let key = (
            SatId::new(SatSystem::GPS, 1),
            ObsCode::from_v2("L1").unwrap(),
        );
        let file = read(&one_epoch_srnx(*b"SOCD", &socd_payload(1, 7)), None).unwrap();
        assert_eq!(
            file.observations[&key],
            [(
                0,
                Obs {
                    value: 7,
                    lli: b' ',
                    ssi: b' '
                }
            )]
        );

        // A chunk may not declare more observations than the satellite has
        // epochs, and the check must precede decoding the chunk.
        assert!(matches!(
            read(&one_epoch_srnx(*b"SOCD", &socd_payload(2, 7)), None),
            Err(SrnxReadError::File(FileError::ObservationCount {
                declared: 2,
                present: 1
            }))
        ));
    }

    #[test]
    fn chunk_digests_are_checked_without_a_file_digest() {
        // Our writer, with per-chunk digests but no file digest.
        let file = read(&one_epoch_srnx(*b"SOCD", &socd_payload(1, 7)), None).unwrap();
        let mut bytes = Vec::new();
        crate::writer::write_with(
            &file,
            &crate::writer::WriteOptions {
                obs_chunk: crate::writer::ObsChunk::Socd,
                chunk_digest: crate::writer::Digest::Crc32c,
                file_digest: crate::writer::Digest::None,
                ..Default::default()
            },
            &mut bytes,
        )
        .unwrap();
        read(&bytes, None).expect("a valid file reads");

        // The file digest covers every byte, so it is the only check made
        // when present; without one, every chunk digest must verify.
        let srnx_digest = {
            let framed = crate::chunk::read_framed(&bytes).unwrap();
            bytes.len() - framed.after_payload.len()
        };
        let mut tampered = bytes.clone();
        tampered[srnx_digest] ^= 0x01;
        assert!(matches!(
            read(&tampered, None),
            Err(SrnxReadError::SrnxChunkDigestMismatch)
        ));

        let socd = SrnxFile::parse(&bytes)
            .unwrap()
            .chunks
            .iter()
            .find(|chunk| chunk.fourcc == *b"SOCD")
            .unwrap()
            .file_offset;
        let framed = crate::chunk::read_framed(&bytes[usize::try_from(socd).unwrap()..]).unwrap();
        let digest_at = bytes.len() - framed.after_payload.len();
        let mut tampered = bytes.clone();
        tampered[digest_at] ^= 0x01;
        assert!(matches!(
            read(&tampered, None),
            Err(SrnxReadError::ChunkDigestMismatch(offset)) if offset == socd
        ));
    }
}
