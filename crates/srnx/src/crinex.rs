//! Reader for Hatanaka-compressed "Compact RINEX" (CRINEX) observation files.
//!
//! Compact RINEX is a text-to-text compression of RINEX observation files.
//! Three versions exist:
//!
//! - **1.0** (Hatanaka 1996, ION GPS-96) compresses RINEX 2.xx.
//! - **3.0** (Hatanaka 2008, Bull. GSI vol. 55) compresses RINEX 3.xx and 4.xx.
//! - **3.1** (RNXCMP 4.2.0, 2025) extends 3.0 for the RINEX 4.02 picosecond
//!   epoch field.
//!
//! [`decompress`] reconstructs the original RINEX text, byte-for-byte as
//! `CRX2RNX` (RNXCMP 4.2.0) writes it, and [`parse`] hands that text to the
//! RINEX v2 or v3/v4 reader to build an [`ObsFile`].
//!
//! ## Format summary
//!
//! The file starts with two CRINEX header lines (`CRINEX VERS   / TYPE`,
//! `CRINEX PROG / DATE`), followed by the original RINEX header verbatim.
//! Each observation epoch is then encoded as:
//!
//! 1. An **epoch line**: the RINEX epoch record (v2: the first 68 columns,
//!    without the clock; v3: the first 41 columns) with the full satellite
//!    list appended, three characters per satellite.
//! 2. A **clock line**: the receiver clock offset as an integer in units of
//!    10⁻⁹ s (v2) or 10⁻¹² s (v3), or empty when absent.  CRINEX 3.1 may
//!    append a space and the five picosecond digits of a RINEX 4.02 epoch.
//! 3. One **data line** per satellite: one field per observation type, each
//!    followed by a space separator, then the LLI and SSI characters of all
//!    types concatenated.  Each field is an integer in units of 10⁻³, or
//!    empty when the RINEX field is blank.
//!
//! Text items (the epoch line, the flags, the picosecond digits) are stored
//! as differences from the previous epoch, see [`repair`].  Numeric items
//! (clock, observations) are stored as `M`-th order differences of the
//! integer series, see [`DiffState`]; a field `M&value` (re)starts an arc
//! with maximum order `M` and the original `value`.
//!
//! An epoch line whose first column is `&` (v1) or `>` (v3) is not
//! differenced and resets all text and observation arcs.  Event records
//! (epoch flag > 1) are written verbatim, in that form, followed by their
//! special records; the next epoch is always a reset.  In CRINEX 3, a line
//! starting with `&` where an epoch line is expected is an "escape line",
//! reserved for future extensions, and is skipped.

use std::collections::HashMap;

use thiserror::Error;

use crate::model::{Catalog, ObsFile, SignalFilter};
use crate::rinex2_obs::{self, Rinex2Error};
use crate::rinex3_obs::{self, Rinex3Error};
use crate::rinex_header;

/// Highest differencing order accepted.  The 1996 paper allows any order
/// below 10; RNXCMP writes 3 and reads up to 5.
pub const MAX_DIFF_ORDER: usize = 9;

const LABEL_CRINEX_VERS: &[u8] = b"CRINEX VERS";
const LABEL_RINEX_VERSION: &[u8] = b"RINEX VERSION / TYPE";
const LABEL_END_OF_HEADER: &[u8] = b"END OF HEADER";
const LABEL_V2_OBS_TYPES: &[u8] = b"# / TYPES OF OBSERV";
const LABEL_V3_OBS_TYPES: &[u8] = b"SYS / # / OBS TYPES";

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum CrinexError {
    #[error("not a Compact RINEX file (no CRINEX VERS / TYPE record)")]
    NotCrinex,
    #[error("unsupported Compact RINEX version {0:?} (expected 1.0, 3.0 or 3.1)")]
    UnsupportedCrinexVersion(String),
    #[error("line {0}: RINEX VERSION / TYPE record missing or not version 2, 3 or 4")]
    UnsupportedRinexVersion(usize),
    #[error("no END OF HEADER line found")]
    NoEndOfHeader,
    #[error("file truncated in the middle of an epoch (after line {0})")]
    Truncated(usize),
    #[error("line {0}: malformed epoch line")]
    BadEpochLine(usize),
    #[error("line {0}: arc should be initialized, but is not")]
    NotInitialized(usize),
    #[error("line {line}: satellite system {system:?} has no observation types in the header")]
    UndefinedSystem { line: usize, system: char },
    #[error("line {0}: malformed numeric field")]
    BadField(usize),
    #[error("line {0}: differencing order exceeds {MAX_DIFF_ORDER}")]
    OrderTooHigh(usize),
    #[error("line {0}: arithmetic overflow while undoing differences")]
    Overflow(usize),
    #[error("line {line}: {what} is out of the range allowed in RINEX")]
    OutOfRange { line: usize, what: &'static str },
    #[error("decompressed RINEX v2: {0}")]
    Rinex2(#[from] Rinex2Error),
    #[error("decompressed RINEX v3/v4: {0}")]
    Rinex3(#[from] Rinex3Error),
}

/// Compact RINEX format version, from the `CRINEX VERS   / TYPE` record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrinexVersion {
    /// 1.0, for RINEX 2.xx.
    V1_0,
    /// 3.0, for RINEX 3.xx and 4.xx.
    V3_0,
    /// 3.1, for RINEX 4.02 and later (adds picosecond epoch digits).
    V3_1,
}

impl CrinexVersion {
    fn from_record(line: &[u8]) -> Result<Self, CrinexError> {
        match line.get(..3) {
            Some(b"1.0") => Ok(Self::V1_0),
            Some(b"3.0") => Ok(Self::V3_0),
            Some(b"3.1") => Ok(Self::V3_1),
            _ => Err(CrinexError::UnsupportedCrinexVersion(
                String::from_utf8_lossy(&line[..line.len().min(20)])
                    .trim()
                    .to_owned(),
            )),
        }
    }
}

/// Header label (columns 61–80) of a RINEX-style header line, if any.
fn label(line: &[u8]) -> &[u8] {
    line.get(60..).unwrap_or_default()
}

/// The first line of `bytes`, cut off after column 80: enough for a header
/// record's content and label, without scanning a file that has no line
/// breaks.
fn first_record(bytes: &[u8]) -> &[u8] {
    let head = bytes.get(..80).unwrap_or(bytes);
    let end = head.iter().position(|&b| b == b'\n').unwrap_or(head.len());
    &head[..end]
}

/// Returns the Compact RINEX version if `bytes` starts with a
/// `CRINEX VERS   / TYPE` record.
#[must_use]
pub fn sniff_version(bytes: &[u8]) -> Option<CrinexVersion> {
    let first = first_record(bytes);
    if !label(first).starts_with(LABEL_CRINEX_VERS) {
        return None;
    }
    CrinexVersion::from_record(first).ok()
}

/// Returns `true` if `bytes` looks like a Compact RINEX file of any version
/// (including versions this module does not support).
#[must_use]
pub fn is_crinex(bytes: &[u8]) -> bool {
    label(first_record(bytes)).starts_with(LABEL_CRINEX_VERS)
}

/// List the satellites and signals declared by the RINEX header inside a
/// Compact RINEX file, without decompressing its observations.
pub fn catalog(input: &[u8]) -> Result<Catalog, CrinexError> {
    if !is_crinex(input) {
        return Err(CrinexError::NotCrinex);
    }
    // The RINEX header follows the CRINEX VERS and CRINEX PROG / DATE lines.
    let mut lines = Lines::new(input);
    lines.next();
    lines.next().ok_or(CrinexError::NoEndOfHeader)?;
    let rinex = lines.rest;
    match rinex_header::sniff_version(rinex) {
        Some(v) if v.is_v3_or_later() => Ok(rinex3_obs::catalog(rinex)?),
        _ => Ok(rinex2_obs::catalog(rinex)?),
    }
}

/// Decompress a Compact RINEX file and parse the resulting RINEX text,
/// keeping only the satellite/signal pairs that `filter` selects (all of
/// them for `None`).
pub fn parse(input: &[u8], filter: Option<SignalFilter<'_>>) -> Result<ObsFile, CrinexError> {
    let text = decompress(input)?;
    match rinex_header::sniff_version(&text) {
        Some(v) if v.is_v3_or_later() => Ok(rinex3_obs::parse(&text, filter)?),
        _ => Ok(rinex2_obs::parse(&text, filter)?),
    }
}

/// Reconstruct the RINEX observation file from a Compact RINEX file.
///
/// The output uses LF line endings and has trailing blanks removed from
/// every line, as `CRX2RNX` writes it.  Two other `CRX2RNX` quirks are
/// preserved because they are part of its output format rather than of the
/// data: values below one are written without a leading zero (`-.123`),
/// and once a CRINEX 3.1 file has carried picosecond digits, every later
/// epoch with a clock offset or picoseconds repeats the last digits.
pub fn decompress(input: &[u8]) -> Result<Vec<u8>, CrinexError> {
    let mut lines = Lines::new(input);
    let first = lines.next().ok_or(CrinexError::NotCrinex)?;
    if !label(first).starts_with(LABEL_CRINEX_VERS) {
        return Err(CrinexError::NotCrinex);
    }
    let version = CrinexVersion::from_record(first)?;
    // CRINEX PROG / DATE: not part of the RINEX file.
    lines.next().ok_or(CrinexError::NoEndOfHeader)?;

    let mut out = Vec::with_capacity(input.len().saturating_mul(4));
    let mut obs_types = ObsTypeCounts::new();

    let rinex_line = lines.next().ok_or(CrinexError::NoEndOfHeader)?;
    let layout = match rinex_line.get(5) {
        Some(b'2') if label(rinex_line).starts_with(LABEL_RINEX_VERSION) => Layout::V2,
        Some(b'3' | b'4') if label(rinex_line).starts_with(LABEL_RINEX_VERSION) => Layout::V3,
        _ => return Err(CrinexError::UnsupportedRinexVersion(lines.line_no)),
    };
    push_line(&mut out, rinex_line);
    loop {
        let line = lines.next().ok_or(CrinexError::NoEndOfHeader)?;
        push_line(&mut out, line);
        obs_types.update(line);
        if label(line).starts_with(LABEL_END_OF_HEADER) {
            break;
        }
    }

    let mut decoder = Decoder {
        version,
        layout,
        lines,
        out,
        obs_types,
        epoch_line: Vec::new(),
        clock: None,
        picosec: Vec::new(),
        prev_sats: Vec::new(),
        prev_index: HashMap::new(),
    };
    decoder.run()?;
    Ok(decoder.out)
}

/// Append `line` without trailing blanks, plus a newline.
fn push_line(out: &mut Vec<u8>, line: &[u8]) {
    out.extend_from_slice(trim_end(line));
    out.push(b'\n');
}

fn trim_end(line: &[u8]) -> &[u8] {
    let end = line.iter().rposition(|&b| b != b' ').map_or(0, |i| i + 1);
    &line[..end]
}

/// Line iterator: splits on LF, strips one trailing CR, counts lines.
struct Lines<'a> {
    rest: &'a [u8],
    line_no: usize,
}

impl<'a> Lines<'a> {
    fn new(input: &'a [u8]) -> Self {
        Self {
            rest: input,
            line_no: 0,
        }
    }

    fn next(&mut self) -> Option<&'a [u8]> {
        if self.rest.is_empty() {
            return None;
        }
        let (line, rest) = match self.rest.iter().position(|&b| b == b'\n') {
            Some(i) => (&self.rest[..i], &self.rest[i + 1..]),
            None => (self.rest, &[][..]),
        };
        self.rest = rest;
        self.line_no += 1;
        Some(line.strip_suffix(b"\r").unwrap_or(line))
    }

    /// Next line where more of the current epoch is required.
    fn require(&mut self) -> Result<&'a [u8], CrinexError> {
        let line_no = self.line_no;
        self.next().ok_or(CrinexError::Truncated(line_no))
    }
}

/// Number of observation types per epoch: one count for RINEX v2, one per
/// satellite system for v3/v4.  Updated from header records, including
/// those inside event records.
struct ObsTypeCounts {
    v2: usize,
    by_system: [Option<usize>; 256],
}

impl ObsTypeCounts {
    fn new() -> Self {
        Self {
            v2: 0,
            by_system: [None; 256],
        }
    }

    fn update(&mut self, line: &[u8]) {
        let label = label(line);
        if label.starts_with(LABEL_V2_OBS_TYPES) && line.get(5).is_some_and(|&b| b != b' ') {
            self.v2 = count_field(line, 0, 6);
        } else if label.starts_with(LABEL_V3_OBS_TYPES) {
            if let Some(&sys) = line.first().filter(|&&b| b != b' ') {
                self.by_system[usize::from(sys)] = Some(count_field(line, 3, 3));
            }
        }
    }
}

/// The unsigned count in `line[start..start + width]`: leading blanks, then
/// digits.  Reading only the count's columns keeps a long digit run in the
/// rest of the line from inflating it.
fn count_field(line: &[u8], start: usize, width: usize) -> usize {
    atoi(line.get(start..start + width).unwrap_or(&[]))
}

/// C `atoi` on an unsigned count: leading blanks, then digits up to the
/// first non-digit.  Saturates instead of overflowing.
fn atoi(s: &[u8]) -> usize {
    s.iter()
        .skip_while(|&&b| b == b' ')
        .take_while(|b| b.is_ascii_digit())
        .fold(0usize, |acc, &b| {
            acc.saturating_mul(10).saturating_add(usize::from(b - b'0'))
        })
}

/// Column layout of the reconstructed RINEX file.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Layout {
    V2,
    V3,
}

impl Layout {
    /// Column of the epoch flag.
    fn flag_col(self) -> usize {
        match self {
            Self::V2 => 28,
            Self::V3 => 31,
        }
    }
    /// Column where the satellite list starts in the CRINEX epoch line.
    fn sat_list_col(self) -> usize {
        match self {
            Self::V2 => 32,
            Self::V3 => 41,
        }
    }
    /// First character of a reset (or event) epoch line in CRINEX.
    fn reset_mark(self) -> u8 {
        match self {
            Self::V2 => b'&',
            Self::V3 => b'>',
        }
    }
    /// First character of the reconstructed RINEX epoch record.
    fn record_mark(self) -> u8 {
        match self {
            Self::V2 => b' ',
            Self::V3 => b'>',
        }
    }
    /// Decimal places of the receiver clock offset (`F12.9` or `F15.12`).
    fn clock_decimals(self) -> u32 {
        match self {
            Self::V2 => 9,
            Self::V3 => 12,
        }
    }
}

/// Apply the text difference `ds` to `s`.
///
/// Each character of `ds` replaces the one at the same column of `s`,
/// except that a blank means "unchanged" and `&` means "changed to blank".
/// Where `ds` is longer than `s`, its tail is appended (with `&` mapped to
/// blank).  `s` never gets shorter: characters past the end of `ds` are
/// unchanged, since the compressor writes `&` for any it blanks out.
pub fn repair(s: &mut Vec<u8>, ds: &[u8]) {
    let common = s.len().min(ds.len());
    for (c, &d) in s[..common].iter_mut().zip(ds) {
        match d {
            b' ' => {}
            b'&' => *c = b' ',
            _ => *c = d,
        }
    }
    s.extend(
        ds[common..]
            .iter()
            .map(|&d| if d == b'&' { b' ' } else { d }),
    );
}

/// One numeric field of a CRINEX clock or data line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct DiffField {
    /// `Some(M)` for `M&value`, which starts a new arc of maximum order `M`.
    reset: Option<usize>,
    value: i64,
}

impl DiffField {
    /// Parse `[M&][-]digits`.  Returns `Ok(None)` for an empty field.
    fn parse(field: &[u8], line: usize) -> Result<Option<Self>, CrinexError> {
        if field.is_empty() {
            return Ok(None);
        }
        let (reset, digits) = match field {
            [m, b'&', value @ ..] => {
                if !m.is_ascii_digit() {
                    return Err(CrinexError::BadField(line));
                }
                let order = usize::from(m - b'0');
                if order > MAX_DIFF_ORDER {
                    return Err(CrinexError::OrderTooHigh(line));
                }
                (Some(order), value)
            }
            _ => (None, field),
        };
        let (negative, digits) = match digits {
            [b'-', magnitude @ ..] => (true, magnitude),
            _ => (false, digits),
        };
        if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
            return Err(CrinexError::BadField(line));
        }
        let mut value: i64 = 0;
        for &d in digits {
            value = value
                .checked_mul(10)
                .and_then(|v| v.checked_add(i64::from(d - b'0')))
                .ok_or(CrinexError::BadField(line))?;
        }
        Ok(Some(Self {
            reset,
            value: if negative { -value } else { value },
        }))
    }
}

/// State of one numeric arc: the current value and its differences.
///
/// `y[k]` is the `k`-th order difference at the current epoch, for
/// `k ≤ order`.  An arc starts at order 0 and gains one order per epoch
/// until it reaches `max_order`; the stored field is always `y[order]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DiffState {
    order: usize,
    max_order: usize,
    y: [i64; MAX_DIFF_ORDER + 1],
}

impl DiffState {
    /// The reconstructed value at this epoch.
    #[must_use]
    pub fn value(&self) -> i64 {
        self.y[0]
    }

    /// Apply one field to the arc's state at the previous epoch.
    fn next(prev: Option<&Self>, field: DiffField, line: usize) -> Result<Self, CrinexError> {
        let mut y = [0; MAX_DIFF_ORDER + 1];
        if let Some(max_order) = field.reset {
            y[0] = field.value;
            return Ok(Self {
                order: 0,
                max_order,
                y,
            });
        }
        let prev = prev.ok_or(CrinexError::NotInitialized(line))?;
        let order = (prev.order + 1).min(prev.max_order);
        y[order] = field.value;
        for k in (0..order).rev() {
            y[k] = y[k + 1]
                .checked_add(prev.y[k])
                .ok_or(CrinexError::Overflow(line))?;
        }
        Ok(Self {
            order,
            max_order: prev.max_order,
            y,
        })
    }
}

/// Per-satellite state carried to the next epoch.
struct SatState {
    id: [u8; 3],
    /// One arc per observation type; `None` where the field was blank.
    obs: Vec<Option<DiffState>>,
    /// LLI and SSI characters, two per observation type.
    flags: Vec<u8>,
}

struct Decoder<'a> {
    version: CrinexVersion,
    layout: Layout,
    lines: Lines<'a>,
    out: Vec<u8>,
    obs_types: ObsTypeCounts,
    /// Reconstructed epoch line (without trailing blanks).
    epoch_line: Vec<u8>,
    clock: Option<DiffState>,
    /// Picosecond digits of the epoch (CRINEX 3.1); never reset.
    picosec: Vec<u8>,
    prev_sats: Vec<SatState>,
    /// Satellite ID → index into `prev_sats` (first occurrence).
    prev_index: HashMap<[u8; 3], usize>,
}

impl<'a> Decoder<'a> {
    fn run(&mut self) -> Result<(), CrinexError> {
        // After an event, the next epoch line must be a reset.
        let mut need_reset = false;
        while let Some(line) = self.next_epoch_line() {
            if line.first() == Some(&0x1A) {
                break; // DOS end-of-file marker
            }
            let line_no = self.lines.line_no;
            if line.first() == Some(&self.layout.reset_mark()) {
                let flag = line.get(self.layout.flag_col()).copied();
                if !matches!(flag, Some(b'0' | b'1')) {
                    self.event(line)?;
                    need_reset = true;
                    continue;
                }
                self.epoch_line.clear();
                self.prev_sats.clear();
                self.prev_index.clear();
                let mut line = line.to_vec();
                line[0] = self.layout.record_mark();
                repair(&mut self.epoch_line, &line);
            } else if need_reset {
                return Err(CrinexError::NotInitialized(line_no));
            } else {
                repair(&mut self.epoch_line, line);
            }
            need_reset = false;
            let len = trim_end(&self.epoch_line).len();
            self.epoch_line.truncate(len);
            self.epoch(line_no)?;
        }
        Ok(())
    }

    /// Next line where an epoch may start, skipping CRINEX 3 escape lines.
    fn next_epoch_line(&mut self) -> Option<&'a [u8]> {
        loop {
            let line = self.lines.next()?;
            if self.version == CrinexVersion::V1_0 || line.first() != Some(&b'&') {
                return Some(line);
            }
        }
    }

    /// Copy an event record (epoch flag 2–6) and its special records.
    fn event(&mut self, line: &[u8]) -> Result<(), CrinexError> {
        let mut line = trim_end(line).to_vec();
        line[0] = self.layout.record_mark();
        push_line(&mut self.out, &line);
        if line.len() > 29 {
            let count = atoi(line.get(self.layout.flag_col() + 1..).unwrap_or_default());
            for _ in 0..count {
                let record = self.lines.require()?;
                push_line(&mut self.out, record);
                self.obs_types.update(record);
            }
        }
        Ok(())
    }

    /// Decode one observation epoch (flag 0 or 1) whose epoch line has
    /// already been reconstructed.
    fn epoch(&mut self, line_no: usize) -> Result<(), CrinexError> {
        let layout = self.layout;
        let flag_col = layout.flag_col();
        let el = &self.epoch_line;
        if el.first() != Some(&layout.record_mark())
            || el.len() <= flag_col
            || el[flag_col - 2] != b' '
            || el[flag_col - 1] != b' '
            || !el[flag_col].is_ascii_digit()
        {
            return Err(CrinexError::BadEpochLine(line_no));
        }
        let num_sats = atoi(el.get(flag_col + 1..flag_col + 4).unwrap_or_default());

        // Clock line (with CRINEX 3.1 picosecond digits).
        let clock_line = self.lines.require()?;
        let clock_line_no = self.lines.line_no;
        let (clock_text, picosec_text) = if self.version == CrinexVersion::V3_1 {
            match clock_line.iter().position(|&b| b == b' ') {
                Some(i) => (&clock_line[..i], &clock_line[i + 1..]),
                None => (clock_line, &[][..]),
            }
        } else {
            (clock_line, &[][..])
        };
        self.clock = match DiffField::parse(clock_text, clock_line_no)? {
            Some(field) => Some(DiffState::next(self.clock.as_ref(), field, clock_line_no)?),
            None => None,
        };
        if picosec_text.len() > 5 {
            return Err(CrinexError::OutOfRange {
                line: clock_line_no,
                what: "picosecond field",
            });
        }
        repair(&mut self.picosec, picosec_text);

        self.write_epoch_record(num_sats, clock_line_no)?;

        let sat_list_col = layout.sat_list_col();
        let mut sats = Vec::with_capacity(num_sats);
        for i in 0..num_sats {
            let col = sat_list_col + 3 * i;
            let mut id = [b' '; 3];
            for (k, b) in id.iter_mut().enumerate() {
                if let Some(&c) = self.epoch_line.get(col + k) {
                    *b = c;
                }
            }
            sats.push(self.satellite(id)?);
        }

        self.prev_index.clear();
        for (i, sat) in sats.iter().enumerate() {
            self.prev_index.entry(sat.id).or_insert(i);
        }
        self.prev_sats = sats;
        Ok(())
    }

    /// Write the RINEX epoch record, including the clock offset, the
    /// picosecond digits (RINEX 4.02), and for v2 the continuation lines of
    /// the satellite list.
    fn write_epoch_record(&mut self, num_sats: usize, line_no: usize) -> Result<(), CrinexError> {
        let clock = match &self.clock {
            Some(state) => Some(format_clock(
                state.value(),
                self.layout.clock_decimals(),
                line_no,
            )?),
            None => None,
        };
        let el = &self.epoch_line;
        let out = &mut self.out;
        match self.layout {
            Layout::V2 => {
                let head = &el[..el.len().min(68)];
                out.extend_from_slice(head);
                if let Some(clock) = clock {
                    out.resize(out.len() + 68 - head.len(), b' ');
                    out.extend_from_slice(clock.as_bytes());
                }
                out.push(b'\n');
                // Satellites 13, 25, ... start continuation lines.
                for k in 0..num_sats.saturating_sub(1) / 12 {
                    let sats = el.get(68 + 36 * k..).unwrap_or_default();
                    out.extend_from_slice(&[b' '; 32]);
                    out.extend_from_slice(&sats[..sats.len().min(36)]);
                    out.push(b'\n');
                }
            }
            Layout::V3 => {
                let head = &el[..el.len().min(41)];
                // CRINEX 3.1 only: CRX2RNX keeps the digits from earlier epochs.
                let picosec = &self.picosec[..self.picosec.len().min(5)];
                if clock.is_none() && picosec.is_empty() {
                    out.extend_from_slice(trim_end(head));
                } else {
                    out.extend_from_slice(head);
                    out.resize(out.len() + 41 - head.len(), b' ');
                    match clock {
                        Some(clock) => out.extend_from_slice(clock.as_bytes()),
                        None => out.extend_from_slice(&[b' '; 15]),
                    }
                    if !picosec.is_empty() {
                        out.push(b' ');
                        out.resize(out.len() + 5 - picosec.len(), b' ');
                        out.extend_from_slice(picosec);
                    }
                }
                out.push(b'\n');
            }
        }
        Ok(())
    }

    /// Decode and write the observation record of one satellite.
    fn satellite(&mut self, id: [u8; 3]) -> Result<SatState, CrinexError> {
        let num_types = match self.layout {
            Layout::V2 => self.obs_types.v2,
            Layout::V3 => self.obs_types.by_system[usize::from(id[0])].ok_or(
                CrinexError::UndefinedSystem {
                    line: self.lines.line_no,
                    system: char::from(id[0]),
                },
            )?,
        };
        let line = self.lines.require()?;
        let line_no = self.lines.line_no;
        let prev = self.prev_index.get(&id).map(|&i| &self.prev_sats[i]);

        // Split the first `num_types` blank-separated fields; the rest of
        // the line is the flag text.  Fields past the end are empty.
        let mut obs = Vec::with_capacity(num_types);
        let mut rest = line;
        for j in 0..num_types {
            let end = rest.iter().position(|&b| b == b' ').unwrap_or(rest.len());
            let field = DiffField::parse(&rest[..end], line_no)?;
            rest = rest.get(end + 1..).unwrap_or_default();
            obs.push(match field {
                None => None,
                Some(field) => {
                    let prev_arc = prev.and_then(|p| p.obs.get(j).copied().flatten());
                    Some(DiffState::next(prev_arc.as_ref(), field, line_no)?)
                }
            });
        }

        let mut flags = match prev {
            Some(p) => p.flags.clone(),
            None => Vec::new(),
        };
        repair(&mut flags, rest);
        flags.resize(2 * num_types, b' ');

        // Write the RINEX record: v3 has the satellite ID and all fields on
        // one line; v2 wraps after every 5 fields.
        let out = &mut self.out;
        let mut line_start = out.len();
        if self.layout == Layout::V3 {
            out.extend_from_slice(&id);
        }
        for (j, arc) in obs.iter().enumerate() {
            let pair = &mut flags[2 * j..2 * j + 2];
            match arc {
                Some(arc) => {
                    out.extend_from_slice(format_obs(arc.value(), line_no)?.as_bytes());
                    out.extend_from_slice(pair);
                }
                // CRINEX 1 drops the flags of blank fields; CRINEX 3 keeps
                // them independently of the data.
                None if self.version == CrinexVersion::V1_0 => {
                    out.extend_from_slice(&[b' '; 16]);
                    pair.fill(b' ');
                }
                None => {
                    out.extend_from_slice(&[b' '; 14]);
                    out.extend_from_slice(pair);
                }
            }
            if j + 1 == num_types || (self.layout == Layout::V2 && (j + 1) % 5 == 0) {
                let len = line_start + trim_end(&out[line_start..]).len();
                out.truncate(len);
                out.push(b'\n');
                line_start = out.len();
            }
        }
        if num_types == 0 && self.layout == Layout::V3 {
            out.push(b'\n');
        }
        Ok(SatState { id, obs, flags })
    }
}

/// Format an observation in units of 10⁻³ as RINEX `F14.3`, the way
/// `CRX2RNX` does: no leading zero before the decimal point.
fn format_obs(value: i64, line: usize) -> Result<String, CrinexError> {
    let s = format_fixed(value, 3);
    if s.len() > 14 {
        return Err(CrinexError::OutOfRange {
            line,
            what: "observation",
        });
    }
    Ok(format!("{s:>14}"))
}

/// Format a receiver clock offset in units of 10⁻ᵈ s as RINEX `F12.9`
/// (`d` = 9) or `F15.12` (`d` = 12), without a leading zero.
fn format_clock(value: i64, decimals: u32, line: usize) -> Result<String, CrinexError> {
    let s = format_fixed(value, decimals);
    let width = decimals as usize + 3;
    if s.len() > width {
        return Err(CrinexError::OutOfRange {
            line,
            what: "receiver clock offset",
        });
    }
    Ok(format!("{s:>width$}"))
}

/// `value / 10^decimals` in fixed point, omitting a zero integer part.
fn format_fixed(value: i64, decimals: u32) -> String {
    let unit = 10_u64.pow(decimals);
    let magnitude = value.unsigned_abs();
    let (int, frac) = (magnitude / unit, magnitude % unit);
    let sign = if value < 0 { "-" } else { "" };
    let width = decimals as usize;
    if int == 0 {
        format!("{sign}.{frac:0width$}")
    } else {
        format!("{sign}{int}.{frac:0width$}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn repaired(old: &str, diff: &str) -> String {
        let mut s = old.as_bytes().to_vec();
        repair(&mut s, diff.as_bytes());
        String::from_utf8(s).unwrap()
    }

    #[test]
    fn repair_text_differences() {
        assert_eq!(repaired("", "AB C"), "AB C");
        assert_eq!(repaired("", "&7&&"), " 7  ");
        assert_eq!(repaired("ABCD", "  X"), "ABXD");
        assert_eq!(repaired("ABCD", " &"), "A CD");
        assert_eq!(repaired("AB", "  CD&E"), "ABCD E");
        assert_eq!(repaired("ABC", ""), "ABC");
    }

    fn field(s: &str) -> Option<DiffField> {
        DiffField::parse(s.as_bytes(), 1).unwrap()
    }

    #[test]
    fn diff_field_parsing() {
        assert_eq!(field(""), None);
        assert_eq!(
            field("3&-1234567"),
            Some(DiffField {
                reset: Some(3),
                value: -1_234_567
            })
        );
        assert_eq!(
            field("0"),
            Some(DiffField {
                reset: None,
                value: 0
            })
        );
        for bad in ["-", "3&", "+5", "1 2", "x&5", "12a", "99999999999999999999"] {
            assert!(DiffField::parse(bad.as_bytes(), 1).is_err(), "{bad:?}");
        }
        assert!(matches!(
            DiffField::parse(b"9&1", 1),
            Ok(Some(DiffField { reset: Some(9), .. }))
        ));
    }

    /// Third-order differencing reproduces a polynomial series, as in the
    /// Appendix 1 algorithm of Hatanaka (2008).
    #[test]
    fn diff_state_undoes_third_order_differences() {
        let series: Vec<i64> = (0..10_i64)
            .map(|t| 7 + 3 * t - 2 * t * t + t * t * t)
            .collect();
        // Forward differencing, ramping the order up to 3.
        let mut fields = vec![DiffField {
            reset: Some(3),
            value: series[0],
        }];
        let mut diffs = vec![series[0]];
        for (i, &v) in series.iter().enumerate().skip(1) {
            let order = i.min(3);
            let mut d = vec![v];
            for k in 0..order {
                d.push(d[k] - diffs[k]);
            }
            fields.push(DiffField {
                reset: None,
                value: d[order],
            });
            diffs = d;
        }
        let mut state: Option<DiffState> = None;
        for (f, &want) in fields.iter().zip(&series) {
            let next = DiffState::next(state.as_ref(), *f, 1).unwrap();
            assert_eq!(next.value(), want);
            state = Some(next);
        }
        assert_eq!(state.unwrap().order, 3);
    }

    #[test]
    fn diff_state_requires_initialization() {
        let f = DiffField {
            reset: None,
            value: 5,
        };
        assert!(matches!(
            DiffState::next(None, f, 7),
            Err(CrinexError::NotInitialized(7))
        ));
    }

    #[test]
    fn obs_formatting_matches_crx2rnx() {
        let f = |v| format_obs(v, 1).unwrap();
        assert_eq!(f(23_629_347_915), "  23629347.915");
        assert_eq!(f(-1_234_567), "     -1234.567");
        assert_eq!(f(12_345), "        12.345");
        assert_eq!(f(-2_345), "        -2.345");
        assert_eq!(f(123), "          .123");
        assert_eq!(f(-5), "         -.005");
        assert_eq!(f(0), "          .000");
        assert_eq!(f(9_999_999_999_999), "9999999999.999");
        assert_eq!(f(-999_999_999_999), "-999999999.999");
        assert!(format_obs(10_000_000_000_000, 1).is_err());
        assert!(format_obs(-1_000_000_000_000, 1).is_err());
    }

    #[test]
    fn clock_formatting_matches_crx2rnx() {
        let f = |v, d| format_clock(v, d, 1).unwrap();
        assert_eq!(f(-123_456_789, 9), " -.123456789");
        assert_eq!(f(500_000_000, 9), "  .500000000");
        assert_eq!(f(12_500_000_000, 9), "12.500000000");
        assert_eq!(f(0, 9), "  .000000000");
        assert_eq!(f(-5, 9), " -.000000005");
        assert_eq!(f(-9_999_999_999, 9), "-9.999999999");
        assert!(format_clock(-10_000_000_000, 9, 1).is_err());
        assert!(format_clock(100_000_000_000, 9, 1).is_err());
        assert_eq!(f(100_001, 12), "  .000000100001");
        assert_eq!(f(-1, 12), " -.000000000001");
    }

    #[test]
    fn sniffing() {
        let v1 =
            b"1.0                 COMPACT RINEX FORMAT                    CRINEX VERS   / TYPE\n";
        let v31 =
            b"3.1                 COMPACT RINEX FORMAT                    CRINEX VERS   / TYPE\r\n";
        let v9 =
            b"9.9                 COMPACT RINEX FORMAT                    CRINEX VERS   / TYPE\n";
        assert_eq!(sniff_version(v1), Some(CrinexVersion::V1_0));
        assert_eq!(sniff_version(v31), Some(CrinexVersion::V3_1));
        assert_eq!(sniff_version(v9), None);
        assert!(is_crinex(v9));
        // Only the first 80 bytes are examined.
        let mut long = v1[..80].to_vec();
        long.resize(1 << 20, b'x');
        assert_eq!(sniff_version(&long), Some(CrinexVersion::V1_0));
        assert!(!is_crinex(&vec![b'x'; 1 << 20]));
        assert!(!is_crinex(b""));
        assert!(!is_crinex(&v1[1..]));
        assert!(!is_crinex(
            b"     2.11           OBSERVATION DATA    G (GPS)             RINEX VERSION / TYPE\n"
        ));
        assert!(matches!(
            decompress(v9),
            Err(CrinexError::UnsupportedCrinexVersion(v)) if v == "9.9"
        ));
    }

    const V1_HEADER: &str = "\
1.0                 COMPACT RINEX FORMAT                    CRINEX VERS   / TYPE
RNX2CRX ver.4.2.0                       14-Oct-26 00:52     CRINEX PROG / DATE
     2.11           OBSERVATION DATA    G (GPS)             RINEX VERSION / TYPE
     2    C1    L1                                          # / TYPES OF OBSERV
                                                            END OF HEADER
";

    /// A CRINEX 1.0 stream exercising a reset epoch, a differenced epoch
    /// with a new satellite and a blank field, an event record, and an
    /// epoch without a clock offset.  Expected output is from `CRX2RNX`.
    #[test]
    fn decompress_v1_small() {
        let crx = format!(
            "{V1_HEADER}\
&24 10  1  0  0  0.0000000  0  2G01G02
3&123456789
3&20001000000 3&100005000000  6 6
3&20002000000 3&-5  7
                 1             3      G03
-3456789
1000 5000
1000
3&20003000000 3&100015000000    7
&24 10  1  0  0  2.0000000  4  1
COMMENT LINE                                                COMMENT
&24 10  1  0  0  3.0000000  0  1G01

3&20001003000
"
        );
        let rnx = String::from_utf8(decompress(crx.as_bytes()).unwrap()).unwrap();
        let want =
            "     2.11           OBSERVATION DATA    G (GPS)             RINEX VERSION / TYPE
     2    C1    L1                                          # / TYPES OF OBSERV
                                                            END OF HEADER
 24 10  1  0  0  0.0000000  0  2G01G02                                .123456789
  20001000.000 6 100005000.000 6
  20002000.000 7         -.005
 24 10  1  0  0  1.0000000  0  3G01G02G03                             .120000000
  20001001.000 6 100005005.000 6
  20002001.000 7
  20003000.000   100015000.000 7
 24 10  1  0  0  2.0000000  4  1
COMMENT LINE                                                COMMENT
 24 10  1  0  0  3.0000000  0  1G01
  20001003.000
";
        assert_eq!(rnx, want);
    }

    #[test]
    fn v1_reset_required_after_event() {
        let crx = format!(
            "{V1_HEADER}\
&24 10  1  0  0  2.0000000  4  0
                 3        1G01
"
        );
        assert!(matches!(
            decompress(crx.as_bytes()),
            Err(CrinexError::NotInitialized(7))
        ));
    }

    #[test]
    fn truncated_epoch_is_an_error() {
        let crx = format!("{V1_HEADER}&24 10  1  0  0  0.0000000  0  2G01G02\n\n");
        assert!(matches!(
            decompress(crx.as_bytes()),
            Err(CrinexError::Truncated(7))
        ));
    }

    const V3_HEADER: &str = "\
3.1                 COMPACT RINEX FORMAT                    CRINEX VERS   / TYPE
RNX2CRX ver.4.2.0                       14-Oct-26 00:52     CRINEX PROG / DATE
     4.02           OBSERVATION DATA    M                   RINEX VERSION / TYPE
G    2 C1C L1C                                              SYS / # / OBS TYPES
E    1 C1C                                                  SYS / # / OBS TYPES
                                                            END OF HEADER
";

    /// CRINEX 3.1: picosecond digits, an escape line, flags without data,
    /// and an epoch without a clock offset.  Expected output is from
    /// `CRX2RNX`, which repeats the picosecond digits once seen.
    #[test]
    fn decompress_v31_small() {
        let crx = format!(
            "{V3_HEADER}\
> 2024 10 01 00 00  0.0000001  0  2      G05E11
3&100001 00001
3&23629347915 3&124173212411 &&&7
3&25180322100 &&
&escape line, skipped
                    1
 99999
1
0
                    2             1      E11&&&
3&-5
3&25180322500  5
"
        );
        let rnx = String::from_utf8(decompress(crx.as_bytes()).unwrap()).unwrap();
        let body = rnx.split("END OF HEADER\n").nth(1).unwrap();
        let want = "\
> 2024 10 01 00 00  0.0000001  0  2        .000000100001 00001
G05  23629347.915   124173212.411 7
E11  25180322.100
> 2024 10 01 00 00  1.0000001  0  2                      99999
G05  23629347.916                 7
E11  25180322.100
> 2024 10 01 00 00  2.0000001  0  1       -.000000000005 99999
E11  25180322.500 5
";
        assert_eq!(body, want);
    }

    #[test]
    fn a_huge_count_field_is_bounded_by_its_columns() {
        // The count field is three columns; digits after it must not inflate
        // it to an absurd value (which would be a huge allocation).
        let mut counts = ObsTypeCounts::new();
        let line = format!(
            "{:<60}SYS / # / OBS TYPES",
            "G  99999999999999999999999 C1C"
        );
        counts.update(line.as_bytes());
        assert_eq!(counts.by_system[usize::from(b'G')], Some(999));

        // The same header still decodes; the data line supplies one field.
        let header = "\
3.1                 COMPACT RINEX FORMAT                    CRINEX VERS   / TYPE
RNX2CRX ver.4.2.0                       14-Oct-26 00:52     CRINEX PROG / DATE
     4.02           OBSERVATION DATA    M                   RINEX VERSION / TYPE
G  99999999999999999999999 C1C                              SYS / # / OBS TYPES
                                                            END OF HEADER
";
        let crx = format!("{header}> 2024 10 01 00 00  0.0000000  0  1      G01\n\n3&1\n");
        let rnx = String::from_utf8(decompress(crx.as_bytes()).unwrap()).unwrap();
        let line = rnx.lines().last().unwrap();
        assert!(line.starts_with("G01"), "{line:?}");
        assert_eq!(line.matches(".001").count(), 1, "{line:?}");
    }

    #[test]
    fn v3_undefined_system_is_an_error() {
        let crx = format!("{V3_HEADER}> 2024 10 01 00 00  0.0000001  0  1      C05\n\n3&1\n");
        assert!(matches!(
            decompress(crx.as_bytes()),
            Err(CrinexError::UndefinedSystem {
                line: 8,
                system: 'C'
            })
        ));
    }
}
