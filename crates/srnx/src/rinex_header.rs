//! RINEX observation file header: raw text with query helpers.
//!
//! RINEX headers are line-oriented, fixed-column text:
//!
//! ```text
//!     2.11           OBSERVATION DATA    M (MIXED)           RINEX VERSION / TYPE
//!     ^^^^---- content (cols 1-60) ----^^^^---- label (cols 61-80) ----^^^^
//! ```
//!
//! The primary public type is [`HeaderText`], which wraps the raw text and
//! provides label-based lookup helpers.  A full semantic parser exists
//! internally for the library's own needs (extracting observation codes for
//! SOCD parsing) but is `#[doc(hidden)]` — end-user code should work with
//! the raw text directly.
//!
//! In SRNX `RHDR` chunks, newlines are LF only (per `Succinct_RINEX.md`),
//! even when the source file used CRLF. We're tolerant of either.

use std::collections::BTreeMap;

use thiserror::Error;

use crate::model::SatId;
use crate::rinex3_obs::parse_sat_id;

// ---------------------------------------------------------------------------
// Public API: HeaderText — raw text with label-based query helpers
// ---------------------------------------------------------------------------

/// Split a RINEX header line into (content, label).  Content is cols 1–60,
/// label is cols 61–80 (trailing whitespace trimmed).  Returns `None` for
/// lines shorter than 60 characters.
fn split_label(line: &str) -> Option<(&str, &str)> {
    if line.len() < 60 || !line.is_char_boundary(60) {
        return None;
    }
    let (content, label_field) = line.split_at(60);
    let label = label_field.trim_end();
    Some((content, label))
}

/// RINEX header stored as raw text, with helpers for label-based lookups.
///
/// The text is normalized to LF newlines.  Labels are extracted from
/// columns 61–80 of each line; trailing whitespace in the label is ignored
/// for matching.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HeaderText {
    text: String,
}

impl HeaderText {
    /// Create from raw bytes.  Normalizes CRLF/CR to LF, and maps any byte
    /// that is not valid UTF-8 to `?`, one per byte, so columns stay put.
    #[must_use]
    pub fn from_bytes(bytes: &[u8]) -> Self {
        Self {
            text: to_text(bytes).replace("\r\n", "\n").replace('\r', "\n"),
        }
    }

    /// Create from a `&str`. Normalizes CRLF/CR to LF.
    #[must_use]
    pub fn from_text(text: &str) -> Self {
        Self {
            text: text.replace("\r\n", "\n").replace('\r', "\n"),
        }
    }

    /// Return the raw header text (LF-normalized).
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.text
    }

    /// Return an iterator over the header lines (excluding newline chars).
    pub fn lines(&self) -> impl Iterator<Item = &str> + '_ {
        self.text.lines()
    }

    /// Return the content (cols 1–60, trimmed) of the last line whose label
    /// (cols 61–80) matches `label`.  Trailing whitespace in both the stored
    /// label and the query is ignored.
    ///
    /// Returns `None` if no matching line is found.
    #[must_use]
    pub fn get_line(&self, label: &str) -> Option<String> {
        let query_label = label.trim_end();
        let mut result = None;
        for line in self.text.lines() {
            if let Some((_, stored_label)) = split_label(line) {
                if stored_label == query_label {
                    result = Some(line[..60.min(line.len())].trim_end().to_string());
                }
            }
        }
        result
    }

    /// Return an iterator over the content (cols 1–60, trimmed) of all lines
    /// whose label matches `label`, in file order.
    ///
    /// Trailing whitespace in both the stored label and the query is ignored.
    pub fn get_lines<'a>(&'a self, label: &'a str) -> impl Iterator<Item = String> + 'a {
        let query_label = label.trim_end();
        self.text.lines().filter_map(move |line| {
            let (_, stored_label) = split_label(line)?;
            if stored_label == query_label {
                Some(line[..60.min(line.len())].trim_end().to_string())
            } else {
                None
            }
        })
    }
}

impl std::ops::Deref for HeaderText {
    type Target = str;

    fn deref(&self) -> &Self::Target {
        &self.text
    }
}

// ---------------------------------------------------------------------------
// Internal parser (#[doc(hidden)]) — used by the library for obs-code
// extraction; end-user code should prefer HeaderText directly.
// ---------------------------------------------------------------------------

#[doc(hidden)]
#[derive(Debug, Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum HeaderError {
    #[error("missing END OF HEADER")]
    NoEndOfHeader,
    #[error("missing RINEX VERSION / TYPE line")]
    MissingVersion,
    #[error("malformed RINEX VERSION / TYPE: {0}")]
    BadVersion(String),
    #[error("unsupported RINEX major version {0} (this reader handles 2 and 3/4)")]
    UnsupportedVersion(u8),
    #[error("file type {got:?} is not 'O' (observation)")]
    NotObservationFile { got: char },
    #[error("malformed {label}: {detail}")]
    Malformed { label: &'static str, detail: String },
    #[error(
        "obs-codes count mismatch for system {system:?}: declared {expected}, parsed {actual}"
    )]
    ObsCountMismatch {
        system: Option<u8>,
        expected: usize,
        actual: usize,
    },
    #[error("obs-codes continuation line with no preceding declaration")]
    OrphanObsContinuation,
    #[error("obs code {code} declared twice for system {system:?}")]
    DuplicateObsCode { system: Option<u8>, code: String },
}

/// A RINEX format version.  `minor` is the number formed by the two digits
/// after the point, so version 2.11 is `{ major: 2, minor: 11 }` and 3.04 is
/// `{ major: 3, minor: 4 }`.  Versions order numerically.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash, Ord, PartialOrd)]
pub struct Version {
    pub major: u8,
    pub minor: u8,
}

impl Version {
    #[must_use]
    pub fn is_v2(self) -> bool {
        self.major == 2
    }
    #[must_use]
    pub fn is_v3_or_later(self) -> bool {
        self.major >= 3
    }
}

/// Single-character RINEX satellite-system code (`G`, `R`, `E`, `C`, `S`,
/// `J`, `I`, or `M` for the file-level "mixed" indicator).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Ord, PartialOrd)]
pub struct SatSystem(pub u8);

impl SatSystem {
    pub const GPS: Self = Self(b'G');
    pub const GLONASS: Self = Self(b'R');
    pub const GALILEO: Self = Self(b'E');
    pub const BEIDOU: Self = Self(b'C');
    pub const SBAS: Self = Self(b'S');
    pub const QZSS: Self = Self(b'J');
    pub const IRNSS: Self = Self(b'I');
    pub const MIXED: Self = Self(b'M');

    #[must_use]
    pub fn as_char(self) -> char {
        self.0 as char
    }
}

/// Observation code. v2 codes are 2 chars (e.g. `L1`); v3+ codes are 3 chars
/// (e.g. `L1C`). Stored as a fixed 3-byte array; the trailing byte is `0`
/// for v2 codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Ord, PartialOrd)]
pub struct ObsCode {
    bytes: [u8; 3],
}

impl ObsCode {
    /// Construct from a 2-char v2 code. Trims whitespace; rejects empty.
    #[must_use]
    pub fn from_v2(s: &str) -> Option<Self> {
        let t = s.trim();
        let b = t.as_bytes();
        if b.is_empty() || b.len() > 2 {
            return None;
        }
        let mut out = [0u8; 3];
        out[..b.len()].copy_from_slice(b);
        Some(Self { bytes: out })
    }

    /// Construct from a 3-char v3+ code. Trims whitespace; rejects empty.
    #[must_use]
    pub fn from_v3(s: &str) -> Option<Self> {
        let t = s.trim();
        let b = t.as_bytes();
        if b.is_empty() || b.len() > 3 {
            return None;
        }
        let mut out = [0u8; 3];
        out[..b.len()].copy_from_slice(b);
        Some(Self { bytes: out })
    }

    #[must_use]
    pub fn as_str(&self) -> &str {
        let end = self.bytes.iter().position(|&b| b == 0).unwrap_or(3);
        // All inputs are ASCII by construction.
        std::str::from_utf8(&self.bytes[..end]).unwrap_or("")
    }
}

/// Per-system observation-code lists.
///
/// In RINEX v2, all satellites share a single list; in v3+, each satellite
/// system declares its own list. We model both with `Vec<(SatSystem, Vec<…>)>`,
/// using the synthetic `SatSystem::MIXED` (`'M'`) entry for v2.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ObsCodes {
    /// In file order. For v2 this contains exactly one entry whose system
    /// is `SatSystem::MIXED`.
    pub by_system: Vec<(SatSystem, Vec<ObsCode>)>,
}

impl ObsCodes {
    #[must_use]
    pub fn for_system(&self, sys: SatSystem) -> Option<&[ObsCode]> {
        self.by_system
            .iter()
            .find(|(s, _)| *s == sys)
            .map(|(_, v)| v.as_slice())
    }

    #[must_use]
    pub fn total_codes(&self) -> usize {
        self.by_system.iter().map(|(_, v)| v.len()).sum()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Receiver {
    pub serial: String,
    pub model: String,
    pub version: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Antenna {
    pub serial: String,
    pub model: String,
}

/// Calendar timestamp as written in RINEX `TIME OF FIRST OBS` /
/// `TIME OF LAST OBS`. Sub-second precision is preserved as
/// `seconds_e7` (seconds × 10⁷), matching the granularity used by the SRNX
/// `EPOC` chunk.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EpochTime {
    pub year: i32,
    pub month: u8,
    pub day: u8,
    pub hour: u8,
    pub minute: u8,
    /// Whole seconds and ten-millionths of a second packed as `s × 10⁷`.
    pub seconds_e7: u64,
    /// E.g. `"GPS"`, `"GLO"`, `"GAL"`, `"BDT"`, `"QZS"`, or empty.
    pub time_system: String,
}

#[doc(hidden)]
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Header {
    pub version: Version,
    /// File type, e.g. `b'O'` for observation.
    pub file_type: u8,
    /// Satellite-system indicator from the version line; `b'M'` for mixed.
    pub sat_system: SatSystem,
    pub obs_codes: ObsCodes,

    pub program: Option<String>,
    pub run_by: Option<String>,
    pub date_str: Option<String>,
    pub marker_name: Option<String>,
    pub marker_number: Option<String>,
    pub marker_type: Option<String>,
    pub observer: Option<String>,
    pub agency: Option<String>,
    pub receiver: Option<Receiver>,
    pub antenna: Option<Antenna>,
    pub approx_position: Option<[f64; 3]>,
    pub antenna_delta_hen: Option<[f64; 3]>,
    pub interval: Option<f64>,
    pub time_of_first_obs: Option<EpochTime>,
    pub time_of_last_obs: Option<EpochTime>,
    pub leap_seconds: Option<i32>,
    pub wavelength_fact_l1l2: Option<(u8, u8)>,
    /// Satellites listed by `PRN / # OF OBS` records, in header order and
    /// without repeats; empty if the header has none.
    pub satellites: Vec<SatId>,

    pub comments: Vec<String>,
    /// Recognized labels with content we didn't promote to a typed field,
    /// plus any unknown labels. Most-recent-wins per label.
    pub unknown_lines: BTreeMap<String, String>,
}

impl Default for SatSystem {
    fn default() -> Self {
        Self(b' ')
    }
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

const LABEL_VERSION: &str = "RINEX VERSION / TYPE";
const LABEL_PGM: &str = "PGM / RUN BY / DATE";
const LABEL_COMMENT: &str = "COMMENT";
const LABEL_MARKER_NAME: &str = "MARKER NAME";
const LABEL_MARKER_NUMBER: &str = "MARKER NUMBER";
const LABEL_MARKER_TYPE: &str = "MARKER TYPE";
const LABEL_OBSERVER: &str = "OBSERVER / AGENCY";
const LABEL_RECEIVER: &str = "REC # / TYPE / VERS";
const LABEL_ANTENNA: &str = "ANT # / TYPE";
const LABEL_APPROX_POS: &str = "APPROX POSITION XYZ";
const LABEL_ANT_DELTA: &str = "ANTENNA: DELTA H/E/N";
const LABEL_OBS_TYPES_V2: &str = "# / TYPES OF OBSERV";
const LABEL_OBS_TYPES_V3: &str = "SYS / # / OBS TYPES";
const LABEL_INTERVAL: &str = "INTERVAL";
const LABEL_TIME_FIRST: &str = "TIME OF FIRST OBS";
const LABEL_TIME_LAST: &str = "TIME OF LAST OBS";
const LABEL_LEAP: &str = "LEAP SECONDS";
const LABEL_WAVELEN: &str = "WAVELENGTH FACT L1/2";
const LABEL_PRN_OBS: &str = "PRN / # OF OBS";
const LABEL_END: &str = "END OF HEADER";

/// Read the format version from the `RINEX VERSION / TYPE` record without
/// parsing the rest of the header.  Returns `None` if no such record appears
/// before `END OF HEADER` or its version field is unreadable.
#[must_use]
pub fn sniff_version(bytes: &[u8]) -> Option<Version> {
    for line in bytes.split(|&b| b == b'\n') {
        let label = line.get(60..).unwrap_or_default();
        let label = std::str::from_utf8(label).unwrap_or_default().trim_end();
        if label == LABEL_END {
            break;
        }
        if label != LABEL_VERSION {
            continue;
        }
        let field = std::str::from_utf8(line.get(..9)?).ok()?.trim();
        let (major, minor) = field.split_once('.').unwrap_or((field, ""));
        let digits: String = minor.chars().take_while(char::is_ascii_digit).collect();
        return Some(Version {
            major: major.parse().ok()?,
            minor: digits.parse().unwrap_or(0),
        });
    }
    None
}

#[doc(hidden)]
pub fn parse(bytes: &[u8]) -> Result<Header, HeaderError> {
    parse_str(&to_text(bytes))
}

/// `bytes` as text for fixed-column parsing.  A file with any invalid UTF-8
/// byte is mapped to one `?` per byte, so its columns (and therefore its
/// header labels) stay where the file put them; a file that is valid UTF-8,
/// the usual case, is borrowed.
pub(crate) fn to_text(bytes: &[u8]) -> std::borrow::Cow<'_, str> {
    match std::str::from_utf8(bytes) {
        Ok(text) => std::borrow::Cow::Borrowed(text),
        Err(_) => std::borrow::Cow::Owned(
            bytes
                .iter()
                .map(|&b| if b < 0x80 { char::from(b) } else { '?' })
                .collect(),
        ),
    }
}

/// [`parse`] for text already known to be UTF-8.
pub(crate) fn parse_str(text: &str) -> Result<Header, HeaderError> {
    let mut state = Parser::default();
    let mut version_seen = false;

    for raw_line in text.lines() {
        // Tolerate lingering CR.
        let line = raw_line.strip_suffix('\r').unwrap_or(raw_line);
        // Fixed-column slicing below assumes one byte per column; map any
        // non-ASCII character (e.g. in a COMMENT) to a single placeholder.
        let ascii;
        let line = if line.is_ascii() {
            line
        } else {
            ascii = line
                .chars()
                .map(|c| if c.is_ascii() { c } else { '?' })
                .collect::<String>();
            ascii.as_str()
        };
        let Some((content, label)) = split_label(line) else {
            // Lines shorter than 60 chars can only be the END OF HEADER
            // sentinel in degenerate inputs; skip silently otherwise.
            continue;
        };

        if label == LABEL_END {
            state.finalize_pending_obs(&state.last_obs_label.clone())?;
            return Ok(state.into_header());
        }

        if label == LABEL_VERSION {
            parse_version_line(content, &mut state.header)?;
            version_seen = true;
            state.last_obs_label = String::new();
            continue;
        }

        if !version_seen {
            // RINEX requires VERSION to be the first line. Be lenient and
            // accept missing-but-not-yet-required state, but skip dispatch
            // on records we can't interpret without knowing the version.
            // We loop and let the END OF HEADER check above catch
            // truly-missing-version cases as MissingVersion.
            continue;
        }

        // If switching away from an obs-codes label, finalize the pending
        // accumulator (catches the obs codes ending exactly at the count).
        if !state.last_obs_label.is_empty() && label != state.last_obs_label {
            let prev = state.last_obs_label.clone();
            state.finalize_pending_obs(&prev)?;
            state.last_obs_label.clear();
        }

        match label {
            LABEL_COMMENT => state.header.comments.push(content.trim_end().to_string()),
            LABEL_PGM => parse_pgm(content, &mut state.header),
            LABEL_MARKER_NAME => state.header.marker_name = Some(trim60(content)),
            LABEL_MARKER_NUMBER => state.header.marker_number = Some(trim60(content)),
            LABEL_MARKER_TYPE => state.header.marker_type = Some(trim60(content)),
            LABEL_OBSERVER => parse_observer(content, &mut state.header),
            LABEL_RECEIVER => parse_receiver(content, &mut state.header),
            LABEL_ANTENNA => parse_antenna(content, &mut state.header),
            LABEL_APPROX_POS => state.header.approx_position = parse_xyz(content),
            LABEL_ANT_DELTA => state.header.antenna_delta_hen = parse_xyz(content),
            LABEL_INTERVAL => {
                state.header.interval = parse_fortran_float(&content[..20.min(content.len())]);
            }
            LABEL_TIME_FIRST => state.header.time_of_first_obs = parse_time_of_obs(content),
            LABEL_TIME_LAST => state.header.time_of_last_obs = parse_time_of_obs(content),
            LABEL_LEAP => {
                state.header.leap_seconds = content[..6.min(content.len())].trim().parse().ok();
            }
            LABEL_WAVELEN => state.header.wavelength_fact_l1l2 = parse_wavelength_fact(content),
            // 3X,A1,I2: the satellite, blank on continuation lines.
            LABEL_PRN_OBS => {
                if let Some(sat) = content.get(3..).and_then(parse_sat_id) {
                    if !state.header.satellites.contains(&sat) {
                        state.header.satellites.push(sat);
                    }
                }
            }

            LABEL_OBS_TYPES_V2 => {
                state.handle_obs_v2_line(content)?;
                state.last_obs_label = LABEL_OBS_TYPES_V2.to_string();
            }
            LABEL_OBS_TYPES_V3 => {
                state.handle_obs_v3_line(content)?;
                state.last_obs_label = LABEL_OBS_TYPES_V3.to_string();
            }

            other if !other.is_empty() => {
                state
                    .header
                    .unknown_lines
                    .insert(other.to_string(), trim60(content));
            }
            _ => {}
        }
    }

    Err(HeaderError::NoEndOfHeader)
}

/// The observation codes declared by the header records of a special event
/// (flags 2–5): its `# / TYPES OF OBSERV` (RINEX 2) or `SYS / # / OBS TYPES`
/// (RINEX 3 and 4) records, which change the codes of the epochs after it.
/// Empty if it has none.
pub(crate) fn event_obs_codes<'a>(
    lines: impl IntoIterator<Item = &'a str>,
) -> Result<ObsCodes, HeaderError> {
    let mut state = Parser::default();
    let mut last = "";
    for line in lines {
        let label = split_label(line).map_or("", |(_, label)| label);
        if !last.is_empty() && label != last {
            state.finalize_pending_obs(last)?;
        }
        last = "";
        if let Some((content, label)) = split_label(line) {
            if label == LABEL_OBS_TYPES_V2 {
                state.handle_obs_v2_line(content)?;
                last = LABEL_OBS_TYPES_V2;
            } else if label == LABEL_OBS_TYPES_V3 {
                state.handle_obs_v3_line(content)?;
                last = LABEL_OBS_TYPES_V3;
            }
        }
    }
    state.finalize_pending_obs(last)?;
    Ok(state.header.obs_codes)
}

/// Check a completed obs-code declaration: it has the declared number of
/// codes, and no code twice.  Every reader keys observations by satellite
/// and code, so a repeated code would give one series two values per epoch
/// (`rnx2srnx` rejects such a file too).
fn check_codes(system: Option<u8>, expected: usize, codes: &[ObsCode]) -> Result<(), HeaderError> {
    if codes.len() != expected {
        return Err(HeaderError::ObsCountMismatch {
            system,
            expected,
            actual: codes.len(),
        });
    }
    for (i, code) in codes.iter().enumerate() {
        if codes[..i].contains(code) {
            return Err(HeaderError::DuplicateObsCode {
                system,
                code: code.as_str().to_owned(),
            });
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Parser state
// ---------------------------------------------------------------------------

#[derive(Default)]
struct Parser {
    header: Header,
    /// In-flight v2 obs-code accumulator: (declared count, codes so far).
    obs_v2: Option<(usize, Vec<ObsCode>)>,
    /// In-flight v3 obs-code accumulator: (system, declared count, codes so far).
    obs_v3_current: Option<(SatSystem, usize, Vec<ObsCode>)>,
    /// Label of the most recently processed obs-codes line; used to detect
    /// end of a multi-line block when we see a different label.
    last_obs_label: String,
}

impl Parser {
    fn into_header(self) -> Header {
        self.header
    }

    fn handle_obs_v2_line(&mut self, content: &str) -> Result<(), HeaderError> {
        // First 6 cols: I6 count (blank on continuation lines).
        let count_field = &content[..6.min(content.len())];
        let count_str = count_field.trim();

        if !count_str.is_empty() {
            // New v2 obs-codes block: finalize any prior one first.
            self.finalize_v2()?;
            let expected: usize = count_str.parse().map_err(|_| HeaderError::Malformed {
                label: "# / TYPES OF OBSERV",
                detail: format!("count {count_str:?}"),
            })?;
            self.obs_v2 = Some((expected, Vec::with_capacity(expected)));
        } else if self.obs_v2.is_none() {
            return Err(HeaderError::OrphanObsContinuation);
        }

        // 9 codes per line, each occupying (4X, A2) starting at col 7 (index 6).
        for slot in 0..9 {
            let start = 6 + slot * 6;
            if start + 6 > content.len() {
                break;
            }
            let field = &content[start..start + 6];
            let trimmed = field.trim();
            if trimmed.is_empty() {
                continue;
            }
            let code = ObsCode::from_v2(trimmed).ok_or_else(|| HeaderError::Malformed {
                label: "# / TYPES OF OBSERV",
                detail: format!("invalid v2 obs code {trimmed:?}"),
            })?;
            self.obs_v2
                .as_mut()
                .expect("v2 accumulator exists")
                .1
                .push(code);
        }
        Ok(())
    }

    fn handle_obs_v3_line(&mut self, content: &str) -> Result<(), HeaderError> {
        // Col 1: A1 system code (blank = continuation).
        // Cols 4-6: I3 count.
        // Cols 7-58: 13 × (1X, A3).
        let bytes = content.as_bytes();
        let leading = *bytes.first().unwrap_or(&b' ');

        if leading != b' ' {
            // New v3 obs-codes block for a new system: finalize any prior.
            self.finalize_v3_current()?;
            let count_field = &content[3..6.min(content.len())];
            let count_str = count_field.trim();
            let expected: usize = count_str.parse().map_err(|_| HeaderError::Malformed {
                label: "SYS / # / OBS TYPES",
                detail: format!("count {count_str:?}"),
            })?;
            self.obs_v3_current =
                Some((SatSystem(leading), expected, Vec::with_capacity(expected)));
        } else if self.obs_v3_current.is_none() {
            return Err(HeaderError::OrphanObsContinuation);
        }

        for slot in 0..13 {
            let start = 6 + slot * 4;
            if start + 4 > content.len() {
                break;
            }
            let field = &content[start..start + 4];
            let trimmed = field.trim();
            if trimmed.is_empty() {
                continue;
            }
            let code = ObsCode::from_v3(trimmed).ok_or_else(|| HeaderError::Malformed {
                label: "SYS / # / OBS TYPES",
                detail: format!("invalid v3 obs code {trimmed:?}"),
            })?;
            self.obs_v3_current
                .as_mut()
                .expect("v3 accumulator exists")
                .2
                .push(code);
        }
        Ok(())
    }

    fn finalize_v2(&mut self) -> Result<(), HeaderError> {
        if let Some((expected, codes)) = self.obs_v2.take() {
            check_codes(None, expected, &codes)?;
            self.header
                .obs_codes
                .by_system
                .push((SatSystem::MIXED, codes));
        }
        Ok(())
    }

    fn finalize_v3_current(&mut self) -> Result<(), HeaderError> {
        if let Some((sys, expected, codes)) = self.obs_v3_current.take() {
            check_codes(Some(sys.0), expected, &codes)?;
            self.header.obs_codes.by_system.push((sys, codes));
        }
        Ok(())
    }

    fn finalize_pending_obs(&mut self, last_label: &str) -> Result<(), HeaderError> {
        match last_label {
            LABEL_OBS_TYPES_V2 => self.finalize_v2()?,
            LABEL_OBS_TYPES_V3 => self.finalize_v3_current()?,
            _ => {}
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Per-line parsers
// ---------------------------------------------------------------------------

fn trim60(content: &str) -> String {
    content.trim_end().to_string()
}

fn parse_version_line(content: &str, header: &mut Header) -> Result<(), HeaderError> {
    if content.len() < 41 {
        return Err(HeaderError::BadVersion(format!(
            "line too short ({} cols)",
            content.len()
        )));
    }
    let version_field = content[..9].trim();
    let (major_str, minor_str) = match version_field.split_once('.') {
        Some(pair) => pair,
        None => (version_field, ""),
    };
    let major: u8 = major_str
        .trim()
        .parse()
        .map_err(|_| HeaderError::BadVersion(format!("major {major_str:?}")))?;
    let minor: u8 = if minor_str.is_empty() {
        0
    } else {
        // Take the leading digits; "11" or "0" or "04" are all valid.
        let digits: String = minor_str.chars().take_while(char::is_ascii_digit).collect();
        digits.parse().unwrap_or(0)
    };
    if major != 2 && major != 3 && major != 4 {
        return Err(HeaderError::UnsupportedVersion(major));
    }

    let file_type = content.as_bytes().get(20).copied().unwrap_or(b' ');
    if file_type != b'O' {
        return Err(HeaderError::NotObservationFile {
            got: file_type as char,
        });
    }
    let sat_system = content.as_bytes().get(40).copied().unwrap_or(b' ');

    header.version = Version { major, minor };
    header.file_type = file_type;
    header.sat_system = SatSystem(sat_system);
    Ok(())
}

fn parse_pgm(content: &str, header: &mut Header) {
    // 3 fields × 20 cols = 60 cols.
    let pgm = content.get(..20).map_or("", str::trim_end);
    let run_by = content.get(20..40).map_or("", str::trim_end);
    let date = content.get(40..60).map_or("", str::trim_end);
    if !pgm.is_empty() {
        header.program = Some(pgm.to_string());
    }
    if !run_by.is_empty() {
        header.run_by = Some(run_by.to_string());
    }
    if !date.is_empty() {
        header.date_str = Some(date.to_string());
    }
}

fn parse_observer(content: &str, header: &mut Header) {
    let observer = content.get(..20).map_or("", str::trim_end);
    let agency = content.get(20..60).map_or("", str::trim_end);
    if !observer.is_empty() {
        header.observer = Some(observer.to_string());
    }
    if !agency.is_empty() {
        header.agency = Some(agency.to_string());
    }
}

fn parse_receiver(content: &str, header: &mut Header) {
    let serial = content.get(..20).map_or("", str::trim_end);
    let model = content.get(20..40).map_or("", str::trim_end);
    let version = content.get(40..60).map_or("", str::trim_end);
    header.receiver = Some(Receiver {
        serial: serial.to_string(),
        model: model.to_string(),
        version: version.to_string(),
    });
}

fn parse_antenna(content: &str, header: &mut Header) {
    let serial = content.get(..20).map_or("", str::trim_end);
    let model = content.get(20..40).map_or("", str::trim_end);
    header.antenna = Some(Antenna {
        serial: serial.to_string(),
        model: model.to_string(),
    });
}

/// Parse three F14.4 fields (cols 1-14, 15-28, 29-42).
fn parse_xyz(content: &str) -> Option<[f64; 3]> {
    let a = parse_fortran_float(content.get(..14)?)?;
    let b = parse_fortran_float(content.get(14..28)?)?;
    let c = parse_fortran_float(content.get(28..42)?)?;
    Some([a, b, c])
}

fn parse_wavelength_fact(content: &str) -> Option<(u8, u8)> {
    let l1 = content.get(..6)?.trim().parse().ok()?;
    let l2 = content.get(6..12)?.trim().parse().ok()?;
    Some((l1, l2))
}

/// Parse `TIME OF FIRST OBS` / `TIME OF LAST OBS`. Format is:
/// `5I6, F13.7, 5X, A3` (year/month/day/hour/minute, seconds, time system).
fn parse_time_of_obs(content: &str) -> Option<EpochTime> {
    let year: i32 = content.get(..6)?.trim().parse().ok()?;
    let month: u8 = content.get(6..12)?.trim().parse().ok()?;
    let day: u8 = content.get(12..18)?.trim().parse().ok()?;
    let hour: u8 = content.get(18..24)?.trim().parse().ok()?;
    let minute: u8 = content.get(24..30)?.trim().parse().ok()?;
    let seconds = parse_fortran_float(content.get(30..43)?)?;
    let time_system = content
        .get(48..51)
        .map_or(String::new(), |s| s.trim().to_string());
    // Convert seconds to seconds_e7 with rounding at the boundary.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let seconds_e7 = (seconds * 1e7).round() as u64;
    Some(EpochTime {
        year,
        month,
        day,
        hour,
        minute,
        seconds_e7,
        time_system,
    })
}

/// Tolerate Fortran-style D-exponents (e.g. `1.234D-5`) in addition to E.
fn parse_fortran_float(s: &str) -> Option<f64> {
    let trimmed = s.trim();
    if trimmed.is_empty() {
        return None;
    }
    if trimmed.bytes().any(|b| b == b'D' || b == b'd') {
        let normalized: String = trimmed
            .chars()
            .map(|c| if c == 'D' || c == 'd' { 'E' } else { c })
            .collect();
        normalized.parse().ok()
    } else {
        trimmed.parse().ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Pad a content string to exactly 60 chars and append a label, producing
    /// one well-formed RINEX header line (no newline).
    fn line(content: &str, label: &str) -> String {
        let mut s = String::with_capacity(80);
        s.push_str(content);
        if s.len() < 60 {
            s.push_str(&" ".repeat(60 - s.len()));
        } else {
            s.truncate(60);
        }
        s.push_str(label);
        s
    }

    fn build(lines: &[&str]) -> String {
        let mut out = String::new();
        for l in lines {
            out.push_str(l);
            out.push('\n');
        }
        out
    }

    #[test]
    fn parse_minimal_v2() {
        let header = build(&[
            &line(
                "     2.11           OBSERVATION DATA    M (MIXED)           ",
                "RINEX VERSION / TYPE",
            ),
            &line(
                "     2    L1    C1                                          ",
                "# / TYPES OF OBSERV",
            ),
            &line("", "END OF HEADER"),
        ]);
        let h = parse(header.as_bytes()).unwrap();
        assert_eq!(
            h.version,
            Version {
                major: 2,
                minor: 11
            }
        );
        assert_eq!(h.file_type, b'O');
        assert_eq!(h.sat_system, SatSystem::MIXED);
        let codes = h.obs_codes.for_system(SatSystem::MIXED).unwrap();
        assert_eq!(codes.len(), 2);
        assert_eq!(codes[0].as_str(), "L1");
        assert_eq!(codes[1].as_str(), "C1");
    }

    #[test]
    fn prn_obs_records_list_satellites() {
        let header = build(&[
            &line(
                "     2.11           OBSERVATION DATA    M (MIXED)           ",
                "RINEX VERSION / TYPE",
            ),
            &line(
                "     2    L1    C1                                          ",
                "# / TYPES OF OBSERV",
            ),
            &line("   E01   971   971", "PRN / # OF OBS"),
            &line("         971", "PRN / # OF OBS"),
            &line("    12   500   500", "PRN / # OF OBS"),
            &line("   E01     1     1", "PRN / # OF OBS"),
            &line("   R05   100      ", "PRN / # OF OBS"),
            &line("", "END OF HEADER"),
        ]);
        let h = parse(header.as_bytes()).unwrap();
        assert_eq!(
            h.satellites,
            [
                SatId::new(SatSystem::GALILEO, 1),
                SatId::new(SatSystem::GPS, 12),
                SatId::new(SatSystem::GLONASS, 5),
            ]
        );
    }

    #[test]
    fn parse_v2_with_continuation() {
        // 12 obs codes = 1 continuation line.
        let header = build(&[
            &line(
                "     2.10           OBSERVATION DATA    M                   ",
                "RINEX VERSION / TYPE",
            ),
            &line(
                "    12    L1    L2    C1    P2    P1    S1    S2    C2    L5",
                "# / TYPES OF OBSERV",
            ),
            &line(
                "          C5    S5    L6                                    ",
                "# / TYPES OF OBSERV",
            ),
            &line("", "END OF HEADER"),
        ]);
        let h = parse(header.as_bytes()).unwrap();
        let codes = h.obs_codes.for_system(SatSystem::MIXED).unwrap();
        assert_eq!(codes.len(), 12);
        let names: Vec<&str> = codes.iter().map(ObsCode::as_str).collect();
        assert_eq!(
            names,
            vec!["L1", "L2", "C1", "P2", "P1", "S1", "S2", "C2", "L5", "C5", "S5", "L6"]
        );
    }

    #[test]
    fn parse_v3_two_systems() {
        let header = build(&[
            &line(
                "     3.04           OBSERVATION DATA    M                   ",
                "RINEX VERSION / TYPE",
            ),
            &line(
                "G    4 C1C L1C D1C S1C                                      ",
                "SYS / # / OBS TYPES",
            ),
            &line(
                "E    4 C1C L1C D1C S1C                                      ",
                "SYS / # / OBS TYPES",
            ),
            &line("", "END OF HEADER"),
        ]);
        let h = parse(header.as_bytes()).unwrap();
        assert_eq!(h.version, Version { major: 3, minor: 4 });
        assert_eq!(h.obs_codes.by_system.len(), 2);
        assert_eq!(h.obs_codes.by_system[0].0, SatSystem::GPS);
        assert_eq!(h.obs_codes.by_system[1].0, SatSystem::GALILEO);
        assert_eq!(h.obs_codes.for_system(SatSystem::GPS).unwrap().len(), 4);
    }

    #[test]
    fn parse_v3_with_continuation() {
        // 18 obs codes for G — needs one continuation line (13 + 5).
        let header = build(&[
            &line(
                "     3.04           OBSERVATION DATA    M                   ",
                "RINEX VERSION / TYPE",
            ),
            &line(
                "G   18 C1C L1C D1C S1C C1W S1W C2W L2W D2W S2W C2L L2L D2L  ",
                "SYS / # / OBS TYPES",
            ),
            &line(
                "       S2L C5Q L5Q D5Q S5Q                                  ",
                "SYS / # / OBS TYPES",
            ),
            &line("", "END OF HEADER"),
        ]);
        let h = parse(header.as_bytes()).unwrap();
        let codes = h.obs_codes.for_system(SatSystem::GPS).unwrap();
        assert_eq!(codes.len(), 18);
        assert_eq!(codes[0].as_str(), "C1C");
        assert_eq!(codes[12].as_str(), "D2L");
        assert_eq!(codes[13].as_str(), "S2L");
        assert_eq!(codes[17].as_str(), "S5Q");
    }

    #[test]
    fn detects_count_mismatch() {
        // Says 5 codes, only gives 2.
        let header = build(&[
            &line(
                "     2.11           OBSERVATION DATA    G                   ",
                "RINEX VERSION / TYPE",
            ),
            &line(
                "     5    L1    C1                                          ",
                "# / TYPES OF OBSERV",
            ),
            &line("", "END OF HEADER"),
        ]);
        let err = parse(header.as_bytes()).unwrap_err();
        assert!(matches!(
            err,
            HeaderError::ObsCountMismatch {
                system: None,
                expected: 5,
                actual: 2
            }
        ));
    }

    #[test]
    fn rejects_duplicate_codes() {
        for (version, types) in [
            (
                "     2.11           OBSERVATION DATA    G                   ",
                line(
                    "     3    L1    C1    L1                                    ",
                    "# / TYPES OF OBSERV",
                ),
            ),
            (
                "     3.04           OBSERVATION DATA    M                   ",
                line("G    2 C1C C1C", "SYS / # / OBS TYPES"),
            ),
        ] {
            let header = build(&[
                &line(version, "RINEX VERSION / TYPE"),
                &types,
                &line("", "END OF HEADER"),
            ]);
            let err = parse(header.as_bytes()).unwrap_err();
            assert!(
                matches!(&err, HeaderError::DuplicateObsCode { code, .. } if code == "L1" || code == "C1C"),
                "{err:?}"
            );
        }
    }

    #[test]
    fn sniffs_version() {
        let v = |first: &str| {
            let text = build(&[
                &line(first, "RINEX VERSION / TYPE"),
                &line("", "END OF HEADER"),
            ]);
            sniff_version(text.as_bytes())
        };
        let want = |major, minor| Some(Version { major, minor });
        assert_eq!(v("     2.11           OBSERVATION DATA    M"), want(2, 11));
        assert_eq!(v("     3.05           OBSERVATION DATA    M"), want(3, 5));
        assert_eq!(v("     4.02           OBSERVATION DATA    M"), want(4, 2));
        assert_eq!(v("   junk"), None);
        assert_eq!(sniff_version(b"SRNX\x01"), None);
        assert_eq!(sniff_version(&[0xff; 200]), None);
    }

    #[test]
    fn rejects_non_observation_files() {
        let header = build(&[
            &line(
                "     2.11           NAVIGATION DATA                         ",
                "RINEX VERSION / TYPE",
            ),
            &line("", "END OF HEADER"),
        ]);
        let err = parse(header.as_bytes()).unwrap_err();
        assert!(matches!(err, HeaderError::NotObservationFile { got: 'N' }));
    }

    #[test]
    fn rejects_missing_end_of_header() {
        let header = build(&[&line(
            "     2.11           OBSERVATION DATA    G                   ",
            "RINEX VERSION / TYPE",
        )]);
        let err = parse(header.as_bytes()).unwrap_err();
        assert!(matches!(err, HeaderError::NoEndOfHeader));
    }

    #[test]
    fn parses_metadata_fields() {
        let header = build(&[
            &line(
                "     2.11           OBSERVATION DATA    G                   ",
                "RINEX VERSION / TYPE",
            ),
            &line(
                "AT01                                                        ",
                "MARKER NAME",
            ),
            &line(
                "AOOS Community      Alaska Ocean Observing System           ",
                "OBSERVER / AGENCY",
            ),
            &line(
                "3025329             SEPT POLARX5        5.2.0               ",
                "REC # / TYPE / VERS",
            ),
            &line(
                "5732337004          TRM159800.00    SCIT                    ",
                "ANT # / TYPE",
            ),
            &line(
                " -2715527.3710  -881994.3010  5684288.5300                  ",
                "APPROX POSITION XYZ",
            ),
            &line(
                "        0.0083        0.0000        0.0000                  ",
                "ANTENNA: DELTA H/E/N",
            ),
            &line(
                "     1.0000                                                 ",
                "INTERVAL",
            ),
            &line(
                "  2020     7    18     0     0    0.0000000     GPS         ",
                "TIME OF FIRST OBS",
            ),
            &line(
                "    18                                                      ",
                "LEAP SECONDS",
            ),
            &line(
                "     2    L1    C1                                          ",
                "# / TYPES OF OBSERV",
            ),
            &line("", "END OF HEADER"),
        ]);
        let h = parse(header.as_bytes()).unwrap();
        assert_eq!(h.marker_name.as_deref(), Some("AT01"));
        assert_eq!(h.observer.as_deref(), Some("AOOS Community"));
        assert_eq!(h.agency.as_deref(), Some("Alaska Ocean Observing System"));
        let r = h.receiver.as_ref().unwrap();
        assert_eq!(r.serial, "3025329");
        assert_eq!(r.model, "SEPT POLARX5");
        assert_eq!(r.version, "5.2.0");
        let a = h.antenna.as_ref().unwrap();
        assert_eq!(a.serial, "5732337004");
        assert_eq!(a.model, "TRM159800.00    SCIT");
        let pos = h.approx_position.unwrap();
        assert!((pos[0] - -2_715_527.371).abs() < 1e-3);
        assert!((pos[1] - -881_994.301).abs() < 1e-3);
        assert!((pos[2] - 5_684_288.53).abs() < 1e-3);
        let dh = h.antenna_delta_hen.unwrap();
        assert!((dh[0] - 0.0083).abs() < 1e-7);
        assert_eq!(h.interval, Some(1.0));
        let t = h.time_of_first_obs.as_ref().unwrap();
        assert_eq!(t.year, 2020);
        assert_eq!(t.month, 7);
        assert_eq!(t.day, 18);
        assert_eq!(t.hour, 0);
        assert_eq!(t.minute, 0);
        assert_eq!(t.seconds_e7, 0);
        assert_eq!(t.time_system, "GPS");
        assert_eq!(h.leap_seconds, Some(18));
    }

    /// The `RHDR` payload of an SRNX fixture.
    fn fixture_rhdr(name: &str) -> Vec<u8> {
        let bytes = crate::test_data(name);
        let file = crate::file::SrnxFile::parse(&bytes).unwrap();
        file.rhdr_payload().unwrap().to_vec()
    }

    #[test]
    fn fixture_at01_v2() {
        let h = parse(&fixture_rhdr("at01_v211_1hz.srnx")).unwrap();
        assert_eq!(
            h.version,
            Version {
                major: 2,
                minor: 11
            }
        );
        assert_eq!(h.sat_system, SatSystem::MIXED);
        let codes = h.obs_codes.for_system(SatSystem::MIXED).unwrap();
        // From the file: "20 L1 L2 C1 P2 P1 S1 S2 C2 L5 C5 S5 L6 C6 S6 L7 C7 S7 L8 C8 S8".
        assert_eq!(codes.len(), 20);
        let names: Vec<&str> = codes.iter().map(ObsCode::as_str).collect();
        assert_eq!(
            names,
            vec![
                "L1", "L2", "C1", "P2", "P1", "S1", "S2", "C2", "L5", "C5", "S5", "L6", "C6", "S6",
                "L7", "C7", "S7", "L8", "C8", "S8"
            ]
        );
        assert_eq!(h.marker_name.as_deref(), Some("AT01"));
        assert_eq!(h.interval, Some(1.0));
        assert_eq!(h.leap_seconds, Some(18));
        let t = h.time_of_first_obs.as_ref().unwrap();
        assert_eq!((t.year, t.month, t.day), (2020, 7, 18));
        assert_eq!(t.time_system, "GPS");
    }

    #[test]
    fn fixture_balm_v3() {
        let h = parse(&fixture_rhdr("balm_v304_1hz.srnx")).unwrap();
        assert_eq!(h.version, Version { major: 3, minor: 4 });
        assert_eq!(h.sat_system, SatSystem::MIXED);
        let systems: Vec<u8> = h.obs_codes.by_system.iter().map(|(s, _)| s.0).collect();
        assert_eq!(systems, b"GREJCI");
        let count = |sys| h.obs_codes.for_system(sys).unwrap().len();
        assert_eq!(count(SatSystem::GPS), 24);
        assert_eq!(count(SatSystem::GLONASS), 20);
        assert_eq!(count(SatSystem::GALILEO), 20);
        assert_eq!(count(SatSystem::QZSS), 24);
        assert_eq!(count(SatSystem::BEIDOU), 20);
        assert_eq!(count(SatSystem::IRNSS), 8);
        assert_eq!(h.marker_name.as_deref(), Some("BALM"));
        assert_eq!(h.interval, Some(1.0));
        assert_eq!(h.leap_seconds, Some(18));
    }

    // ── HeaderText tests ─────────────────────────────────────────────────

    fn hline(content: &str, label: &str) -> String {
        let mut s = content.to_string();
        if s.len() < 60 {
            s.push_str(&" ".repeat(60 - s.len()));
        } else {
            s.truncate(60);
        }
        s.push_str(label);
        s
    }

    #[test]
    fn header_text_from_bytes_normalizes_crlf() {
        let bytes = b"line one\r\nline two\r\n";
        let h = HeaderText::from_bytes(bytes);
        assert_eq!(h.as_str(), "line one\nline two\n");
    }

    #[test]
    fn header_text_get_line() {
        let text = format!(
            "{}\n{}\n{}\n",
            hline(
                "     2.11           OBSERVATION DATA    M                   ",
                "RINEX VERSION / TYPE"
            ),
            hline(
                "AT01                                                        ",
                "MARKER NAME"
            ),
            hline("", "END OF HEADER")
        );
        let h = HeaderText::from_text(&text);
        assert_eq!(h.get_line("MARKER NAME"), Some("AT01".to_string()));
        assert_eq!(h.get_line("NONEXISTENT"), None);
    }

    #[test]
    fn header_text_get_line_trailing_ws_invariant() {
        // Label stored with extra trailing spaces should still match.
        let text = format!(
            "{}\n{}\n",
            hline(
                "AT01                                                        ",
                "MARKER NAME   "
            ),
            hline("", "END OF HEADER")
        );
        let h = HeaderText::from_text(&text);
        assert_eq!(h.get_line("MARKER NAME"), Some("AT01".to_string()));
    }

    #[test]
    fn header_text_get_line_last_wins() {
        // When a label appears multiple times, last occurrence wins.
        let text = format!(
            "{}\n{}\n{}\n",
            hline(
                "OLD NAME                                                      ",
                "MARKER NAME"
            ),
            hline(
                "NEW NAME                                                      ",
                "MARKER NAME"
            ),
            hline("", "END OF HEADER")
        );
        let h = HeaderText::from_text(&text);
        assert_eq!(h.get_line("MARKER NAME"), Some("NEW NAME".to_string()));
    }

    #[test]
    fn header_text_get_lines() {
        let text = format!(
            "{}\n{}\n{}\n",
            hline(
                "first comment                                                   ",
                "COMMENT"
            ),
            hline(
                "second comment                                                  ",
                "COMMENT"
            ),
            hline("", "END OF HEADER")
        );
        let h = HeaderText::from_text(&text);
        let lines: Vec<String> = h.get_lines("COMMENT").collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], "first comment");
        assert_eq!(lines[1], "second comment");
    }

    #[test]
    fn header_text_deref() {
        let text = "some header\nEND OF HEADER\n";
        let h = HeaderText::from_text(text);
        // Deref to str works.
        assert!(h.contains("header"));
        assert!(h.starts_with("some"));
    }

    #[test]
    fn header_text_fixture_at01() {
        let h = HeaderText::from_bytes(&fixture_rhdr("at01_v211_1hz.srnx"));
        assert_eq!(h.get_line("MARKER NAME"), Some("AT01".to_string()));
        // Multiple COMMENT lines.
        let comments: Vec<String> = h.get_lines("COMMENT").collect();
        assert_ne!(comments, [] as [String; 0]);
    }
}
