//! Parser for RINEX v2.x observation data files.
//!
//! Reads the fixed-column text format defined in `docs/rinex211.txt`,
//! Table A2.  Populates the shared [`ObsFile`] semantic model.
//!
//! ## Format summary
//!
//! ### Epoch / satellite header line (80 chars)
//!
//! ```text
//!  YY MM DD HH MM SS.SSSSSSS  F NNN S01S02...S12  CCCCCCCCCCCC
//!  ^col1                                            ^clock offset
//! ```
//!
//! Fortran: `1X, I2.2, 4(1X,I2), F11.7, 2X, I1, I3, 12(A1,I2), F12.9`
//!
//! Continuation line for >12 satellites:
//! `32X, 12(A1,I2)`
//!
//! ### Observation records (one set per satellite in epoch order)
//!
//! `m(F14.3, I1, I1)` per line, up to 5 observations (80 chars total).
//! If a satellite has more than 5 observation codes the record wraps to a
//! second (or more) 80-char line.  Lines shorter than the full 80 chars
//! have been right-truncated; treat missing characters as spaces.
//!
//! ### Event records (epoch flag 2–5)
//!
//! The "number of satellites" field instead gives the number of verbatim
//! lines that follow the epoch header.
//!
//! ### Power-failure and cycle-slip records (epoch flags 1 and 6)
//!
//! Same layout as a flag-0 epoch.  As in SRNX (which stores them only in
//! `EVTF`), they become events rather than observation epochs.

use std::collections::HashSet;

use thiserror::Error;

use crate::model::{
    Catalog, EpochRecord, EventRecord, Obs, ObsFile, SatId, SeriesBuilder, SignalFilter,
    SignalMasks,
};
use crate::rinex3_obs::parse_sat_id;
use crate::rinex_header::{self, Header, HeaderError, HeaderText, ObsCode, SatSystem};

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Rinex2Error {
    #[error("header parse error: {0}")]
    Header(#[from] HeaderError),
    #[error("input is not valid UTF-8")]
    NotUtf8,
    #[error("no END OF HEADER line found")]
    NoEndOfHeader,
    #[error("epoch line is too short to contain a timestamp (line {0})")]
    EpochLineTooShort(usize),
    #[error("invalid epoch flag or satellite count on line {0}")]
    BadEpochHeader(usize),
    #[error("observation file has no obs codes in header")]
    NoObsCodes,
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// List the satellites and signals declared by a RINEX v2 observation
/// file's header.  Every satellite may carry every observation code.
pub fn catalog(input: &[u8]) -> Result<Catalog, Rinex2Error> {
    let header = rinex_header::parse(input)?;
    let codes = obs_codes(&header)?;
    Ok(Catalog::from_header(&header, |_| codes))
}

/// Parse a complete RINEX v2 observation file, keeping only the
/// satellite/signal pairs that `filter` selects (all of them for `None`).
///
/// Both LF and CRLF line endings are accepted.  Malformed epochs and stray
/// lines are skipped and counted in [`ObsFile::skipped_records`]; the next
/// plausible epoch header resumes parsing.
#[allow(clippy::too_many_lines)]
pub fn parse(input: &[u8], filter: Option<SignalFilter<'_>>) -> Result<ObsFile, Rinex2Error> {
    let text = std::str::from_utf8(input).map_err(|_| Rinex2Error::NotUtf8)?;

    // Split into lines (strips the newline character(s)).
    let lines = split_lines(text);

    // Locate the END OF HEADER line.
    let eoh_idx = lines
        .iter()
        .position(|l| {
            let trimmed = if l.len() >= 60 {
                l.get(60..).unwrap_or("").trim()
            } else {
                l.trim()
            };
            trimmed.starts_with("END OF HEADER")
        })
        .ok_or(Rinex2Error::NoEndOfHeader)?;

    // Capture raw header text (all lines through END OF HEADER), normalized to LF.
    // `lines()` strips newlines, so re-add the trailing LF to match RHDR semantics.
    let mut raw_header_text = lines[..=eoh_idx].join("\n");
    raw_header_text.push('\n');

    // Parse the header section (pass the whole text — the header parser
    // tolerates trailing data because it stops at END OF HEADER).
    let header = rinex_header::parse_str(text).map_err(Rinex2Error::Header)?;

    // An event's `# / TYPES OF OBSERV` records change these.
    let mut obs_codes = obs_codes(&header)?.to_vec();
    let mut num_obs = obs_codes.len();
    let mut masks = SignalMasks::new(filter);
    for &sat in &header.satellites {
        masks.get(sat, &obs_codes);
    }
    let mut obs_lines_per_sat = obs_lines_count(num_obs);

    // Data starts on the line after END OF HEADER.
    let data_lines = &lines[eoh_idx + 1..];
    let base_line = eoh_idx + 1; // for error reporting

    // Sparse per-(sat, code) → [(epoch_index, Obs)] list, appended in order.
    let mut series = SeriesBuilder::default();
    let mut epochs: Vec<EpochRecord> = Vec::new();
    let mut events: Vec<(usize, EventRecord)> = Vec::new();
    let mut skipped_records = 0_usize;
    let mut seen_sats: HashSet<SatId> = HashSet::new();

    let mut pos = 0; // index into data_lines
    while pos < data_lines.len() {
        let line = data_lines[pos];
        let line_no = base_line + pos + 1; // 1-based line number for errors

        // Parse the epoch / satellite header line.
        let Ok((epoch_rec, flag, num_entries)) = parse_epoch_header(line, line_no) else {
            // Not a valid epoch header — skip to the next plausible one.
            if !line.trim().is_empty() {
                skipped_records += 1;
            }
            pos += 1;
            continue;
        };
        pos += 1;
        // Lines after the epoch header: continuation lines and/or records.
        let body_start = pos;

        // Epoch flags 2–5: event record.  The num_entries field is the count of
        // verbatim special-record lines that follow; there is NO satellite PRN
        // list on the header line.  Must check flag BEFORE calling collect_sats.
        if (2u8..=5).contains(&flag) {
            let epoch_idx = epochs.len();
            let end = (pos + num_entries).min(data_lines.len());
            if end - pos < num_entries {
                // Truncated by the end of the file.
                skipped_records += 1;
            }
            let event_lines: Vec<String> = data_lines[pos..end]
                .iter()
                .map(|l| (*l).to_string())
                .collect();
            pos = end;
            let changed = rinex_header::event_obs_codes(event_lines.iter().map(String::as_str))?;
            if let Some(codes) = changed.for_system(SatSystem::MIXED) {
                obs_codes = codes.to_vec();
                num_obs = obs_codes.len();
                obs_lines_per_sat = obs_lines_count(num_obs);
                masks.reset();
                series.reset();
            }
            events.push((
                epoch_idx,
                EventRecord {
                    flag,
                    epoch_line: line.to_owned(),
                    lines: event_lines,
                },
            ));
            continue;
        }

        // Flags 0, 1, 6: epoch header carries a satellite PRN list, continued
        // on following lines for more than 12 satellites.
        let mut sats: Vec<SatId> = Vec::with_capacity(num_entries);
        // 'line' still refers to the header line.
        let mut list_ok = collect_sats(line, num_entries.min(12), &mut sats).is_some();
        let mut cont_malformed = false;
        while list_ok && sats.len() < num_entries && pos < data_lines.len() {
            // Each continuation line must actually look like one.
            if !is_continuation_line(data_lines[pos]) {
                cont_malformed = true;
                break;
            }
            let want = (num_entries - sats.len()).min(12);
            list_ok = collect_sats(data_lines[pos], want, &mut sats).is_some();
            pos += 1;
        }

        // --- Validate epoch structure before committing ---
        //
        // A well-formed epoch has:
        //   1. Exactly num_entries valid satellite IDs, as rnx2srnx requires.
        //   2. No epoch-header line embedded in the observation block.
        //   3. The line after the observation block is an epoch header (or EOF).
        //
        // If any check fails, discard this epoch and scan forward to the
        // next epoch header.

        let obs_lines_total = num_entries * obs_lines_per_sat;

        // Check 1: a corrupt satellite ID (such as `e03` or `G1x`), a
        // continuation line that does not look like one, or a satellite list
        // longer or shorter than the claimed count.
        let sat_list_ok = list_ok && !cont_malformed && sats.len() == num_entries;

        // Check 2: is there an epoch or event header line within the expected
        // observation block?  If so, the epoch structure is broken — skip to
        // that header.  This also finds where to resume after a bad
        // satellite list, whose count may be wrong.
        let mut malformed = !sat_list_ok;
        let mut skip_to: Option<usize> = None;
        for i in pos..pos + obs_lines_total {
            if i >= data_lines.len() {
                malformed = true;
                break;
            }
            if is_any_header_line(data_lines[i]) {
                skip_to = Some(i);
                malformed = true;
                break;
            }
        }

        // Check 3: the line after the observation block should be an epoch
        // header, event record header, or EOF.  Non-conforming files
        // (e.g., teqc-spliced) may have extra lines (comments, blanks) between
        // epochs.  Scan forward to find the next epoch or event header; if
        // found nearby, skip the extra lines after parsing observations.
        // Only mark malformed if no header appears within a reasonable distance.
        let mut skip_extra_to: Option<usize> = None;
        if !malformed {
            let after_obs = pos + obs_lines_total;
            if after_obs < data_lines.len() && !is_any_header_line(data_lines[after_obs]) {
                // Scan forward for the next epoch or event header.
                let mut scan_pos = after_obs;
                let mut found_next = None;
                while scan_pos < data_lines.len() {
                    if is_any_header_line(data_lines[scan_pos]) {
                        found_next = Some(scan_pos);
                        break;
                    }
                    scan_pos += 1;
                    // Give up after scanning too far — something is genuinely wrong.
                    if scan_pos - after_obs > 100 {
                        break;
                    }
                }
                if let Some(next_header) = found_next {
                    // Extra lines between epochs — remember to skip them.
                    skip_extra_to = Some(next_header);
                } else {
                    // No header found — genuinely malformed.
                    malformed = true;
                }
            }
            // At EOF with exactly the right number of lines is fine.
        }

        if malformed {
            skipped_records += 1;
            // Already found an epoch header mid-block? Jump there.
            if let Some(target) = skip_to {
                pos = target;
                continue;
            }
            // Otherwise scan forward from after the expected block.
            let scan_from = pos + obs_lines_total.min(data_lines.len().saturating_sub(pos));
            pos = scan_from;
            while pos < data_lines.len() {
                if is_any_header_line(data_lines[pos]) {
                    break;
                }
                pos += 1;
            }
            continue;
        }

        // --- Valid epoch — parse observations ---

        // Flags 1 (power failure) and 6 (cycle slip) share the observation
        // layout, but SRNX keeps such records only in EVTF, never as an
        // observation epoch.  Keep them as an event (continuation and
        // observation lines), not in `epochs` / `observations`.
        if flag == 1 || flag == 6 {
            pos += obs_lines_total;
            events.push((
                epochs.len(),
                EventRecord {
                    flag,
                    epoch_line: line.to_owned(),
                    lines: data_lines[body_start..pos]
                        .iter()
                        .map(|l| (*l).to_string())
                        .collect(),
                },
            ));
            // Skip any extra lines between epochs.
            if let Some(target) = skip_extra_to {
                pos = target;
            }
            continue;
        }

        // Normal epoch (flag 0): parse observation records.
        //
        // A satellite listed twice in one epoch is malformed; keep its first
        // record, as the RINEX v3 reader and rnx2srnx do, and report the rest.
        seen_sats.clear();
        let duplicated: Vec<bool> = sats.iter().map(|&sat| !seen_sats.insert(sat)).collect();
        skipped_records += duplicated.iter().filter(|&&d| d).count();

        let epoch_idx = epochs.len();
        epochs.push(epoch_rec);

        for (&sat, &duplicate) in sats.iter().zip(&duplicated) {
            if duplicate {
                pos += obs_lines_per_sat;
                continue;
            }
            let keep = masks.get(sat, &obs_codes);
            let handle = series.satellite(sat);
            if keep.map_or(true, |keep| keep.contains(&true)) {
                series.list(handle, epoch_idx);
            }
            for (obs_line_idx, obs_chunk) in
                data_lines[pos..pos + obs_lines_per_sat].iter().enumerate()
            {
                let base_obs = obs_line_idx * 5;
                let end_obs = (base_obs + 5).min(num_obs);
                for (slot, &code) in obs_codes[base_obs..end_obs].iter().enumerate() {
                    if keep.is_some_and(|keep| !keep[base_obs + slot]) {
                        continue;
                    }
                    if let Some(obs) = parse_obs_field(obs_chunk, slot) {
                        series
                            .series(handle, base_obs + slot, code)
                            .push((epoch_idx, obs));
                    }
                }
            }
            pos += obs_lines_per_sat;
        }

        // Skip any extra lines between epochs (e.g., teqc splice comments).
        if let Some(target) = skip_extra_to {
            pos = target;
        }
    }

    let (observations, satellite_epochs) = series.finish();
    Ok(ObsFile {
        header: HeaderText::from_text(&raw_header_text),
        epochs,
        events,
        observations,
        satellite_epochs,
        skipped_records,
    })
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// The observation-code list.  For v2, all systems share a single list
/// stored under `SatSystem::MIXED`.  Fall back to the first list if MIXED is
/// absent (should not happen for conforming v2 files).
fn obs_codes(header: &Header) -> Result<&[ObsCode], Rinex2Error> {
    header
        .obs_codes
        .for_system(SatSystem::MIXED)
        .or_else(|| {
            header
                .obs_codes
                .by_system
                .first()
                .map(|(_, v)| v.as_slice())
        })
        .ok_or(Rinex2Error::NoObsCodes)
}

/// The lines of `text`, split as [`str::lines`] does but faster: a line
/// ends at LF, and a CR just before the LF is dropped.
pub(crate) fn split_lines(text: &str) -> Vec<&str> {
    let bytes = text.as_bytes();
    let mut lines = Vec::new();
    let mut start = 0;
    while start < bytes.len() {
        let Some(end) = find_newline(bytes, start) else {
            lines.push(&text[start..]);
            break;
        };
        let content = if end > start && bytes[end - 1] == b'\r' {
            end - 1
        } else {
            end
        };
        lines.push(&text[start..content]);
        start = end + 1;
    }
    lines
}

/// The index of the first LF in `bytes` at or after `from`, found eight
/// bytes at a time.
fn find_newline(bytes: &[u8], from: usize) -> Option<usize> {
    const ONES: u64 = u64::from_le_bytes([1; 8]);
    let mut at = from;
    while let Some(chunk) = bytes.get(at..at + 8) {
        // A zero byte of `x` sets its high bit in `found`; bits above the
        // first zero byte may be wrong, but the first is right.
        let x = u64::from_le_bytes(chunk.try_into().expect("8 bytes")) ^ (ONES * u64::from(b'\n'));
        let found = x.wrapping_sub(ONES) & !x & (ONES << 7);
        if found != 0 {
            return Some(at + found.trailing_zeros() as usize / 8);
        }
        at += 8;
    }
    bytes[at..].iter().position(|&b| b == b'\n').map(|i| at + i)
}

/// Number of 80-char observation lines required for `num_obs` obs codes.
#[inline]
fn obs_lines_count(num_obs: usize) -> usize {
    num_obs.div_ceil(5)
}

/// Check if a line could be an epoch header or event record header.
/// Uses a dummy line number since we only care about Ok/Err.
///
/// Fast-path: check the separator spaces at positions 26–27 and the epoch
/// flag digit at position 28 before calling the full parser.  Observation
/// lines have the hundredths/tenths of a seconds value at 26–27 (digits or
/// `.`), so a header line must have spaces there and a digit at 28.
fn is_any_header_line(line: &str) -> bool {
    let b = line.as_bytes();
    if b.len() < 32 || b[26] != b' ' || b[27] != b' ' || !b[28].is_ascii_digit() {
        return false;
    }
    parse_epoch_header(line, 0).is_ok()
}

/// Heuristic: does this line look like an epoch-header continuation?
///
/// Continuation lines start with 32 spaces then a letter (system id).
fn is_continuation_line(line: &str) -> bool {
    line.len() > 32
        && line.as_bytes()[..32].iter().all(|&b| b == b' ')
        && line
            .as_bytes()
            .get(32)
            .is_some_and(|&b| b.is_ascii_alphabetic())
}

/// Return a character at byte offset `i` in `line`, or `b' '` if the line is
/// shorter.  RINEX v2 lines may have trailing whitespace stripped by writers.
#[inline]
pub(crate) fn byte_at(line: &str, i: usize) -> u8 {
    line.as_bytes().get(i).copied().unwrap_or(b' ')
}

/// Parse a 2-digit year to a 4-digit year using the RINEX v2 convention:
/// 80–99 → 1980–1999, 00–79 → 2000–2079.
#[inline]
fn expand_year(y2: u8) -> u16 {
    if y2 >= 80 {
        1900 + u16::from(y2)
    } else {
        2000 + u16::from(y2)
    }
}

/// Parse a fixed-width decimal integer from `line[start..start+width]`,
/// interpreting blank/missing as zero.
fn parse_int_field(line: &str, start: usize, width: usize) -> i64 {
    let mut acc: i64 = 0;
    let mut any_digit = false;
    let mut negative = false;
    for i in start..start + width {
        let b = byte_at(line, i);
        match b {
            b'-' => negative = true,
            b'0'..=b'9' => {
                acc = acc * 10 + i64::from(b - b'0');
                any_digit = true;
            }
            _ => {}
        }
    }
    if !any_digit {
        return 0;
    }
    if negative {
        -acc
    } else {
        acc
    }
}

/// Parse a fixed-width floating-point field (already sliced out of the
/// line), interpreting blank/missing as zero.  Returns `None` if the field
/// is non-blank but not a finite number.
fn parse_float_field(field: &str) -> Option<f64> {
    let field = field.trim();
    if field.is_empty() {
        return Some(0.0);
    }
    field.parse::<f64>().ok().filter(|v| v.is_finite())
}

/// Parse the epoch / satellite header line.
///
/// Returns `(EpochRecord, flag, num_entries)` where `num_entries` is the
/// number of satellite IDs to read (or special-record lines for events).
fn parse_epoch_header(line: &str, line_no: usize) -> Result<(EpochRecord, u8, usize), Rinex2Error> {
    if line.len() < 32 {
        return Err(Rinex2Error::EpochLineTooShort(line_no));
    }

    // Byte layout (0-indexed):
    //  0     : space (1X)
    //  1-2   : year (I2.2)
    //  3     : space
    //  4-5   : month (I2)
    //  6     : space
    //  7-8   : day (I2)
    //  9     : space
    //  10-11 : hour (I2)
    //  12    : space
    //  13-14 : minute (I2)
    //  15-25 : seconds (F11.7)
    //  26-27 : spaces (2X)
    //  28    : epoch flag (I1)
    //  29-31 : num satellites (I3)
    //  32-67 : up to 12 satellite IDs  (12 × A1,I2)
    //  68-79 : receiver clock offset (F12.9, optional)

    // Verify separator positions are spaces.  This rejects observation data
    // lines that happen to have digits at timestamp positions.
    for &sep_pos in &[3, 6, 9, 12, 26, 27] {
        if byte_at(line, sep_pos) != b' ' {
            return Err(Rinex2Error::BadEpochHeader(line_no));
        }
    }

    let y2 = u8::try_from(parse_int_field(line, 1, 2))
        .map_err(|_| Rinex2Error::BadEpochHeader(line_no))?;
    let month = u8::try_from(parse_int_field(line, 4, 2))
        .map_err(|_| Rinex2Error::BadEpochHeader(line_no))?;
    let day = u8::try_from(parse_int_field(line, 7, 2))
        .map_err(|_| Rinex2Error::BadEpochHeader(line_no))?;
    let hour = u8::try_from(parse_int_field(line, 10, 2))
        .map_err(|_| Rinex2Error::BadEpochHeader(line_no))?;
    let minute = u8::try_from(parse_int_field(line, 13, 2))
        .map_err(|_| Rinex2Error::BadEpochHeader(line_no))?;

    // Seconds: F11.7 at positions 15-25.  Parse as f64 then convert to e7.
    let sec_str = line.get(15..line.len().min(26)).unwrap_or("");
    let sec_f64 = parse_float_field(sec_str)
        .filter(|&s| s >= 0.0)
        .ok_or(Rinex2Error::BadEpochHeader(line_no))?;
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let seconds_e7 = (sec_f64 * 1e7).round() as u64;

    let flag_raw = parse_int_field(line, 28, 1);
    if !(0..=6).contains(&flag_raw) {
        return Err(Rinex2Error::BadEpochHeader(line_no));
    }
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let flag = flag_raw as u8;

    // Validate timestamp ranges to reject observation data lines that
    // happen to have digits at epoch-header positions.  Events (flags 2-5)
    // may have missing timestamps in non-conforming files.
    let timestamp_present = month != 0 || day != 0 || hour != 0 || minute != 0 || sec_f64 != 0.0;
    if timestamp_present {
        if !(1..=12).contains(&month) || !(1..=31).contains(&day) || hour > 23 || minute > 59 {
            return Err(Rinex2Error::BadEpochHeader(line_no));
        }
    } else if !(2..=5).contains(&flag) {
        // Normal epochs must have a valid timestamp.
        return Err(Rinex2Error::BadEpochHeader(line_no));
    }

    let num_entries_raw = parse_int_field(line, 29, 3);
    if num_entries_raw < 0 {
        return Err(Rinex2Error::BadEpochHeader(line_no));
    }
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let num_entries = num_entries_raw as usize;

    // Receiver clock offset: F12.9 at positions 68-79.
    let clock_offset_e12 = if line.len() > 68 {
        let clk_str = line.get(68..line.len().min(80)).unwrap_or("");
        let clk_sec = parse_float_field(clk_str).ok_or(Rinex2Error::BadEpochHeader(line_no))?;
        let clk_ps = (clk_sec * 1e12).round();
        // Reject values (e.g. in exponent notation) that would saturate.
        #[allow(clippy::cast_precision_loss)]
        if !(i64::MIN as f64..i64::MAX as f64).contains(&clk_ps) {
            return Err(Rinex2Error::BadEpochHeader(line_no));
        }
        #[allow(clippy::cast_possible_truncation)]
        let raw = clk_ps as i64;
        raw
    } else {
        0
    };

    let rec = EpochRecord {
        year: expand_year(y2),
        month,
        day,
        hour,
        minute,
        seconds_e7,
        extra_ps: 0,
        flag,
        clock_offset_e12,
    };

    Ok((rec, flag, num_entries))
}

/// Append the first `want` satellite IDs in `line[32..68]` (12 × A1,I2) to
/// `sats`.  Epoch lines and continuation lines both put them there.
///
/// Returns `None`, leaving `sats` partly extended, if one of those slots is
/// not a valid satellite ID (see [`parse_sat_id`]) or a later slot is not
/// blank, as when the satellite count disagrees with the list.
fn collect_sats(line: &str, want: usize, sats: &mut Vec<SatId>) -> Option<()> {
    for slot in 0..12 {
        let start = 32 + slot * 3;
        let field = [start, start + 1, start + 2].map(|i| byte_at(line, i));
        if slot < want {
            sats.push(parse_sat_id(std::str::from_utf8(&field).ok()?)?);
        } else if field != *b"   " {
            return None;
        }
    }
    Some(())
}

/// Parse one observation field from `line` at the given slot (0–4).
///
/// Each slot is 16 bytes: 14-byte F14.3 value + 1-byte LLI + 1-byte SSI.
/// Returns `None` for an absent observation (see [`parse_obs_at`]).
fn parse_obs_field(line: &str, slot: usize) -> Option<Obs> {
    parse_obs_at(line, slot * 16)
}

/// Parse the 16-byte `F14.3, I1, I1` observation field starting at byte
/// offset `start`.  Shared with the RINEX v3/v4 reader, whose fields follow
/// a 3-byte satellite identifier.
///
/// Returns `None` if the value has no digits, dropping any LLI/SSI with it
/// as `rnx2srnx` does.  RINEX says "missing observations are written as 0.0
/// or blanks", but zero is often a real value (an LLI, or the Doppler of a
/// geostationary satellite), so SRNX keeps a written zero as zero.
///
/// F14.3 values have exactly 3 decimal places.  We parse them directly as
/// integers by skipping the decimal point — no `String` allocation and no
/// `f64` round-trip.
pub(crate) fn parse_obs_at(line: &str, start: usize) -> Option<Obs> {
    let bytes = line.as_bytes();
    if let Some(field) = bytes.get(start..start + 16) {
        parse_obs(field.try_into().expect("16 bytes"))
    } else {
        let mut field = [b' '; 16];
        let rest = bytes.get(start..).unwrap_or(&[]);
        for (to, &from) in field.iter_mut().zip(rest) {
            *to = from;
        }
        parse_obs(&field)
    }
}

const BLANKS: u64 = u64::from_le_bytes([b' '; 8]);

/// [`parse_obs_at`] for a field of exactly 16 bytes.
#[inline]
fn parse_obs(field: &[u8; 16]) -> Option<Obs> {
    let word = |at: usize| u64::from_le_bytes(field[at..at + 8].try_into().expect("8 bytes"));
    if word(0) == BLANKS && word(6) == BLANKS {
        return None;
    }
    let value = match parse_canonical(field) {
        Some(value) => value,
        None => parse_value(field)?,
    };
    Some(Obs {
        value,
        lli: field[14],
        ssi: field[15],
    })
}

/// The value of an `F14.3` field written as usual: blanks, an optional
/// `-`, digits, `.` and three more digits.  Returns `None` for any other
/// field, which [`parse_value`] handles, and gives the same value as it
/// for these.
#[inline]
fn parse_canonical(field: &[u8; 16]) -> Option<i64> {
    const HIGH_BITS: u64 = 0x8080_8080_8080_8080;
    if field[10] != b'.' {
        return None;
    }
    let frac = [field[11], field[12], field[13]].map(|b| b.wrapping_sub(b'0'));
    if frac.iter().any(|&d| d > 9) {
        return None;
    }
    let frac = i64::from(frac[0]) * 100 + i64::from(frac[1]) * 10 + i64::from(frac[2]);

    // Columns 3–10 of the integer part, eight at a time.  For ASCII bytes,
    // `b ^ b'0'` is below 10 just for digits, and adding 0x76 sets the high
    // bit of the others without carrying into the next byte.
    let word = u64::from_le_bytes(field[2..10].try_into().expect("8 bytes"));
    if word & HIGH_BITS != 0 {
        return None;
    }
    let values = word ^ u64::from_le_bytes([b'0'; 8]);
    let others = ((values + 0x7676_7676_7676_7676) & HIGH_BITS) >> 7;
    let others = others * 0xFF;
    // The digits must follow every other byte (the low bytes come first).
    if others & others.wrapping_add(1) != 0 {
        return None;
    }
    let low = i64::try_from(eight_digits(values & !others)).expect("eight digits");
    let leading = others.trailing_ones() / 8;

    // Columns 1–2, and the sign.
    let digit = |b: u8| b.is_ascii_digit().then(|| i64::from(b - b'0'));
    let mut negative = false;
    let mut high = 0;
    if leading == 0 {
        match (field[0], field[1]) {
            (b' ', b' ') => {}
            (b' ', b'-') => negative = true,
            (b' ', d) => high = digit(d)?,
            (b'-', d) => {
                negative = true;
                high = digit(d)?;
            }
            (d, e) => high = digit(d)? * 10 + digit(e)?,
        }
    } else {
        if field[0] != b' ' || field[1] != b' ' {
            return None;
        }
        let blanks = BLANKS & others;
        if word & others != blanks {
            let minus = blanks ^ (u64::from(b' ' ^ b'-') << (8 * (leading - 1)));
            if word & others != minus {
                return None;
            }
            negative = true;
        }
    }

    let acc = (high * 100_000_000 + low) * 1000 + frac;
    Some(if negative { -acc } else { acc })
}

/// The number written by eight decimal digit values, one per byte, the
/// most significant in the low byte.
///
/// Each multiplication overflows `u64` by design: the windows kept by the
/// masks never receive a carry, and the discarded high part holds only the
/// most significant byte/word of the product, so the low 64 bits are exact.
/// The multiplies are `wrapping_mul` to keep debug builds from panicking.
#[inline]
fn eight_digits(values: u64) -> u64 {
    let pairs = values.wrapping_mul(10 << 8 | 1) >> 8 & 0x00FF_00FF_00FF_00FF;
    let quads = pairs.wrapping_mul(100 << 16 | 1) >> 16 & 0x0000_FFFF_0000_FFFF;
    quads.wrapping_mul(10_000 << 32 | 1) >> 32
}

/// The value of any `F14.3` field: its digits as an integer (×1000
/// milliunits), negated if it has a `-`.  Returns `None` if the field has
/// no digits or a character other than a digit, blank, `.`, `+` or `-`.
fn parse_value(field: &[u8; 16]) -> Option<i64> {
    let mut acc: i64 = 0;
    let mut negative = false;
    let mut digits = false;
    for &b in &field[..14] {
        match b {
            b'-' => negative = true,
            b'0'..=b'9' => {
                acc = acc * 10 + i64::from(b - b'0');
                digits = true;
            }
            b'.' | b' ' | b'+' => {}
            _ => return None, // unexpected character → treat as absent
        }
    }
    digits.then_some(if negative { -acc } else { acc })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- unit tests -----------------------------------------------

    #[test]
    fn obs_lines_count_cases() {
        assert_eq!(obs_lines_count(0), 0);
        assert_eq!(obs_lines_count(1), 1);
        assert_eq!(obs_lines_count(5), 1);
        assert_eq!(obs_lines_count(6), 2);
        assert_eq!(obs_lines_count(10), 2);
        assert_eq!(obs_lines_count(11), 3);
        assert_eq!(obs_lines_count(20), 4);
        assert_eq!(obs_lines_count(21), 5);
    }

    #[test]
    fn lines_split_like_str_lines() {
        for text in [
            "",
            "\n",
            "a",
            "a\n",
            "a\r\n",
            "a\n\n",
            "\r\n\r\n",
            "one line that is longer than eight bytes\r\nand another\nlast",
            "x\ry\r\n\rz",
            "0123456\n89abcdef\n0123456789abcdef0\r\n",
        ] {
            assert_eq!(
                split_lines(text),
                text.lines().collect::<Vec<_>>(),
                "{text:?}"
            );
        }
    }

    #[test]
    fn expand_year_cases() {
        assert_eq!(expand_year(0), 2000);
        assert_eq!(expand_year(20), 2020);
        assert_eq!(expand_year(79), 2079);
        assert_eq!(expand_year(80), 1980);
        assert_eq!(expand_year(99), 1999);
    }

    #[test]
    fn parse_obs_field_present() {
        // " 128250890.432 5" = F14.3 value + LLI=' ' + SSI='5'
        let line = " 128250890.432 5  99935750.10842";
        let obs = parse_obs_field(line, 0).unwrap();
        assert_eq!(obs.value, 128_250_890_432);
        assert_eq!(obs.lli, b' ');
        assert_eq!(obs.ssi, b'5');

        let obs2 = parse_obs_field(line, 1).unwrap();
        assert_eq!(obs2.value, 99_935_750_108);
        assert_eq!(obs2.lli, b'4');
        assert_eq!(obs2.ssi, b'2');

        // Largest canonical field: ten integer digits, the first a digit
        // (some writers put a digit in the spec's sign column).
        let obs4 = parse_obs_field("9999999999.999  ", 0).unwrap();
        assert_eq!(obs4.value, 9_999_999_999_999);
    }

    #[test]
    fn parse_obs_field_blank() {
        // Entirely blank field = absent.
        let line = "                                ";
        assert!(parse_obs_field(line, 0).is_none());
        assert!(parse_obs_field(line, 1).is_none());
    }

    #[test]
    fn parse_obs_field_zero_is_present() {
        // A written zero is a value, with its LLI/SSI; a blank value is
        // absent even with flags.
        let line = [
            "         0.000 7",
            "         -.0000 ",
            "             .  ",
            "              81",
        ]
        .concat();
        let line = line.as_str();
        let zero = |lli, ssi| Some(Obs { value: 0, lli, ssi });
        assert_eq!(parse_obs_field(line, 0), zero(b' ', b'7'));
        assert_eq!(parse_obs_field(line, 1), zero(b'0', b' '));
        assert!(parse_obs_field(line, 2).is_none());
        assert!(parse_obs_field(line, 3).is_none());
    }

    #[test]
    fn parse_obs_field_short_line() {
        // Line shorter than 80 chars: missing fields treated as blank.
        let line = "        35.250";
        let obs = parse_obs_field(line, 0).unwrap();
        assert_eq!(obs.value, 35_250);
        // LLI and SSI beyond end of line → spaces.
        assert_eq!(obs.lli, b' ');
        assert_eq!(obs.ssi, b' ');
        // Second slot is entirely beyond the line.
        assert!(parse_obs_field(line, 1).is_none());
    }

    /// The fast path agrees with the general parser on every field it
    /// takes, and takes the usual layouts.
    #[test]
    fn canonical_values_match_the_general_parser() {
        let mut taken = 0;
        let mut state = 0x9E37_79B9_7F4A_7C15_u64;
        let alphabet = b"0123456789     -.+x";
        for i in 0..400_000 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let mut field = [b' '; 16];
            if i % 2 == 0 {
                // A usual layout: a number of up to 14 places.
                let digits = 1 + usize::try_from(state % 10).unwrap();
                let mut n = state >> 8;
                for at in (0..14).rev() {
                    if at == 10 {
                        field[at] = b'.';
                    } else if at > 10 || 10 - at <= digits {
                        field[at] = b'0' + u8::try_from(n % 10).unwrap();
                        n /= 10;
                    }
                }
                let first = field.iter().position(|&b| b != b' ').unwrap();
                if state & 0x10 != 0 && first > 0 {
                    field[first - 1] = b'-';
                }
            } else {
                for (at, b) in field[..14].iter_mut().enumerate() {
                    let pick = (state >> (at * 4)) % alphabet.len() as u64;
                    *b = alphabet[usize::try_from(pick).unwrap()];
                }
            }
            if let Some(value) = parse_canonical(&field) {
                taken += 1;
                assert_eq!(
                    Some(value),
                    parse_value(&field),
                    "{}",
                    String::from_utf8_lossy(&field)
                );
            } else {
                assert!(i % 2 == 1, "{}", String::from_utf8_lossy(&field));
            }
        }
        assert!(taken >= 200_000, "{taken}");
    }

    #[test]
    fn parse_obs_field_negative() {
        // Doppler: -1010.755.  F14.3 = "     -1010.755" (5 spaces + 9 chars = 14).
        // LLI = ' ', SSI = ' '.  Total field = 16 chars.
        let line = "     -1010.755  ";
        let obs = parse_obs_field(line, 0).unwrap();
        assert_eq!(obs.value, -1_010_755);
        assert_eq!(obs.lli, b' ');
        assert_eq!(obs.ssi, b' ');
    }

    /// 12-satellite epoch header from at012000, up to (not including) the
    /// receiver clock offset field at positions 68-79.
    const EPOCH_HEADER: &str =
        " 20  7 18  0  0 30.0000000  0 29G13R10E12E01R03E26G16G11R02G08E31G27";

    #[test]
    fn parse_epoch_header_seconds_and_clock() {
        let line = format!("{EPOCH_HEADER}-0.123456789");
        let (rec, flag, num_entries) = parse_epoch_header(&line, 7).unwrap();
        assert_eq!(rec.seconds_e7, 300_000_000);
        assert_eq!(rec.clock_offset_e12, -123_456_789_000);
        assert_eq!(flag, 0);
        assert_eq!(num_entries, 29);

        // Blank clock field: no offset.
        let line = format!("{EPOCH_HEADER}            ");
        let (rec, _, _) = parse_epoch_header(&line, 7).unwrap();
        assert_eq!(rec.clock_offset_e12, 0);
    }

    #[test]
    fn parse_epoch_header_garbled_seconds() {
        for garbled in ["30.00x0000", "30.0.00000", "        NaN", "       -1.0"] {
            let line = format!(
                "{}{garbled:>11}{}",
                &EPOCH_HEADER[..15],
                &EPOCH_HEADER[26..]
            );
            assert_eq!(line.len(), EPOCH_HEADER.len(), "{garbled:?}");
            assert!(
                matches!(
                    parse_epoch_header(&line, 7),
                    Err(Rinex2Error::BadEpochHeader(7))
                ),
                "{garbled:?}"
            );
        }
    }

    #[test]
    fn parse_epoch_header_garbled_clock() {
        for garbled in ["-0.12x456789", "        inf", "     1e300"] {
            let line = format!("{EPOCH_HEADER}{garbled:>12}");
            assert!(
                matches!(
                    parse_epoch_header(&line, 7),
                    Err(Rinex2Error::BadEpochHeader(7))
                ),
                "{garbled:?}"
            );
        }
    }

    #[test]
    fn collect_sats_first_line() {
        // First 12 sats from the at012000 epoch line (after pos 32).
        let line = " 20  7 18  0  0  0.0000000  0 29G13R10E12E01R03E26G16G11R02G08E31G27";
        let mut sats = Vec::new();
        assert_eq!(collect_sats(line, 12, &mut sats), Some(()));
        assert_eq!(sats.len(), 12);
        assert_eq!(sats[0], SatId::new(SatSystem::GPS, 13));
        assert_eq!(sats[1], SatId::new(SatSystem::GLONASS, 10));
        assert_eq!(sats[2], SatId::new(SatSystem::GALILEO, 12));
        assert_eq!(sats[11], SatId::new(SatSystem::GPS, 27));
    }

    #[test]
    fn collect_sats_continuation() {
        let line = "                                G20G30E13E24R16G10E07R09G21R19R18E03";
        let mut sats = Vec::new();
        assert_eq!(collect_sats(line, 12, &mut sats), Some(()));
        assert_eq!(sats.len(), 12);
        assert_eq!(sats[0], SatId::new(SatSystem::GPS, 20));
        assert_eq!(sats[11], SatId::new(SatSystem::GALILEO, 3));
    }

    #[test]
    fn collect_sats_rejects_corrupt_and_extra_ids() {
        let line = |ids: &str| format!(" 24 10  1  0  0  0.0000000  0  2{ids}");
        let mut sats = Vec::new();
        // A blank system letter means GPS, and a blank tens digit is zero.
        assert_eq!(collect_sats(&line("G01  5"), 2, &mut sats), Some(()));
        let gps = |prn| SatId::new(SatSystem::GPS, prn);
        assert_eq!(sats, [gps(1), gps(5)]);
        // Lower-case system, letter or blank in the PRN, PRN 0, too few IDs,
        // and too many IDs.
        for bad in ["G01e02", "G01G1x", "G01G  ", "G01G00", "G01", "G01G02G03"] {
            assert_eq!(
                collect_sats(&line(bad), 2, &mut Vec::new()),
                None,
                "{bad:?}"
            );
        }
    }

    /// A one-satellite-system, one-observation file with epochs at 0, 30 and
    /// 60 s, the middle one given by `middle`.
    fn three_epochs(middle: &str) -> String {
        let h = |content: &str, label: &str| format!("{content:<60}{label}\n");
        let epoch = |sec: u32, n: usize, ids: &str| {
            format!(
                " 24 10  1  0  0{:>11}  0{n:>3}{ids}\n",
                format!("{sec}.0000000")
            )
        };
        let mut text = h(
            "     2.11           OBSERVATION DATA    G (GPS)",
            "RINEX VERSION / TYPE",
        );
        text += &h("     1    C1", "# / TYPES OF OBSERV");
        text += &h("", "END OF HEADER");
        text += &epoch(0, 2, "G01G02");
        text += "  20000000.000\n  20000001.000\n";
        text += &middle.replace("EPOCH", &epoch(30, 0, "")[..29]);
        text += &epoch(60, 2, "G01G02");
        text += "  20000060.000\n  20000061.000\n";
        text
    }

    #[test]
    fn epochs_with_bad_satellite_lists_are_dropped() {
        let obs = "  20000030.000\n  20000031.000\n";
        let good = three_epochs(&format!("EPOCH  2G01G02\n{obs}"));
        assert_eq!(parse(good.as_bytes(), None).unwrap().epochs.len(), 3);
        assert_eq!(
            parse(good.as_bytes(), None).unwrap().skipped_records,
            0,
            "a well-formed file skips nothing"
        );

        for middle in [
            "EPOCH  2G01e02\n", // lower-case system letter
            "EPOCH  2G01G  \n", // blank PRN
            "EPOCH  2G01G0x\n", // letter in the PRN
            "EPOCH  2G01G00\n", // PRN 0
            "EPOCH  3G01G02\n", // count larger than the list
            "EPOCH  1G01G02\n", // count smaller than the list
        ] {
            let text = three_epochs(&format!("{middle}{obs}"));
            let f = parse(text.as_bytes(), None).unwrap();
            let secs: Vec<u64> = f.epochs.iter().map(|e| e.seconds_e7).collect();
            assert_eq!(secs, [0, 600_000_000], "{middle:?}");
            assert!(f.skipped_records > 0, "{middle:?}");
            let c1 = &f.observations[&(
                SatId::new(SatSystem::GPS, 1),
                ObsCode::from_v2("C1").unwrap(),
            )];
            let got: Vec<(usize, i64)> = c1.iter().map(|&(i, o)| (i, o.value)).collect();
            assert_eq!(
                got,
                [(0, 20_000_000_000), (1, 20_000_060_000)],
                "{middle:?}"
            );
        }
    }

    #[test]
    fn a_truncated_event_is_counted() {
        let h = |content: &str, label: &str| format!("{content:<60}{label}\n");
        let text = h(
            "     2.11           OBSERVATION DATA    G (GPS)",
            "RINEX VERSION / TYPE",
        ) + &h("     1    C1", "# / TYPES OF OBSERV")
            + &h("", "END OF HEADER")
            + " 24 10  1  0  0  0.0000000  0  1G01\n  20000000.000\n"
            + " 24 10  1  0  0  1.0000000  4  2\n"
            + &h("ONLY ONE OF TWO", "COMMENT");
        let f = parse(text.as_bytes(), None).unwrap();
        assert_eq!(f.epochs.len(), 1);
        assert_eq!(f.events.len(), 1);
        assert_eq!(f.skipped_records, 1);
    }

    #[test]
    fn a_duplicate_satellite_keeps_its_first_record() {
        let h = |content: &str, label: &str| format!("{content:<60}{label}\n");
        let text = h(
            "     2.11           OBSERVATION DATA    G (GPS)",
            "RINEX VERSION / TYPE",
        ) + &h("     1    C1", "# / TYPES OF OBSERV")
            + &h("", "END OF HEADER")
            + " 24 10  1  0  0  0.0000000  0  2G01G01\n  20000000.000\n  20000001.000\n";
        let f = parse(text.as_bytes(), None).unwrap();
        assert_eq!(f.epochs.len(), 1);
        assert_eq!(f.skipped_records, 1);
        let key = (
            SatId::new(SatSystem::GPS, 1),
            ObsCode::from_v2("C1").unwrap(),
        );
        let got: Vec<(usize, i64)> = f.observations[&key]
            .iter()
            .map(|&(i, o)| (i, o.value))
            .collect();
        assert_eq!(got, [(0, 20_000_000_000)]);
    }

    #[test]
    fn flag1_and_flag6_records_are_events() {
        let text = include_str!("../tests/data/rinex211_flag1_flag6.rnx");
        let f = parse(text.as_bytes(), None).unwrap();

        // Only the two flag-0 epochs are observation epochs.
        let secs: Vec<u64> = f.epochs.iter().map(|e| e.seconds_e7).collect();
        assert_eq!(secs, vec![0, 20_000_000]);
        assert!(f.epochs.iter().all(|e| e.flag == 0));

        // Flag 1 (with a satellite continuation line) and flag 6 are kept as
        // events positioned between the two epochs.
        let events: Vec<(usize, u8, usize)> = f
            .events
            .iter()
            .map(|(i, e)| (*i, e.flag, e.lines.len()))
            .collect();
        assert_eq!(events, vec![(1, 1, 14), (1, 6, 1)]);
        assert_eq!(f.events[0].1.lines[0].trim(), "G13");

        // The flag-1 values are not observations; G03..G13 appear only there.
        let g01_c1 = &f.observations[&(
            SatId::new(SatSystem::GPS, 1),
            ObsCode::from_v2("C1").unwrap(),
        )];
        let got: Vec<(usize, i64)> = g01_c1.iter().map(|&(i, o)| (i, o.value)).collect();
        assert_eq!(got, vec![(0, 20_001_000_000), (1, 20_001_002_000)]);
        assert!(f.observations.keys().all(|(s, _)| s.prn <= 2));

        // Written zeros are observations, with their flags.
        let g02 = |code| {
            let key = (
                SatId::new(SatSystem::GPS, 2),
                ObsCode::from_v2(code).unwrap(),
            );
            f.observations[&key].last().copied()
        };
        let zero = |lli| {
            Some((
                1,
                Obs {
                    value: 0,
                    lli,
                    ssi: b' ',
                },
            ))
        };
        assert_eq!(g02("C1"), zero(b' '));
        assert_eq!(g02("L1"), zero(b'0'));
    }

    // ---- fixture tests -----------------------------------------------

    fn fixture(stem: &str) -> ObsFile {
        parse(&crate::test_rinex(stem), None).unwrap()
    }

    #[test]
    fn fixture_at01_header_and_epoch_count() {
        let obs = fixture("at01_v211_1hz");
        assert_eq!(obs.header.get_line("MARKER NAME"), Some("AT01".to_string()));
        // Two minutes at 1 s.
        assert_eq!(obs.epochs.len(), 120);
    }

    #[test]
    fn fixture_at01_first_epoch_timestamp() {
        let obs = fixture("at01_v211_1hz");
        let e0 = &obs.epochs[0];
        assert_eq!(e0.year, 2020);
        assert_eq!(e0.month, 7);
        assert_eq!(e0.day, 18);
        assert_eq!(e0.hour, 0);
        assert_eq!(e0.minute, 0);
        assert_eq!(e0.seconds_e7, 0);
        assert_eq!(e0.flag, 0);
    }

    #[test]
    fn fixture_at01_g13_l1_first_epoch() {
        let obs = fixture("at01_v211_1hz");
        // In the first epoch, G13 L1 = 128250890.432, SSI=5.
        let l1 = ObsCode::from_v2("L1").unwrap();
        let g13 = SatId::new(SatSystem::GPS, 13);
        let series = obs.observations.get(&(g13, l1)).expect("G13/L1 present");
        // Sparse series: first entry should be (epoch_idx=0, obs).
        let &(ep_idx, first) = series.first().expect("at least one G13/L1 observation");
        assert_eq!(ep_idx, 0, "G13/L1 observed at epoch 0");
        assert_eq!(first.value, 128_250_890_432);
        assert_eq!(first.lli, b' ');
        assert_eq!(first.ssi, b'5');
    }

    #[test]
    fn fixture_at01_s1_first_epoch() {
        let obs = fixture("at01_v211_1hz");
        // G13 S1 = 35.250 in epoch 0 (obs index 5 on obs line 2).
        let s1 = ObsCode::from_v2("S1").unwrap();
        let g13 = SatId::new(SatSystem::GPS, 13);
        let series = obs.observations.get(&(g13, s1)).expect("G13/S1 present");
        let &(ep_idx, first) = series.first().expect("at least one G13/S1 observation");
        assert_eq!(ep_idx, 0);
        assert_eq!(first.value, 35_250);
    }

    #[test]
    fn fixture_observation_series_sorted() {
        for stem in ["at01_v211_1hz", "n302_v211_30s", "mty2_v211_15s"] {
            let obs = fixture(stem);
            // The fixtures carry blank lines between epochs, which the
            // reader skips without counting them.
            assert_eq!(obs.skipped_records, 0, "{stem}");
            // Every observation series must have epoch indices in ascending
            // order, within the valid range, and contain at least one entry.
            let num_epochs = obs.epochs.len();
            for (key, series) in &obs.observations {
                let name = format!("{stem}: {}/{}", key.0, key.1.as_str());
                assert!(!series.is_empty(), "{name}: empty series");
                assert!(
                    series.windows(2).all(|w| w[0].0 < w[1].0),
                    "{name}: not sorted"
                );
                assert!(
                    series.last().map(|&(i, _)| i).unwrap() < num_epochs,
                    "{name}: epoch index out of range"
                );
            }
        }
    }

    #[test]
    fn fixture_mty2_event() {
        // 80 records at 15 s; record 40 is a flag-4 event, not an epoch.
        let obs = fixture("mty2_v211_15s");
        assert_eq!(obs.epochs.len(), 79);
        assert_eq!(obs.events.len(), 1);
        assert_eq!((obs.events[0].0, obs.events[0].1.flag), (40, 4));
    }

    /// A file with `before` observation types, a flag-4 event that changes
    /// them to `after`, and one epoch of `values` (one line of up to five).
    fn types_change(before: &[&str], after: &[&str], values: &str) -> String {
        let h = |content: &str, label: &str| format!("{content:<60}{label}\n");
        let types = |codes: &[&str]| {
            let list = codes
                .iter()
                .map(|c| format!("{c:>6}"))
                .collect::<Vec<_>>()
                .concat();
            h(&format!("{:>6}{list}", codes.len()), "# / TYPES OF OBSERV")
        };
        let mut text = h(
            "     2.11           OBSERVATION DATA    G (GPS)",
            "RINEX VERSION / TYPE",
        );
        text += &types(before);
        text += &h("", "END OF HEADER");
        text += " 24 10  1  0  0  0.0000000  0  1G01\n";
        for line in 0..before.len().div_ceil(5) {
            let n = (before.len() - 5 * line).min(5);
            text += &"  20001000.000 6".repeat(n);
            text += "\n";
        }
        text += " 24 10  1  0  0 30.0000000  4  1\n";
        text += &types(after);
        text += " 24 10  1  0  1  0.0000000  0  1G01\n";
        text += values;
        text += "\n";
        text
    }

    #[test]
    fn events_change_the_observation_types() {
        let g01 = SatId::new(SatSystem::GPS, 1);
        let value = |f: &ObsFile, code: &str| {
            f.observations[&(g01, ObsCode::from_v2(code).unwrap())]
                .iter()
                .find(|&&(i, _)| i == 1)
                .map(|&(_, o)| o.value)
        };
        // Fewer lines per satellite after the change.
        let seven = ["C1", "L1", "L2", "P2", "D1", "S1", "S2"];
        let text = types_change(&seven, &["C1", "L1"], "  11111111.000 7  22222222.000 7");
        let f = parse(text.as_bytes(), None).unwrap();
        assert_eq!(f.epochs.len(), 2);
        assert_eq!(value(&f, "C1"), Some(11_111_111_000));
        assert_eq!(value(&f, "L1"), Some(22_222_222_000));
        assert_eq!(value(&f, "L2"), None);
        // Reordered: values follow the new order.
        let text = types_change(
            &["C1", "L1"],
            &["L1", "C1"],
            "  11111111.000 7  22222222.000 7",
        );
        let f = parse(text.as_bytes(), None).unwrap();
        assert_eq!(value(&f, "L1"), Some(11_111_111_000));
        assert_eq!(value(&f, "C1"), Some(22_222_222_000));

        // A filter is asked once per pair, and applies to the new order.
        let mut asked = Vec::new();
        let mut keep_c1 = |_: SatId, code: ObsCode| {
            asked.push(code);
            code.as_str() == "C1"
        };
        let f = parse(text.as_bytes(), Some(&mut keep_c1)).unwrap();
        assert_eq!(value(&f, "C1"), Some(22_222_222_000));
        assert!(!f
            .observations
            .contains_key(&(g01, ObsCode::from_v2("L1").unwrap())));
        assert_eq!(asked.len(), 2);
    }
}
