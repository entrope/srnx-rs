//! Parser for RINEX v3.xx and v4.xx observation data files.
//!
//! Reads the text format defined in `docs/rinex_4.02.pdf`, Table A3.  The
//! observation data records are unchanged from RINEX 3.00 through 4.02
//! except for the optional picosecond extension of the epoch record added in
//! 4.02.  Populates the shared [`ObsFile`] semantic model.
//!
//! ## Format summary
//!
//! ### Epoch record
//!
//! ```text
//! > YYYY MM DD HH MM SS.SSSSSSS  F NNN      CCCCCCCCCCCCCCC PPPPP
//! ```
//!
//! Fortran: `A1, 1X,I4, 4(1X,I2.2), F11.7, 2X,I1, I3, 6X, F15.12, 1X,I5.5`
//!
//! The receiver clock offset (`F15.12`) is optional.  The trailing `I5.5`
//! field (RINEX 4.02+) holds five more digits of seconds, extending the
//! resolution to picoseconds; when it is present the `F11.7` seconds field
//! is floored rather than rounded.
//!
//! ### Observation records (epoch flags 0, 1 and 6)
//!
//! One line per satellite: `A1,I2.2, m(F14.3, I1, I1)`, where `m` is the
//! number of observation codes declared for the satellite's system in the
//! `SYS / # / OBS TYPES` header records.  Unlike RINEX v2 there is no
//! 80-column wrapping, and trailing blanks may be truncated.
//!
//! Only flag-0 epochs become [`ObsFile::epochs`].  Flag-1 (power failure)
//! and flag-6 (cycle slip) records are kept as events, matching the SRNX
//! rule that such records live only in `EVTF`.
//!
//! ### Event records (epoch flags 2–5)
//!
//! The "number of satellites" field instead gives the number of verbatim
//! lines that follow.  The epoch fields may be blank when the event has no
//! significant epoch.
//!
//! Every record that starts a new epoch begins with `>`, so a reader can
//! always resynchronize after a malformed or truncated epoch by scanning for
//! the next such line.

use std::collections::HashSet;

use thiserror::Error;

use crate::model::{
    Catalog, EpochRecord, EventRecord, ObsFile, SatId, SeriesBuilder, SignalFilter, SignalMasks,
};
use crate::rinex2_obs::{byte_at, parse_obs_at, split_lines};
use crate::rinex_header::{self, Header, HeaderError, HeaderText, ObsCode};

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Rinex3Error {
    #[error("header parse error: {0}")]
    Header(#[from] HeaderError),
    #[error("no END OF HEADER line found")]
    NoEndOfHeader,
    #[error("RINEX version {major}.{minor:02} is not 3.xx or 4.xx")]
    UnsupportedVersion { major: u8, minor: u8 },
    #[error("observation file has no SYS / # / OBS TYPES records")]
    NoObsCodes,
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// List the satellites and signals declared by a RINEX v3 or v4
/// observation file's header.  Every satellite may carry every observation
/// code of its system.
pub fn catalog(input: &[u8]) -> Result<Catalog, Rinex3Error> {
    let header = rinex_header::parse(input)?;
    let codes_by_sys = codes_by_system(&header)?;
    Ok(Catalog::from_header(&header, |sys| {
        codes_by_sys[usize::from(sys.0)]
    }))
}

/// Parse a complete RINEX v3 or v4 observation file, keeping only the
/// satellite/signal pairs that `filter` selects (all of them for `None`).
///
/// Both LF and CRLF line endings are accepted.  Malformed epoch records are
/// skipped and counted in [`ObsFile::skipped_records`]; parsing resumes at
/// the next line starting with `>`.
#[allow(clippy::too_many_lines)]
pub fn parse(input: &[u8], filter: Option<SignalFilter<'_>>) -> Result<ObsFile, Rinex3Error> {
    let text = rinex_header::to_text(input);
    let lines = split_lines(&text);

    let eoh_idx = lines
        .iter()
        .position(|l| {
            let label = if l.len() >= 60 { l.get(60..) } else { Some(*l) };
            label.is_some_and(|s| s.trim().starts_with("END OF HEADER"))
        })
        .ok_or(Rinex3Error::NoEndOfHeader)?;

    // Raw header text through END OF HEADER, LF-normalized as in RHDR.
    let mut raw_header_text = lines[..=eoh_idx].join("\n");
    raw_header_text.push('\n');

    let header = rinex_header::parse_str(&text)?;
    // An event's `SYS / # / OBS TYPES` records change these.
    let mut codes_by_sys: Vec<Vec<ObsCode>> = codes_by_system(&header)?
        .iter()
        .map(|c| c.to_vec())
        .collect();
    let mut masks = SignalMasks::new(filter);
    for &sat in &header.satellites {
        masks.get(sat, &codes_by_sys[usize::from(sat.system.0)]);
    }

    let data = &lines[eoh_idx + 1..];
    let mut series = SeriesBuilder::default();
    let mut epochs: Vec<EpochRecord> = Vec::new();
    let mut events: Vec<(usize, EventRecord)> = Vec::new();
    let mut skipped_records = 0_usize;
    let mut seen_sats: HashSet<SatId> = HashSet::new();

    let mut pos = 0;
    while pos < data.len() {
        let line = data[pos];
        pos += 1;
        // Stray or malformed lines: skip until the next valid epoch record.
        let Some(epoch) = parse_epoch_line(line) else {
            if !line.trim().is_empty() {
                skipped_records += 1;
            }
            continue;
        };

        match epoch.record.flag {
            0 => {
                seen_sats.clear();
                let epoch_idx = epochs.len();
                epochs.push(epoch.record);
                for _ in 0..epoch.count {
                    // A `>` line inside the block, or the end of the file,
                    // means this epoch was truncated; keep what we have,
                    // report it, and resume there.
                    let Some(&sat_line) = data.get(pos).filter(|l| !l.starts_with('>')) else {
                        skipped_records += 1;
                        break;
                    };
                    pos += 1;
                    let Some(sat) = parse_sat_id(sat_line) else {
                        skipped_records += 1;
                        continue;
                    };
                    // A satellite listed twice in one epoch is malformed; keep
                    // its first record, as the RINEX v2 reader does.
                    if !seen_sats.insert(sat) {
                        skipped_records += 1;
                        continue;
                    }
                    let codes = &codes_by_sys[usize::from(sat.system.0)];
                    let keep = masks.get(sat, codes);
                    let handle = series.satellite(sat);
                    if keep.map_or(true, |keep| keep.contains(&true)) {
                        series.list(handle, epoch_idx);
                    }
                    for (slot, &code) in codes.iter().enumerate() {
                        if keep.is_some_and(|keep| !keep[slot]) {
                            continue;
                        }
                        let Some(obs) = parse_obs_at(sat_line, 3 + slot * 16) else {
                            continue;
                        };
                        series.series(handle, slot, code).push((epoch_idx, obs));
                    }
                }
            }
            flag => {
                // Flags 2–5: `count` special records follow verbatim.
                // Flags 1 and 6: `count` satellite records follow; like SRNX
                // (which stores them only in EVTF), keep them as an event
                // rather than as an observation epoch.
                let start = pos;
                if flag == 1 || flag == 6 {
                    while pos < data.len() && pos - start < epoch.count {
                        if data[pos].starts_with('>') {
                            break;
                        }
                        pos += 1;
                    }
                } else {
                    pos = (pos + epoch.count).min(data.len());
                    if apply_types_change(&data[start..pos], &mut codes_by_sys)? {
                        masks.reset();
                        series.reset();
                    }
                }
                if pos - start < epoch.count {
                    // Truncated by the end of the file or the next epoch.
                    skipped_records += 1;
                }
                let lines = data[start..pos].iter().map(|l| (*l).to_string()).collect();
                events.push((
                    epochs.len(),
                    EventRecord {
                        flag,
                        epoch_line: line.to_owned(),
                        lines,
                    },
                ));
            }
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

/// Apply the `SYS / # / OBS TYPES` records of an event's `lines` to the
/// code lists of later epochs.  Returns whether any list changed.
fn apply_types_change(
    lines: &[&str],
    codes_by_sys: &mut [Vec<ObsCode>],
) -> Result<bool, Rinex3Error> {
    let mut changed = false;
    for (sys, codes) in rinex_header::event_obs_codes(lines.iter().copied())?.by_system {
        if let Some(slot) = codes_by_sys.get_mut(usize::from(sys.0)) {
            *slot = codes;
            changed = true;
        }
    }
    Ok(changed)
}

/// Check the header's version and index its observation-code lists by
/// system letter.  The first declaration wins, matching
/// `ObsCodes::for_system` (used by the SRNX reader).
fn codes_by_system(header: &Header) -> Result<[&[ObsCode]; 128], Rinex3Error> {
    let version = header.version;
    if !matches!(version.major, 3 | 4) {
        return Err(Rinex3Error::UnsupportedVersion {
            major: version.major,
            minor: version.minor,
        });
    }
    let mut codes_by_sys: [&[ObsCode]; 128] = [&[]; 128];
    for (sys, codes) in &header.obs_codes.by_system {
        if let Some(slot) = codes_by_sys.get_mut(usize::from(sys.0)) {
            if slot.is_empty() {
                *slot = codes;
            }
        }
    }
    if codes_by_sys.iter().all(|c| c.is_empty()) {
        return Err(Rinex3Error::NoObsCodes);
    }
    Ok(codes_by_sys)
}

/// A decoded epoch record line.
#[derive(Debug)]
struct EpochLine {
    record: EpochRecord,
    /// Number of satellite (flags 0, 1, 6) or special-event (flags 2–5)
    /// records that follow.
    count: usize,
}

/// Parse an epoch record line.  Returns `None` if the line does not start
/// with `>` or is malformed.
///
/// Byte layout (0-indexed):
///
/// ```text
///  0     : '>'
///  2-5   : year (I4)
///  7-8   : month (I2.2)
///  10-11 : day
///  13-14 : hour
///  16-17 : minute
///  18-28 : seconds (F11.7)
///  31    : epoch flag (I1)
///  32-34 : number of satellites / special records (I3)
///  35-40 : reserved
///  41-55 : receiver clock offset (F15.12, optional)
///  57-61 : picoseconds beyond 1e-7 s (I5.5, optional, RINEX 4.02+)
/// ```
fn parse_epoch_line(line: &str) -> Option<EpochLine> {
    if !line.starts_with('>') {
        return None;
    }
    if byte_at(line, 29) != b' ' || byte_at(line, 30) != b' ' {
        return None;
    }
    let flag = match byte_at(line, 31) {
        d @ b'0'..=b'6' => d - b'0',
        _ => return None,
    };
    let count = usize::try_from(parse_uint(line, 32, 3)?).ok()?;

    let blank_time = (1..29).all(|i| byte_at(line, i) == b' ');
    let mut record = EpochRecord {
        year: 0,
        month: 0,
        day: 0,
        hour: 0,
        minute: 0,
        seconds_e7: 0,
        extra_ps: 0,
        flag,
        clock_offset_e12: 0,
    };
    if blank_time {
        // Only events may omit the epoch.
        if !(2..=5).contains(&flag) {
            return None;
        }
    } else {
        if [1, 6, 9, 12, 15].iter().any(|&i| byte_at(line, i) != b' ') {
            return None;
        }
        record.year = u16::try_from(parse_uint(line, 2, 4)?).ok()?;
        record.month = u8::try_from(parse_uint(line, 7, 2)?).ok()?;
        record.day = u8::try_from(parse_uint(line, 10, 2)?).ok()?;
        record.hour = u8::try_from(parse_uint(line, 13, 2)?).ok()?;
        record.minute = u8::try_from(parse_uint(line, 16, 2)?).ok()?;
        record.seconds_e7 = u64::try_from(parse_fixed(line, 18, 11, 7)?).ok()?;
        // Seconds may reach 60.x during a leap second.
        if !(1..=12).contains(&record.month)
            || !(1..=31).contains(&record.day)
            || record.hour > 23
            || record.minute > 59
            || record.seconds_e7 >= 610_000_000
        {
            return None;
        }
    }

    record.clock_offset_e12 = parse_fixed(line, 41, 15, 12)?;

    if (57..62).any(|i| byte_at(line, i) != b' ') {
        if byte_at(line, 56) != b' ' || !(57..62).all(|i| byte_at(line, i).is_ascii_digit()) {
            return None;
        }
        record.extra_ps = u32::try_from(parse_uint(line, 57, 5)?).ok()?;
    }

    Some(EpochLine { record, count })
}

/// Parse an unsigned fixed-width integer field at `line[start..start+width]`.
/// Blanks are ignored (so an all-blank field is zero); any other non-digit
/// makes the field invalid.
fn parse_uint(line: &str, start: usize, width: usize) -> Option<u64> {
    let mut acc: u64 = 0;
    for i in start..start + width {
        match byte_at(line, i) {
            b' ' => {}
            d @ b'0'..=b'9' => acc = acc * 10 + u64::from(d - b'0'),
            _ => return None,
        }
    }
    Some(acc)
}

/// Parse a Fortran `Fw.d` field at `line[start..start+width]` exactly, as
/// an integer scaled by `10^decimals`.  An all-blank field is zero.  Fewer
/// fractional digits than `decimals` are accepted; more are rejected, as is
/// anything that is not `[sign] digits [. digits]` surrounded by blanks.
fn parse_fixed(line: &str, start: usize, width: usize, decimals: u32) -> Option<i64> {
    let field: Vec<u8> = (start..start + width).map(|i| byte_at(line, i)).collect();
    let Some(first) = field.iter().position(|&b| b != b' ') else {
        return Some(0);
    };
    let last = field.iter().rposition(|&b| b != b' ').unwrap_or(first);
    let trimmed = &field[first..=last];
    let (negative, body) = match trimmed.split_first() {
        Some((b'-', rest)) => (true, rest),
        Some((b'+', rest)) => (false, rest),
        _ => (false, trimmed),
    };
    let (int_part, frac_part) = match body.iter().position(|&b| b == b'.') {
        Some(dot) => (&body[..dot], &body[dot + 1..]),
        None => (body, &[][..]),
    };
    let frac_digits = u32::try_from(frac_part.len()).ok()?;
    if frac_digits > decimals || (int_part.is_empty() && frac_part.is_empty()) {
        return None;
    }
    let mut acc: i64 = 0;
    for &b in int_part.iter().chain(frac_part) {
        if !b.is_ascii_digit() {
            return None;
        }
        acc = acc.checked_mul(10)?.checked_add(i64::from(b - b'0'))?;
    }
    acc = acc.checked_mul(10_i64.checked_pow(decimals - frac_digits)?)?;
    Some(if negative { -acc } else { acc })
}

/// Parse the `A1,I2.2` satellite identifier at the start of a record, as
/// [`SatId::parse_name`] does for the three bytes.
pub(crate) fn parse_sat_id(line: &str) -> Option<SatId> {
    SatId::parse_name([byte_at(line, 0), byte_at(line, 1), byte_at(line, 2)])
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rinex_header::SatSystem;

    // ---- unit tests -----------------------------------------------

    #[test]
    fn epoch_line_basic() {
        let e = parse_epoch_line("> 2020 07 07 00 00  0.0000000  0 24").unwrap();
        let r = &e.record;
        assert_eq!(
            (r.year, r.month, r.day, r.hour, r.minute),
            (2020, 7, 7, 0, 0)
        );
        assert_eq!(r.seconds_e7, 0);
        assert_eq!(r.extra_ps, 0);
        assert_eq!(r.flag, 0);
        assert_eq!(r.clock_offset_e12, 0);
        assert_eq!(e.count, 24);
    }

    #[test]
    fn epoch_line_clock_offset() {
        let e =
            parse_epoch_line("> 2006 03 24 13 10 54.0000000  0  7      -0.123456789210").unwrap();
        assert_eq!(e.record.seconds_e7, 540_000_000);
        assert_eq!(e.record.clock_offset_e12, -123_456_789_210);
        assert_eq!(e.count, 7);

        // Explicit zero written without a leading digit.
        let e =
            parse_epoch_line("> 2026 10 07 23 45  0.0000000  0 45        .000000000000").unwrap();
        assert_eq!(e.record.clock_offset_e12, 0);
    }

    #[test]
    fn epoch_line_picoseconds() {
        let e = parse_epoch_line("> 2006 03 24 13 10 36.1234567  0  5      -0.123456789012 12345")
            .unwrap();
        assert_eq!(e.record.seconds_e7, 361_234_567);
        assert_eq!(e.record.extra_ps, 12_345);
        assert_eq!(e.record.clock_offset_e12, -123_456_789_012);

        // Picoseconds without a clock offset.
        let e = parse_epoch_line("> 2006 03 24 13 10 36.0000000  0  5                      00042")
            .unwrap();
        assert_eq!(e.record.extra_ps, 42);
        assert_eq!(e.record.clock_offset_e12, 0);

        // Must be exactly five digits.
        assert!(
            parse_epoch_line("> 2006 03 24 13 10 36.0000000  0  5      -0.123456789012   345",)
                .is_none()
        );
    }

    #[test]
    fn epoch_line_blank_month_padding() {
        // Table A4 writes one month as " 3" rather than "03".
        let e =
            parse_epoch_line("> 2006  3 24 13 11 12.0000000  0  4      -0.123456789876").unwrap();
        assert_eq!(e.record.month, 3);
    }

    #[test]
    fn epoch_line_leap_second() {
        let e = parse_epoch_line("> 2016 12 31 23 59 60.5000000  0  1").unwrap();
        assert_eq!(e.record.seconds_e7, 605_000_000);
        assert!(parse_epoch_line("> 2016 12 31 23 59 61.0000000  0  1").is_none());
    }

    #[test]
    fn epoch_line_events() {
        let e = parse_epoch_line(">                              3  4").unwrap();
        assert_eq!(e.record.flag, 3);
        assert_eq!(e.count, 4);
        assert_eq!(e.record.year, 0);

        let e = parse_epoch_line("> 2006 03 24 13 13  1.2345678  5  0").unwrap();
        assert_eq!(e.record.flag, 5);
        assert_eq!(e.record.seconds_e7, 12_345_678);
        assert_eq!(e.count, 0);

        // A blank epoch is only legal for events.
        assert!(parse_epoch_line(">                              0  4").is_none());
    }

    #[test]
    fn epoch_line_rejects_garbage() {
        assert!(parse_epoch_line("G06  23629347.915            .300 8").is_none());
        assert!(parse_epoch_line("").is_none());
        assert!(parse_epoch_line(">").is_none());
        assert!(parse_epoch_line("> 2006 13 24 13 10 36.0000000  0  5").is_none());
        assert!(parse_epoch_line("> 2006 03 24 13 10 36.0000000  7  5").is_none());
        assert!(
            parse_epoch_line("> 2006 03 24 13 10 36.0000000  0  5      -0.12345678901x").is_none()
        );
        // Multi-byte UTF-8 must not cause a slicing panic.
        assert!(parse_epoch_line("> 2006 03 24 13 10 36.00000°0  0  5").is_none());
    }

    #[test]
    fn fixed_point_fields() {
        assert_eq!(parse_fixed("  36.0000000", 1, 11, 7), Some(360_000_000));
        assert_eq!(parse_fixed("  36.5", 0, 6, 7), Some(365_000_000));
        assert_eq!(parse_fixed("      ", 0, 6, 7), Some(0));
        assert_eq!(parse_fixed("-.5", 0, 3, 3), Some(-500));
        assert_eq!(parse_fixed("1.2345", 0, 6, 3), None);
        assert_eq!(parse_fixed("1 2", 0, 3, 3), None);
        assert_eq!(parse_fixed("-", 0, 1, 3), None);
        assert_eq!(parse_fixed("99999999999999.", 0, 15, 12), None);
    }

    #[test]
    fn sat_ids() {
        assert_eq!(
            parse_sat_id("G06  1.0"),
            Some(SatId::new(SatSystem::GPS, 6))
        );
        assert_eq!(parse_sat_id("S20"), Some(SatId::new(SatSystem::SBAS, 20)));
        assert_eq!(parse_sat_id("E 5"), Some(SatId::new(SatSystem::GALILEO, 5)));
        assert_eq!(parse_sat_id(" 12"), Some(SatId::new(SatSystem::GPS, 12)));
        assert_eq!(parse_sat_id("G0"), None);
        assert_eq!(parse_sat_id("g06"), None);
        assert_eq!(parse_sat_id("G00"), None);
        assert_eq!(parse_sat_id("G  "), None);
        assert_eq!(parse_sat_id("G1x"), None);
    }

    // ---- synthetic file tests -----------------------------------------

    /// RINEX 4.02 Table A4 example, with its satellite counts made
    /// consistent, plus a flag-6 epoch, a flag-1 epoch, a truncated epoch, a
    /// zero value, a record for an undeclared system and a satellite with no
    /// observations.
    const TABLE_A4: &str = include_str!("../tests/data/rinex402_table_a4_variants.rnx");

    /// Conforming RINEX 4.02 file using the picosecond epoch extension.
    const PICOSECONDS: &str = include_str!("../tests/data/rinex402_picoseconds.rnx");

    fn obs(file: &ObsFile, sat: &str, code: &str) -> Vec<(usize, i64, u8, u8)> {
        let sat = parse_sat_id(sat).unwrap();
        file.observations
            .get(&(sat, ObsCode::from_v3(code).unwrap()))
            .map(|s| s.iter().map(|&(i, o)| (i, o.value, o.lli, o.ssi)).collect())
            .unwrap_or_default()
    }

    #[test]
    fn table_a4_epochs_and_events() {
        let f = parse(TABLE_A4.as_bytes(), None).unwrap();
        assert_eq!(f.header.get_line("MARKER NAME"), Some("A 9080".to_string()));
        assert!(f.header.ends_with("END OF HEADER\n"));
        // The 13:14:48 epoch lists two of its three satellites.
        assert_eq!(f.skipped_records, 1);

        // Flag 0 only; flag-1 and flag-6 records become events.
        assert_eq!(f.epochs.len(), 5);
        assert!(f.epochs.iter().all(|e| e.flag == 0));
        let e0 = &f.epochs[0];
        assert_eq!((e0.hour, e0.minute, e0.seconds_e7), (13, 10, 360_000_000));
        assert_eq!(e0.extra_ps, 12_345);
        assert_eq!(e0.clock_offset_e12, -123_456_789_012);
        assert_eq!(f.epochs[1].extra_ps, 0);

        let flags: Vec<(usize, u8, usize)> = f
            .events
            .iter()
            .map(|(i, e)| (*i, e.flag, e.lines.len()))
            .collect();
        assert_eq!(
            flags,
            vec![
                (2, 2, 2),
                (3, 3, 4),
                (3, 5, 0),
                (3, 4, 2),
                (3, 6, 1),
                (3, 1, 2)
            ]
        );
        assert!(f.events[1].1.lines[0].starts_with("A 9081"));
        assert!(f.events[5].1.lines[1].starts_with("G09"));
    }

    #[test]
    fn table_a4_observations() {
        let f = parse(TABLE_A4.as_bytes(), None).unwrap();
        assert_eq!(
            obs(&f, "G06", "C1C"),
            vec![
                (0, 23_629_347_915, b' ', b' '),
                (1, 23_619_095_450, b' ', b' '),
                (2, 21_110_991_756, b' ', b' '),
            ]
        );
        // LLI blank, SSI 8.
        assert_eq!(obs(&f, "G06", "L1W")[0], (0, 300, b' ', b'8'));
        assert_eq!(obs(&f, "G06", "S2W")[0], (0, 24_158, b' ', b' '));
        assert_eq!(
            obs(&f, "E11", "L5I"),
            vec![(0, 178, b' ', b'7'), (1, 48_861_586, b' ', b'7')]
        );
        assert_eq!(obs(&f, "R21", "L1C"), vec![(1, 12_345_567, b' ', b'5')]);
        assert_eq!(obs(&f, "G16", "L2W")[0], (2, 57_505_177, b' ', b'2'));

        // Truncated epoch 3 (declares 3 satellites, has 2).
        assert_eq!(
            obs(&f, "G09", "C1C").last(),
            Some(&(3, 23_507_272_372, b' ', b' '))
        );
        // A written zero is an observation, with its SSI.
        assert_eq!(obs(&f, "G09", "L1W").last(), Some(&(3, 0, b' ', b'7')));
        // Satellite line with no observations at all.
        assert!(obs(&f, "G12", "C1C").iter().all(|&(i, ..)| i < 4));
        // No obs types declared for BeiDou: the record is ignored.
        assert!(f
            .observations
            .keys()
            .all(|(s, _)| s.system != SatSystem::BEIDOU));
        // Flag-1 and flag-6 values are not observations.
        assert_eq!(obs(&f, "G06", "L1W").len(), 3);
    }

    #[test]
    fn picosecond_epochs() {
        let f = parse(PICOSECONDS.as_bytes(), None).unwrap();
        // (minute, seconds_e7, extra_ps, clock_offset_e12, satellites)
        let got: Vec<(u8, u64, u32, i64)> = f
            .epochs
            .iter()
            .map(|e| (e.minute, e.seconds_e7, e.extra_ps, e.clock_offset_e12))
            .collect();
        assert_eq!(
            got,
            vec![
                (0, 1, 1, 100_001),
                (0, 10_000_001, 23_456, 123_456),
                (0, 20_000_000, 99_999, 0),
                (0, 599_999_999, 99_999, -1),
                (1, 0, 0, 0),
                (1, 10_000_001, 0, 100_000),
                (1, 20_000_001, 7, 100_007),
            ]
        );
        let fmt: Vec<String> = f.epochs.iter().map(crate::diff::fmt_epoch).collect();
        assert_eq!(fmt[0], "2024-10-01 00:00:00.000000100001");
        assert_eq!(fmt[3], "2024-10-01 00:00:59.999999999999");
        assert_eq!(fmt[4], "2024-10-01 00:01:00.0000000");

        // The flag-5 event carries picoseconds too, but events keep only text.
        assert_eq!(f.events.len(), 1);
        assert_eq!((f.events[0].0, f.events[0].1.flag), (2, 5));

        let per_epoch: Vec<usize> = (0..f.epochs.len())
            .map(|i| {
                f.observations
                    .iter()
                    .filter(|((_, c), s)| c.as_str() == "C1C" && s.iter().any(|&(j, _)| j == i))
                    .count()
            })
            .collect();
        assert_eq!(per_epoch, vec![4, 4, 3, 4, 4, 2, 4]);
        assert_eq!(obs(&f, "G12", "L1C")[3], (3, 109_782_559_341, b'1', b'7'));
        assert_eq!(obs(&f, "R03", "D1C")[0], (0, 2_345_678, b' ', b'7'));

        // Written zeros are observations, with their flags.
        let last = |sat| obs(&f, sat, "D1C").last().copied();
        assert_eq!(last("G05"), Some((6, 0, b' ', b' ')));
        assert_eq!(last("G12"), Some((6, 0, b' ', b'7')));
        assert_eq!(last("E11"), Some((6, 0, b'0', b' ')));
    }

    #[test]
    fn crlf_and_short_lines() {
        let crlf = TABLE_A4.replace('\n', "\r\n");
        let a = parse(TABLE_A4.as_bytes(), None).unwrap();
        let b = parse(crlf.as_bytes(), None).unwrap();
        assert_eq!(a.header, b.header);
        assert_eq!(a.epochs, b.epochs);
        assert_eq!(a.observations, b.observations);
    }

    #[test]
    fn rejects_v2() {
        let v2 = "     2.11           OBSERVATION DATA    M (MIXED)           RINEX VERSION / TYPE
     2    L1    C1                                          # / TYPES OF OBSERV
                                                            END OF HEADER
";
        assert!(matches!(
            parse(v2.as_bytes(), None),
            Err(Rinex3Error::UnsupportedVersion {
                major: 2,
                minor: 11
            })
        ));
    }

    #[test]
    fn rejects_missing_obs_types() {
        let v3 = "     3.05           OBSERVATION DATA    M                   RINEX VERSION / TYPE
                                                            END OF HEADER
";
        assert!(matches!(
            parse(v3.as_bytes(), None),
            Err(Rinex3Error::NoObsCodes)
        ));
    }

    #[test]
    fn non_ascii_comment_does_not_panic() {
        let text = TABLE_A4.replacen("A 9080      ", "A 9080 Zürich", 1);
        let f = parse(text.as_bytes(), None).unwrap();
        assert_eq!(f.epochs.len(), 5);

        // A multi-byte character straddling the column-60 label boundary.
        let comment = format!("{}\u{fc}COMMENT\n", "x".repeat(59));
        let text = TABLE_A4.replacen("A 9080", &format!("{comment}A 9080"), 1);
        let f = parse(text.as_bytes(), None).unwrap();
        assert_eq!(f.epochs.len(), 5);
    }

    #[test]
    fn non_utf8_header_text_is_read_lossily() {
        // Latin-1 in a comment must not reject the file or move the labels.
        let line = |content: &[u8], label: &str| {
            let mut line = content.to_vec();
            line.resize(60, b' ');
            line.extend_from_slice(label.as_bytes());
            line.push(b'\n');
            line
        };
        let mut text = line(
            b"     3.04           OBSERVATION DATA    M",
            "RINEX VERSION / TYPE",
        );
        text.extend(line(b"caf\xe9", "COMMENT"));
        text.extend(line(b"G    1 C1C", "SYS / # / OBS TYPES"));
        text.extend(line(b"", "END OF HEADER"));
        text.extend(b"> 2024 10 01 00 00  0.0000000  0  1\nG01  20000000.000\n");
        let f = parse(&text, None).unwrap();
        assert_eq!(f.epochs.len(), 1);
        assert_eq!(f.header.get_line("COMMENT"), Some("caf?".to_string()));
    }

    #[test]
    fn stray_and_rejected_records_are_counted() {
        let h = |content: &str, label: &str| format!("{content:<60}{label}\n");
        let mut text = h(
            "     3.04           OBSERVATION DATA    M",
            "RINEX VERSION / TYPE",
        );
        text += &h("G    1 C1C", "SYS / # / OBS TYPES");
        text += &h("", "END OF HEADER");
        text += "not a record\n";
        text += "> 2024 10 01 00 00  0.0000000  0  2\n";
        text += "G01  20000000.000\n";
        text += "g02  20000001.000\n"; // lower-case system: rejected
        let f = parse(text.as_bytes(), None).unwrap();
        assert_eq!(f.epochs.len(), 1);
        assert_eq!(f.observations.len(), 1);
        assert_eq!(f.skipped_records, 2);
    }

    #[test]
    fn truncated_epochs_and_events_are_counted() {
        let h = |content: &str, label: &str| format!("{content:<60}{label}\n");
        let mut head = h(
            "     3.04           OBSERVATION DATA    M",
            "RINEX VERSION / TYPE",
        );
        head += &h("G    1 C1C", "SYS / # / OBS TYPES");
        head += &h("", "END OF HEADER");
        let short = "> 2024 10 01 00 00  0.0000000  0  2\nG01  20000000.000\n";
        let next = "> 2024 10 01 00 00  1.0000000  0  1\nG01  20000001.000\n";
        // An epoch missing a satellite record, at the end of the file and
        // before the next epoch.
        for (text, epochs) in [
            (format!("{head}{short}"), 1),
            (format!("{head}{short}{next}"), 2),
        ] {
            let f = parse(text.as_bytes(), None).unwrap();
            assert_eq!(f.epochs.len(), epochs);
            assert_eq!(f.skipped_records, 1, "{text}");
        }
        // Events missing records: a cycle slip before the next epoch and a
        // comment event at the end of the file.
        let slip = "> 2024 10 01 00 00  0.0000000  6  2\nG01  20000000.000\n";
        let comment = "> 2024 10 01 00 00  1.0000000  4  2\nA COMMENT\n";
        for text in [
            format!("{head}{slip}{next}"),
            format!("{head}{next}{comment}"),
        ] {
            let f = parse(text.as_bytes(), None).unwrap();
            assert_eq!(f.events.len(), 1);
            assert_eq!(f.skipped_records, 1, "{text}");
        }
    }

    #[test]
    fn a_duplicate_satellite_keeps_its_first_record() {
        let h = |content: &str, label: &str| format!("{content:<60}{label}\n");
        let mut text = h(
            "     3.04           OBSERVATION DATA    M",
            "RINEX VERSION / TYPE",
        );
        text += &h("G    1 C1C", "SYS / # / OBS TYPES");
        text += &h("", "END OF HEADER");
        text += "> 2024 10 01 00 00  0.0000000  0  2\n";
        text += "G01  20000000.000\nG01  20000001.000\n";
        let f = parse(text.as_bytes(), None).unwrap();
        assert_eq!(f.epochs.len(), 1);
        assert_eq!(f.skipped_records, 1);
        assert_eq!(obs(&f, "G01", "C1C"), vec![(0, 20_000_000_000, b' ', b' ')]);
    }

    #[test]
    fn a_duplicate_code_is_an_error() {
        // Both slots would feed one (satellite, code) series, giving it two
        // values for the epoch.
        let h = |content: &str, label: &str| format!("{content:<60}{label}\n");
        let mut text = h(
            "     3.04           OBSERVATION DATA    M",
            "RINEX VERSION / TYPE",
        );
        text += &h("G    2 C1C C1C", "SYS / # / OBS TYPES");
        text += &h("", "END OF HEADER");
        text += "> 2024 10 01 00 00  0.0000000  0  1\n";
        text += "G01  20000000.000  20000001.000\n";
        assert!(matches!(
            parse(text.as_bytes(), None),
            Err(Rinex3Error::Header(HeaderError::DuplicateObsCode { .. }))
        ));
    }

    // ---- fixture tests -----------------------------------------------

    #[test]
    fn fixture_balm_v304() {
        let f = parse(&crate::test_rinex("balm_v304_1hz"), None).unwrap();
        assert_eq!(f.header.get_line("MARKER NAME"), Some("BALM".to_string()));
        // Two minutes at 1 s.
        assert_eq!(f.epochs.len(), 120);
        let num_epochs = f.epochs.len();
        for series in f.observations.values() {
            assert!(series.windows(2).all(|w| w[0].0 < w[1].0));
            assert!(series.last().map_or(0, |&(i, _)| i) < num_epochs);
        }
        let systems: std::collections::BTreeSet<u8> =
            f.observations.keys().map(|(s, _)| s.system.0).collect();
        // BALM lists Galileo satellites but reports no observations for them.
        assert_eq!(systems, b"CGJR".iter().copied().collect());
        let e0 = &f.epochs[0];
        assert_eq!((e0.year, e0.month, e0.day), (2026, 10, 8));
        // First epoch, first satellite: G05 C1C = 25053222.110, no flags.
        let g05 = SatId::new(SatSystem::GPS, 5);
        let c1c = ObsCode::from_v3("C1C").unwrap();
        let &(i, o) = f.observations[&(g05, c1c)].first().unwrap();
        assert_eq!((i, o.value, o.lli, o.ssi), (0, 25_053_222_110, b' ', b' '));
    }

    #[test]
    fn events_change_the_observation_types() {
        let h = |content: &str, label: &str| format!("{content:<60}{label}\n");
        let mut text = h(
            "     3.04           OBSERVATION DATA    M",
            "RINEX VERSION / TYPE",
        );
        text += &h("G    3 C1C L1C S1C", "SYS / # / OBS TYPES");
        text += &h("E    1 C1C", "SYS / # / OBS TYPES");
        text += &h("", "END OF HEADER");
        text += "> 2024 10 01 00 00  0.0000000  0  2\n";
        text += "G01  20001000.000 6 100005000.000 6        45.000\n";
        text += "E02  30001000.000\n";
        text += "> 2024 10 01 00 00  0.5000000  4  1\n";
        text += &h("G    2 S1C C1C", "SYS / # / OBS TYPES");
        text += "> 2024 10 01 00 00  1.0000000  0  2\n";
        text += "G01        46.000    20001001.000 6\n";
        text += "E02  30001001.000\n";
        let f = parse(text.as_bytes(), None).unwrap();
        let value = |sat: SatId, code: &str| {
            f.observations[&(sat, ObsCode::from_v3(code).unwrap())]
                .iter()
                .find(|&&(i, _)| i == 1)
                .map(|&(_, o)| o.value)
        };
        let g01 = SatId::new(SatSystem::GPS, 1);
        assert_eq!(value(g01, "S1C"), Some(46_000));
        assert_eq!(value(g01, "C1C"), Some(20_001_001_000));
        assert_eq!(value(g01, "L1C"), None);
        // Other systems keep their codes.
        assert_eq!(
            value(SatId::new(SatSystem::GALILEO, 2), "C1C"),
            Some(30_001_001_000)
        );
    }
}
