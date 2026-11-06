//! Parser for the `EPOC` and `EPO2` (Epoch Records) chunks' payloads.
//!
//! Per `docs/Succinct_RINEX.md` §EPOC and §EPO2:
//!
//! ```text
//! ULEB128 total_epoch_count
//! one or more spans:
//!   SLEB128 interval
//!   ULEB128 (epochs_in_span - 1)
//!   ULEB128 packed_date  =  year*10000 + month*100 + day
//!   ULEB128 packed_time  =  (hour*10000 + minute*100 + second) * U
//! optional RLE-encoded receiver clock offsets:
//!   pairs of  (SLEB128 value*1e12, ULEB128 (run_count - 1))
//! ```
//!
//! The two chunks differ only in their time unit `U`: 10⁷ per second for
//! `EPOC` and 10¹² per second (picoseconds) for `EPO2`.  At most one of the
//! two is present in a file; the spec calls whichever it is the *epoch
//! chunk*.
//!
//! ## Interval encoding
//!
//! The SLEB128 interval is asymmetric for compactness:
//! * **negative** values mean integer-second intervals — interval is
//!   `-stored_value` seconds (so `SLEB128(-1)` is 1 s, `SLEB128(-30)` is 30 s).
//! * **positive** values are intervals expressed as `seconds × U` (so
//!   `SLEB128(5_000_000)` is 0.5 s in `EPOC`).
//!
//! ## Time advancement
//!
//! Within a span, each successive epoch is at `prev + interval`. When the
//! resulting seconds value is **at least 60**, the excess carries into the
//! next higher unit (minutes, then hours).  The sole exception is the
//! leap-second range `23:59:60` up to (but excluding) `23:59:61`, which is
//! preserved as-is.  A timestamp `>= 23:59:61` is illegal per the spec.
//! Encoders MUST start a new span at each day boundary (the first epoch of
//! each day), so within a single span hours stay in `0..=23`.
//!
//! ## Clock offsets
//!
//! Stored as `(value, run_count - 1)` pairs. Trailing epochs not covered by
//! any pair are implicitly zero. RINEX 2.x files SHOULD omit clock offsets
//! entirely, in which case the section is empty and all offsets are zero.

use thiserror::Error;

use crate::leb128::{read_sleb128, read_uleb128, LebError};

#[derive(Debug, Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum EpocError {
    #[error("invalid LEB128 in epoch chunk: {0}")]
    Leb(#[from] LebError),
    #[error("epoch chunk declared {declared} epochs but spans cover {covered}")]
    SpanCountMismatch { declared: u64, covered: u64 },
    #[error("epoch chunk date encoding {0} is invalid")]
    BadDate(u64),
    #[error("epoch chunk time encoding {0} is invalid")]
    BadTime(u64),
    #[error("epoch chunk interval {0} is out of range")]
    BadInterval(i64),
    #[error("epoch chunk epoch count {0} exceeds addressable memory")]
    TooManyEpochs(u64),
}

/// Which epoch chunk a payload came from, which fixes its time unit.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum EpochChunkKind {
    /// `EPOC`: times in units of 10⁻⁷ s.
    #[default]
    Epoc,
    /// `EPO2` (revision 1.2): times in units of 10⁻¹² s.
    Epo2,
}

impl EpochChunkKind {
    /// The chunk's time unit `U`, in units per second.
    #[must_use]
    pub const fn units_per_second(self) -> u64 {
        match self {
            Self::Epoc => 10_000_000,
            Self::Epo2 => 1_000_000_000_000,
        }
    }

    /// The chunk's FOURCC tag.
    #[must_use]
    pub const fn fourcc(self) -> [u8; 4] {
        match self {
            Self::Epoc => *b"EPOC",
            Self::Epo2 => *b"EPO2",
        }
    }
}

/// One materialized epoch record. Time fields are kept verbatim from the
/// span expansion: `seconds_e7` may exceed `60 × 10⁷` during leap seconds,
/// and `hour` may grow past 23 if a (non-conformant) span crosses a day.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Epoch {
    pub year: u16,
    pub month: u8,
    pub day: u8,
    pub hour: u8,
    pub minute: u8,
    /// Seconds within the minute, as `s × 10⁷`, floored.
    pub seconds_e7: u64,
    /// Picoseconds beyond `seconds_e7` (0–99 999).  Always zero for `EPOC`.
    pub extra_ps: u32,
    /// Receiver clock offset, as `value × 10¹²`. Zero if not specified.
    pub clock_offset_e12: i64,
}

impl Epoch {
    /// Seconds within the minute, in units of `1 / units_per_second` s.
    /// `units_per_second` must be a multiple of 10⁷ that divides 10¹².
    #[must_use]
    pub fn seconds_in(&self, units_per_second: u64) -> u64 {
        let ps_per_unit = 1_000_000_000_000 / units_per_second;
        self.seconds_e7 * (units_per_second / 10_000_000) + u64::from(self.extra_ps) / ps_per_unit
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Epoc {
    pub kind: EpochChunkKind,
    pub epochs: Vec<Epoch>,
}

pub fn parse(payload: &[u8], kind: EpochChunkKind) -> Result<Epoc, EpocError> {
    let units = kind.units_per_second();
    let mut cursor: &[u8] = payload;
    let total = read_uleb128(&mut cursor)?;
    let total_usize = usize::try_from(total).map_err(|_| EpocError::TooManyEpochs(total))?;

    if total == 0 {
        return Ok(Epoc {
            kind,
            epochs: Vec::new(),
        });
    }

    let mut epochs: Vec<Epoch> = Vec::with_capacity(total_usize);
    while u64::try_from(epochs.len()).expect("len fits in u64") < total {
        let interval = read_sleb128(&mut cursor)?;
        let span_count_m1 = read_uleb128(&mut cursor)?;
        let span_count = span_count_m1
            .checked_add(1)
            .ok_or(EpocError::Leb(LebError::Overflow))?;
        let date = read_uleb128(&mut cursor)?;
        let time = read_uleb128(&mut cursor)?;

        let (year, month, day) = decode_date(date)?;
        let (mut hour, mut minute, mut seconds) = decode_time(time, units)?;
        let interval_units = interval_to_units(interval, units)?;

        for i in 0..span_count {
            if i > 0 {
                advance_time(&mut hour, &mut minute, &mut seconds, interval_units, units);
            }
            let (seconds_e7, extra_ps) = split_seconds(seconds, units);
            epochs.push(Epoch {
                year,
                month,
                day,
                hour,
                minute,
                seconds_e7,
                extra_ps,
                clock_offset_e12: 0,
            });
            if u64::try_from(epochs.len()).expect("len fits in u64") == total {
                break;
            }
        }
    }

    let covered = u64::try_from(epochs.len()).expect("len fits in u64");
    if covered != total {
        return Err(EpocError::SpanCountMismatch {
            declared: total,
            covered,
        });
    }

    // RLE clock offsets, if any. Trailing zeros are implicit.
    let mut idx: usize = 0;
    while !cursor.is_empty() && idx < epochs.len() {
        let value = read_sleb128(&mut cursor)?;
        let run_count_m1 = read_uleb128(&mut cursor)?;
        let run_count = usize::try_from(run_count_m1.saturating_add(1)).unwrap_or(usize::MAX);
        let end = idx.saturating_add(run_count).min(epochs.len());
        for slot in &mut epochs[idx..end] {
            slot.clock_offset_e12 = value;
        }
        idx = end;
    }

    Ok(Epoc { kind, epochs })
}

/// Split seconds-within-minute in chunk units into `(seconds_e7, extra_ps)`,
/// flooring to 10⁻⁷ s as RINEX 4.02 does when the extra digits are present.
fn split_seconds(seconds: u64, units: u64) -> (u64, u32) {
    let units_per_e7 = units / 10_000_000;
    let ps_per_unit = 1_000_000_000_000 / units;
    let extra = (seconds % units_per_e7) * ps_per_unit;
    (
        seconds / units_per_e7,
        u32::try_from(extra).expect("sub-1e-7 s remainder fits in u32"),
    )
}

fn decode_date(packed: u64) -> Result<(u16, u8, u8), EpocError> {
    let day_u64 = packed % 100;
    let month_u64 = (packed / 100) % 100;
    let year_part = packed / 10_000;

    let year: u16 = if year_part < 100 {
        // RINEX 2 two-digit year: 80..=99 → 1980..=1999, 0..=79 → 2000..=2079.
        let y = if year_part >= 80 {
            year_part + 1900
        } else {
            year_part + 2000
        };
        u16::try_from(y).map_err(|_| EpocError::BadDate(packed))?
    } else {
        u16::try_from(year_part).map_err(|_| EpocError::BadDate(packed))?
    };

    if !(1..=12).contains(&month_u64) || !(1..=31).contains(&day_u64) {
        return Err(EpocError::BadDate(packed));
    }

    // Casts safe by the range checks above.
    #[allow(clippy::cast_possible_truncation)]
    Ok((year, month_u64 as u8, day_u64 as u8))
}

fn decode_time(packed: u64, units: u64) -> Result<(u8, u8, u64), EpocError> {
    let per_minute = 100 * units;
    let per_hour = 100 * per_minute;

    let hour_u64 = packed / per_hour;
    let after_hour = packed % per_hour;
    let minute_u64 = after_hour / per_minute;
    let seconds = after_hour % per_minute;

    if hour_u64 > 23 || minute_u64 > 59 {
        return Err(EpocError::BadTime(packed));
    }

    #[allow(clippy::cast_possible_truncation)]
    Ok((hour_u64 as u8, minute_u64 as u8, seconds))
}

fn interval_to_units(interval: i64, units: u64) -> Result<i64, EpocError> {
    if interval < 0 {
        interval
            .checked_neg()
            .and_then(|n| n.checked_mul(i64::try_from(units).ok()?))
            .ok_or(EpocError::BadInterval(interval))
    } else {
        Ok(interval)
    }
}

#[allow(clippy::cast_possible_truncation)]
pub(crate) fn advance_time(
    hour: &mut u8,
    minute: &mut u8,
    seconds: &mut u64,
    interval: i64,
    units: u64,
) {
    let new_s = i128::from(*seconds) + i128::from(interval);
    if new_s < 0 {
        *seconds = 0;
    } else {
        let new_s64 = u64::try_from(new_s).unwrap_or(u64::MAX);
        let per_minute = 60 * units;
        // Leap-second exception: 23:59:60 up to 23:59:61 is preserved
        // as-is.  The spec forbids >= 23:59:61 and requires new spans at
        // day boundaries, so we don't need to carry past this.
        if *hour == 23 && *minute == 59 && new_s64 < 61 * units {
            *seconds = new_s64;
        } else if new_s64 >= per_minute {
            let carry_minutes = new_s64 / per_minute;
            *seconds = new_s64 % per_minute;
            let new_minute = u64::from(*minute) + carry_minutes;
            if new_minute >= 60 {
                *hour = hour.saturating_add(u8::try_from(new_minute / 60).unwrap_or(u8::MAX));
                *minute = (new_minute % 60) as u8;
            } else {
                *minute = new_minute as u8;
            }
        } else {
            *seconds = new_s64;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn enc_uleb(mut v: u64) -> Vec<u8> {
        let mut out = vec![];
        loop {
            let b = (v & 0x7F) as u8;
            v >>= 7;
            if v == 0 {
                out.push(b);
                return out;
            }
            out.push(b | 0x80);
        }
    }

    fn enc_sleb(v: i64) -> Vec<u8> {
        // ZigZag encode: (v << 1) ^ (v >> 63), as u64.
        #[allow(clippy::cast_sign_loss)]
        let zz = ((v << 1) ^ (v >> 63)) as u64;
        enc_uleb(zz)
    }

    /// Build an EPOC payload from one or more spans, optionally followed by
    /// (value, run-count-minus-1) clock offset pairs.
    fn build_payload(
        total: u64,
        spans: &[(i64, u64, u64, u64)], // (interval, count_m1, date, time)
        clock: &[(i64, u64)],
    ) -> Vec<u8> {
        let mut out = enc_uleb(total);
        for (interval, count_m1, date, time) in spans {
            out.extend(enc_sleb(*interval));
            out.extend(enc_uleb(*count_m1));
            out.extend(enc_uleb(*date));
            out.extend(enc_uleb(*time));
        }
        for (value, count_m1) in clock {
            out.extend(enc_sleb(*value));
            out.extend(enc_uleb(*count_m1));
        }
        out
    }

    #[test]
    fn empty_epoc() {
        let payload = enc_uleb(0);
        let e = parse(&payload, EpochChunkKind::Epoc).unwrap();
        assert_eq!(e.epochs, []);
    }

    #[test]
    fn single_span_integer_interval() {
        // 10 epochs, 1 s interval (encoded as -1), starting 2020-07-18 00:00:00.
        let p = build_payload(10, &[(-1, 9, 20_200_718, 0)], &[]);
        let e = parse(&p, EpochChunkKind::Epoc).unwrap();
        assert_eq!(e.epochs.len(), 10);
        for (i, ep) in e.epochs.iter().enumerate() {
            assert_eq!(ep.year, 2020);
            assert_eq!(ep.month, 7);
            assert_eq!(ep.day, 18);
            assert_eq!(ep.hour, 0);
            assert_eq!(ep.minute, 0);
            assert_eq!(ep.seconds_e7, (i as u64) * 10_000_000);
            assert_eq!(ep.clock_offset_e12, 0);
        }
    }

    #[test]
    fn single_span_fractional_interval() {
        // 4 epochs at 0.5 s starting at 00:00:00 → seconds 0.0, 0.5, 1.0, 1.5.
        let p = build_payload(4, &[(5_000_000, 3, 20_200_718, 0)], &[]);
        let e = parse(&p, EpochChunkKind::Epoc).unwrap();
        let secs: Vec<u64> = e.epochs.iter().map(|ep| ep.seconds_e7).collect();
        assert_eq!(secs, vec![0, 5_000_000, 10_000_000, 15_000_000]);
    }

    #[test]
    fn minute_rollover_on_exact_60() {
        // 70 epochs at 1 s starting at 00:00:00 → epoch 60 is 00:01:00.
        let p = build_payload(70, &[(-1, 69, 20_200_718, 0)], &[]);
        let e = parse(&p, EpochChunkKind::Epoc).unwrap();
        assert_eq!(e.epochs[59].minute, 0);
        assert_eq!(e.epochs[59].seconds_e7, 59 * 10_000_000);
        assert_eq!(e.epochs[60].minute, 1);
        assert_eq!(e.epochs[60].seconds_e7, 0);
        assert_eq!(e.epochs[69].minute, 1);
        assert_eq!(e.epochs[69].seconds_e7, 9 * 10_000_000);
    }

    #[test]
    fn steps_carry_several_minutes_and_hours() {
        // 5-minute and 5-hour intervals in one span each: every step
        // carries more than one minute, as the spec requires readers to do.
        let p = build_payload(3, &[(-300, 2, 20_241_001, 0)], &[]);
        let e = parse(&p, EpochChunkKind::Epoc).unwrap();
        let hm: Vec<(u8, u8)> = e.epochs.iter().map(|e| (e.hour, e.minute)).collect();
        assert_eq!(hm, [(0, 0), (0, 5), (0, 10)]);
        let start = 55 * 1_000_000_000 + 30 * 10_000_000; // 00:55:30
        let p = build_payload(3, &[(-5 * 3600 - 90, 2, 20_241_001, start)], &[]);
        let e = parse(&p, EpochChunkKind::Epoc).unwrap();
        let hms: Vec<(u8, u8, u64)> = e
            .epochs
            .iter()
            .map(|e| (e.hour, e.minute, e.seconds_e7))
            .collect();
        assert_eq!(
            hms,
            [(0, 55, 300_000_000), (5, 57, 0), (10, 58, 300_000_000)]
        );
    }

    #[test]
    fn hour_rollover_on_exact_60_minute() {
        // Start at 00:59:50, advance 1 s for 20 epochs.
        // At i=10 we hit 01:00:00 (m += 1 → 60 → reset to 0, h += 1).
        let start_time = 59 * 1_000_000_000 + 50 * 10_000_000;
        let p = build_payload(20, &[(-1, 19, 20_200_718, start_time)], &[]);
        let e = parse(&p, EpochChunkKind::Epoc).unwrap();
        assert_eq!(
            (e.epochs[9].hour, e.epochs[9].minute, e.epochs[9].seconds_e7),
            (0, 59, 59 * 10_000_000)
        );
        assert_eq!(
            (
                e.epochs[10].hour,
                e.epochs[10].minute,
                e.epochs[10].seconds_e7
            ),
            (1, 0, 0)
        );
    }

    #[test]
    fn leap_second_preserved_through_60_dot_9() {
        // Start at 23:59:60.0 (leap second), 0.5 s interval, 4 epochs.
        // The leap-second range 23:59:60.0 – 23:59:60.9999999 is preserved
        // without carrying.  23:59:61.0 is illegal per the spec, and
        // conformant encoders MUST split spans at leap seconds, but we
        // still handle it gracefully by carrying.
        let start = 23 * 100_000_000_000 + 59 * 1_000_000_000 + 60 * 10_000_000;
        let p = build_payload(4, &[(5_000_000, 3, 20_200_718, start)], &[]);
        let e = parse(&p, EpochChunkKind::Epoc).unwrap();
        // Epoch 0: 23:59:60.0 (leap second, preserved)
        assert_eq!(e.epochs[0].hour, 23);
        assert_eq!(e.epochs[0].minute, 59);
        assert_eq!(e.epochs[0].seconds_e7, 600_000_000);
        // Epoch 1: 23:59:60.5 (still in leap-second range, preserved)
        assert_eq!(e.epochs[1].hour, 23);
        assert_eq!(e.epochs[1].minute, 59);
        assert_eq!(e.epochs[1].seconds_e7, 605_000_000);
        // Epoch 2: 61.0 s >= 61 → carry → 24:00:01.0
        assert_eq!(e.epochs[2].hour, 24);
        assert_eq!(e.epochs[2].minute, 0);
        assert_eq!(e.epochs[2].seconds_e7, 10_000_000);
        // Epoch 3: 24:00:01.5
        assert_eq!(e.epochs[3].hour, 24);
        assert_eq!(e.epochs[3].minute, 0);
        assert_eq!(e.epochs[3].seconds_e7, 15_000_000);
    }

    #[test]
    fn multi_span() {
        // 5 epochs at 1s, then day boundary, then 5 more epochs at 1s.
        let span1 = (
            -1,
            4,
            20_200_718,
            23 * 100_000_000_000 + 59 * 1_000_000_000 + 55 * 10_000_000,
        );
        let span2 = (-1, 4, 20_200_719, 0);
        let p = build_payload(10, &[span1, span2], &[]);
        let e = parse(&p, EpochChunkKind::Epoc).unwrap();
        assert_eq!(e.epochs.len(), 10);
        assert_eq!(
            (
                e.epochs[0].day,
                e.epochs[0].hour,
                e.epochs[0].minute,
                e.epochs[0].seconds_e7
            ),
            (18, 23, 59, 55 * 10_000_000)
        );
        assert_eq!(
            (
                e.epochs[4].day,
                e.epochs[4].hour,
                e.epochs[4].minute,
                e.epochs[4].seconds_e7
            ),
            (18, 23, 59, 59 * 10_000_000)
        );
        assert_eq!(
            (
                e.epochs[5].day,
                e.epochs[5].hour,
                e.epochs[5].minute,
                e.epochs[5].seconds_e7
            ),
            (19, 0, 0, 0)
        );
    }

    #[test]
    fn clock_offsets_basic() {
        // 5 epochs, two clock-offset runs: value=100 ×3, value=-50 ×2.
        let p = build_payload(5, &[(-1, 4, 20_200_718, 0)], &[(100, 2), (-50, 1)]);
        let e = parse(&p, EpochChunkKind::Epoc).unwrap();
        let offsets: Vec<i64> = e.epochs.iter().map(|ep| ep.clock_offset_e12).collect();
        assert_eq!(offsets, vec![100, 100, 100, -50, -50]);
    }

    #[test]
    fn clock_offsets_padded_with_zeros() {
        // 5 epochs but only one offset run covering 2 epochs → trailing 3 are zero.
        let p = build_payload(5, &[(-1, 4, 20_200_718, 0)], &[(42, 1)]);
        let e = parse(&p, EpochChunkKind::Epoc).unwrap();
        let offsets: Vec<i64> = e.epochs.iter().map(|ep| ep.clock_offset_e12).collect();
        assert_eq!(offsets, vec![42, 42, 0, 0, 0]);
    }

    #[test]
    fn clock_offsets_omitted_means_all_zero() {
        // RINEX 2 case: no clock offsets at all.
        let p = build_payload(3, &[(-1, 2, 20_200_718, 0)], &[]);
        let e = parse(&p, EpochChunkKind::Epoc).unwrap();
        for ep in &e.epochs {
            assert_eq!(ep.clock_offset_e12, 0);
        }
    }

    #[test]
    fn rinex2_year_disambiguation() {
        // year_part=99 → 1999, year_part=20 → 2020.
        let p = build_payload(1, &[(-1, 0, 99 * 10_000 + 100 + 2, 0)], &[]);
        let e = parse(&p, EpochChunkKind::Epoc).unwrap();
        assert_eq!(e.epochs[0].year, 1999);
        assert_eq!(e.epochs[0].month, 1);
        assert_eq!(e.epochs[0].day, 2);

        let p = build_payload(1, &[(-1, 0, 20 * 10_000 + 7 * 100 + 18, 0)], &[]);
        let e = parse(&p, EpochChunkKind::Epoc).unwrap();
        assert_eq!(e.epochs[0].year, 2020);
    }

    #[test]
    fn span_count_mismatch_detected() {
        // Declare 10 epochs but span only covers 3.
        let p = build_payload(10, &[(-1, 2, 20_200_718, 0)], &[]);
        // First span supplies 3 epochs; loop will try to read another span
        // header from the (likely empty) cursor → either Leb error or
        // SpanCountMismatch on subsequent iterations. The exact failure
        // depends on the truncation point. Assert that we get *some* error.
        assert!(parse(&p, EpochChunkKind::Epoc).is_err());
    }

    const PS_PER_MINUTE: u64 = 100_000_000_000_000;
    const PS_PER_HOUR: u64 = 100 * PS_PER_MINUTE;

    #[test]
    fn epo2_time_of_day_and_fractional_interval() {
        // 13:10:36.123456712345, then +1.000000023455 s → 13:10:37.123456735800.
        let start = 13 * PS_PER_HOUR + 10 * PS_PER_MINUTE + 36_123_456_712_345;
        let p = build_payload(2, &[(1_000_000_023_455, 1, 20_241_001, start)], &[]);
        let e = parse(&p, EpochChunkKind::Epo2).unwrap();
        assert_eq!(e.kind, EpochChunkKind::Epo2);
        let t: Vec<_> = e
            .epochs
            .iter()
            .map(|ep| (ep.hour, ep.minute, ep.seconds_e7, ep.extra_ps))
            .collect();
        assert_eq!(
            t,
            [(13, 10, 361_234_567, 12_345), (13, 10, 371_234_567, 35_800)]
        );
        assert_eq!(
            e.epochs[1].seconds_in(1_000_000_000_000),
            37_123_456_735_800
        );
        assert_eq!(e.epochs[1].seconds_in(10_000_000), 371_234_567);
    }

    #[test]
    fn epo2_negative_interval_is_whole_seconds() {
        // 3 epochs at 30 s from 00:00:59.000000000001 carry into minutes.
        let p = build_payload(3, &[(-30, 2, 20_241_001, 59_000_000_000_001)], &[]);
        let e = parse(&p, EpochChunkKind::Epo2).unwrap();
        let t: Vec<_> = e
            .epochs
            .iter()
            .map(|ep| (ep.minute, ep.seconds_e7, ep.extra_ps))
            .collect();
        assert_eq!(
            t,
            [
                (0, 590_000_000, 1),
                (1, 290_000_000, 1),
                (1, 590_000_000, 1)
            ]
        );
    }

    #[test]
    fn epo2_leap_second_preserved_below_61() {
        let start = 23 * PS_PER_HOUR + 59 * PS_PER_MINUTE + 60_000_000_000_000;
        let p = build_payload(2, &[(999_999_999_999, 1, 20_161_231, start)], &[]);
        let e = parse(&p, EpochChunkKind::Epo2).unwrap();
        let last = e.epochs[1];
        assert_eq!(
            (last.hour, last.minute, last.seconds_e7, last.extra_ps),
            (23, 59, 609_999_999, 99_999)
        );
    }

    #[test]
    fn epoc_never_has_extra_ps() {
        let p = build_payload(2, &[(1, 1, 20_241_001, 0)], &[]);
        let e = parse(&p, EpochChunkKind::Epoc).unwrap();
        assert!(e.epochs.iter().all(|ep| ep.extra_ps == 0));
        assert_eq!(e.epochs[1].seconds_e7, 1);
    }

    #[test]
    fn overflowing_whole_second_interval_is_rejected() {
        // 10^7 s is fine as 10^14 units of 1e-7 s but overflows i64 picoseconds.
        let p = build_payload(2, &[(-10_000_000, 1, 20_241_001, 0)], &[]);
        assert!(parse(&p, EpochChunkKind::Epoc).is_ok());
        assert_eq!(
            parse(&p, EpochChunkKind::Epo2),
            Err(EpocError::BadInterval(-10_000_000))
        );
    }

    fn fixture_epoch_chunk(name: &str) -> Vec<u8> {
        let bytes = crate::test_data(name);
        let file = crate::file::SrnxFile::parse(&bytes).unwrap();
        let sdir = file.sdir().unwrap().unwrap();
        file.chunk_at_offset(sdir.epoc_offset)
            .unwrap()
            .payload
            .to_vec()
    }

    #[test]
    fn fixture_at01_epoc() {
        let e = parse(
            &fixture_epoch_chunk("at01_v211_1hz.srnx"),
            EpochChunkKind::Epoc,
        )
        .unwrap();
        assert_eq!(e.epochs.len(), 120);
        let first = &e.epochs[0];
        assert_eq!((first.year, first.month, first.day), (2020, 7, 18));
        assert_eq!((first.hour, first.minute, first.seconds_e7), (0, 0, 0));
        let last = e.epochs.last().unwrap();
        assert_eq!(
            (last.hour, last.minute, last.seconds_e7),
            (0, 1, 59 * 10_000_000)
        );
        // RINEX 2: clock offsets all zero.
        assert!(e.epochs.iter().all(|ep| ep.clock_offset_e12 == 0));
    }

    #[test]
    fn fixture_mty2_epoc_skips_event() {
        // 80 records at 15 s, one of which is a flag-4 event.
        let e = parse(
            &fixture_epoch_chunk("mty2_v211_15s.srnx"),
            EpochChunkKind::Epoc,
        )
        .unwrap();
        assert_eq!(e.epochs.len(), 79);
        let first = &e.epochs[0];
        assert_eq!((first.hour, first.minute, first.seconds_e7), (0, 50, 0));
        let last = e.epochs.last().unwrap();
        assert_eq!(
            (last.hour, last.minute, last.seconds_e7),
            (1, 9, 30 * 10_000_000)
        );
    }

    #[test]
    fn fixture_balm_epo2() {
        let e = parse(
            &fixture_epoch_chunk("balm_v402_1hz_ps.srnx"),
            EpochChunkKind::Epo2,
        )
        .unwrap();
        assert_eq!(e.epochs.len(), 120);
        let first = &e.epochs[0];
        assert_eq!((first.year, first.month, first.day), (2026, 10, 8));
        assert_eq!((first.seconds_e7, first.extra_ps), (0, 37_123));
        assert_eq!(
            (e.epochs[1].seconds_e7, e.epochs[1].extra_ps),
            (10_000_000, 74_246)
        );
    }
}
