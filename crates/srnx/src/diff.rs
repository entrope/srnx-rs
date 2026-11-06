//! Semantic comparison of two [`ObsFile`] instances.
//!
//! Entry point: [`compare`], which returns a `Vec<String>` of human-readable
//! difference descriptions (empty when the files are equivalent).
//!
//! This module is used by both the `srnx-diff` CLI and the corpus integration
//! tests so the comparison logic lives in exactly one place.

use std::collections::BTreeSet;
use std::fmt::Write;

use crate::model::{EpochRecord, Obs, ObsFile};

// ---------------------------------------------------------------------------
// Public helpers (also useful standalone)
// ---------------------------------------------------------------------------

/// Format an epoch record as an ISO-like timestamp string.  Picosecond
/// digits (RINEX 4.02) are appended only when nonzero.
#[must_use]
pub fn fmt_epoch(e: &EpochRecord) -> String {
    let whole = e.seconds_e7 / 10_000_000;
    let frac = e.seconds_e7 % 10_000_000;
    let mut out = format!(
        "{:04}-{:02}-{:02} {:02}:{:02}:{:02}.{:07}",
        e.year, e.month, e.day, e.hour, e.minute, whole, frac
    );
    if e.extra_ps != 0 {
        let _ = write!(out, "{:05}", e.extra_ps);
    }
    out
}

/// Compare two epoch records on timestamp fields only (ignores flag and
/// clock offset).
#[must_use]
pub fn epochs_eq(a: &EpochRecord, b: &EpochRecord) -> bool {
    a.year == b.year
        && a.month == b.month
        && a.day == b.day
        && a.hour == b.hour
        && a.minute == b.minute
        && a.seconds_e7 == b.seconds_e7
        && a.extra_ps == b.extra_ps
}

/// Format an observation as `value=<milli-units> lli='<c>' ssi='<c>'`.
fn fmt_obs(o: Obs) -> String {
    format!(
        "value={} lli={:?} ssi={:?}",
        o.value,
        char::from(o.lli),
        char::from(o.ssi)
    )
}

/// Compare two observations on value, LLI, and SSI.  An indicator of `0`
/// differs from a blank one: SRNX keeps them apart.
#[must_use]
pub fn obs_eq(a: Obs, b: Obs) -> bool {
    a == b
}

// ---------------------------------------------------------------------------
// Full-file comparison
// ---------------------------------------------------------------------------

/// Maximum number of individual differences to report before stopping.
pub const MAX_DIFFS: usize = 50;

/// Compare two observation files and return a list of human-readable
/// difference descriptions.
///
/// Returns an empty `Vec` when the files are semantically equivalent.
/// At most [`MAX_DIFFS`] differences are collected.
#[must_use]
pub fn compare(a: &ObsFile, b: &ObsFile, label_a: &str, label_b: &str) -> Vec<String> {
    let mut out = Vec::new();

    compare_headers(a, b, label_a, label_b, &mut out);
    compare_epochs(a, b, label_a, label_b, &mut out);
    compare_events(a, b, label_a, label_b, &mut out);
    compare_satellite_epochs(a, b, label_a, label_b, &mut out);
    compare_observations(a, b, label_a, label_b, &mut out);
    out
}

// ---------------------------------------------------------------------------
// Header comparison
// ---------------------------------------------------------------------------

fn compare_headers(a: &ObsFile, b: &ObsFile, label_a: &str, label_b: &str, out: &mut Vec<String>) {
    if out.len() >= MAX_DIFFS {
        return;
    }

    // HeaderText is already LF-normalized.  RINEX v2 lines are 80-char
    // fixed-width (trailing spaces); SRNX RHDR strips trailing whitespace.
    // Normalize both to trailing-space-free before comparing.
    let a_text = norm_header(a.header.as_str());
    let b_text = norm_header(b.header.as_str());

    if a_text != b_text {
        let a_lines = a_text.lines().count();
        let b_lines = b_text.lines().count();
        out.push(format!(
            "header text differs: {label_a}={a_lines} lines {label_b}={b_lines} lines"
        ));
    }
}

/// Strip trailing whitespace from each line so that 80-column RINEX v2
/// headers compare equal to SRNX RHDR headers (which drop column padding).
fn norm_header(text: &str) -> String {
    text.lines()
        .map(str::trim_end)
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

fn compare_epochs(a: &ObsFile, b: &ObsFile, label_a: &str, label_b: &str, out: &mut Vec<String>) {
    let na = a.epochs.len();
    let nb = b.epochs.len();
    if na != nb {
        out.push(format!("epoch count: {label_a}={na}  {label_b}={nb}"));
    }

    let common = na.min(nb);
    for i in 0..common {
        if out.len() >= MAX_DIFFS {
            break;
        }
        let ea = &a.epochs[i];
        let eb = &b.epochs[i];
        if !epochs_eq(ea, eb) {
            out.push(format!(
                "epoch[{i}] timestamp: {} vs {}",
                fmt_epoch(ea),
                fmt_epoch(eb)
            ));
        } else if ea.clock_offset_e12 != eb.clock_offset_e12 {
            out.push(format!(
                "epoch[{i}] clock offset: {} vs {} (×10⁻¹²)",
                ea.clock_offset_e12, eb.clock_offset_e12
            ));
        }
    }
}

/// Events must sit at the same positions, with the same flags and lines
/// (ignoring trailing whitespace, which SRNX does not keep).
fn compare_events(a: &ObsFile, b: &ObsFile, label_a: &str, label_b: &str, out: &mut Vec<String>) {
    let na = a.events.len();
    let nb = b.events.len();
    if na != nb {
        out.push(format!("event count: {label_a}={na}  {label_b}={nb}"));
    }

    for (k, ((ia, ea), (ib, eb))) in a.events.iter().zip(&b.events).enumerate() {
        if out.len() >= MAX_DIFFS {
            break;
        }
        if ia != ib {
            out.push(format!("event[{k}] position: epoch[{ia}] vs epoch[{ib}]"));
        } else if ea.flag != eb.flag {
            out.push(format!("event[{k}] flag: {} vs {}", ea.flag, eb.flag));
        } else if !epoch_lines_eq(&ea.epoch_line, &eb.epoch_line) {
            out.push(format!(
                "event[{k}] (flag {}) epoch line differs: {label_a}={:?} {label_b}={:?}",
                ea.flag,
                ea.epoch_line.trim_end(),
                eb.epoch_line.trim_end()
            ));
        } else if !ea
            .lines
            .iter()
            .map(|l| l.trim_end())
            .eq(eb.lines.iter().map(|l| l.trim_end()))
        {
            out.push(format!(
                "event[{k}] (flag {}) text differs: {label_a}={} lines {label_b}={} lines",
                ea.flag,
                ea.lines.len(),
                eb.lines.len()
            ));
        }
    }
}

/// Satellites listed in different epochs, including epochs where all of a
/// satellite's observations are blank.
fn compare_satellite_epochs(
    a: &ObsFile,
    b: &ObsFile,
    label_a: &str,
    label_b: &str,
    out: &mut Vec<String>,
) {
    let sats: BTreeSet<_> = a
        .satellite_epochs
        .keys()
        .chain(b.satellite_epochs.keys())
        .collect();
    for sat in sats {
        let none = Vec::new();
        let ea = a.satellite_epochs.get(sat).unwrap_or(&none);
        let eb = b.satellite_epochs.get(sat).unwrap_or(&none);
        let (mut ia, mut ib) = (0, 0);
        while out.len() < MAX_DIFFS && (ia < ea.len() || ib < eb.len()) {
            let (i, file, label) = match (ea.get(ia), eb.get(ib)) {
                (Some(x), Some(y)) if x == y => {
                    ia += 1;
                    ib += 1;
                    continue;
                }
                (Some(&x), y) if y.map_or(true, |&y| x < y) => {
                    ia += 1;
                    (x, a, label_a)
                }
                (_, Some(&y)) => {
                    ib += 1;
                    (y, b, label_b)
                }
                _ => unreachable!(),
            };
            let time = file.epochs.get(i).map_or_else(String::new, fmt_epoch);
            out.push(format!(
                "{sat} listed in epoch[{i}] {time}: only in {label}"
            ));
        }
    }
}

/// Event epoch lines compare field by field, so the leading zero that
/// `crx2rnx` leaves off values below one (`-.123`) does not matter.
fn epoch_lines_eq(a: &str, b: &str) -> bool {
    fn field(t: &str) -> String {
        if let Some(rest) = t.strip_prefix("0.") {
            format!(".{rest}")
        } else if let Some(rest) = t.strip_prefix("-0.") {
            format!("-.{rest}")
        } else {
            t.to_owned()
        }
    }
    a.split_whitespace()
        .map(field)
        .eq(b.split_whitespace().map(field))
}

fn compare_observations(
    a: &ObsFile,
    b: &ObsFile,
    label_a: &str,
    label_b: &str,
    out: &mut Vec<String>,
) {
    let keys_a: BTreeSet<_> = a.observations.keys().collect();
    let keys_b: BTreeSet<_> = b.observations.keys().collect();

    for key in keys_a.union(&keys_b) {
        if out.len() >= MAX_DIFFS {
            break;
        }
        match (a.observations.get(key), b.observations.get(key)) {
            (None, _) => out.push(format!("{}/{}: only in {label_b}", key.0, key.1.as_str())),
            (_, None) => out.push(format!("{}/{}: only in {label_a}", key.0, key.1.as_str())),
            (Some(sa), Some(sb)) => {
                compare_series(**key, sa, sb, a, b, label_a, label_b, out);
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn compare_series(
    key: (crate::model::SatId, crate::rinex_header::ObsCode),
    sa: &[(usize, Obs)],
    sb: &[(usize, Obs)],
    file_a: &ObsFile,
    file_b: &ObsFile,
    label_a: &str,
    label_b: &str,
    out: &mut Vec<String>,
) {
    let prefix = format!("{}/{}", key.0, key.1.as_str());
    let mut ia = 0;
    let mut ib = 0;

    while ia < sa.len() && ib < sb.len() && out.len() < MAX_DIFFS {
        let (ea, oa) = sa[ia];
        let (eb, ob) = sb[ib];

        match ea.cmp(&eb) {
            std::cmp::Ordering::Less => {
                out.push(format!(
                    "{prefix} epoch[{ea}] {}: only in {label_a}",
                    fmt_epoch(&file_a.epochs[ea])
                ));
                ia += 1;
            }
            std::cmp::Ordering::Greater => {
                out.push(format!(
                    "{prefix} epoch[{eb}] {}: only in {label_b}",
                    fmt_epoch(&file_b.epochs[eb])
                ));
                ib += 1;
            }
            std::cmp::Ordering::Equal => {
                if !obs_eq(oa, ob) {
                    out.push(format!(
                        "{prefix} epoch[{ea}] {}: mismatch: {label_a}({})  {label_b}({})",
                        fmt_epoch(&file_a.epochs[ea]),
                        fmt_obs(oa),
                        fmt_obs(ob),
                    ));
                }
                ia += 1;
                ib += 1;
            }
        }
    }

    while ia < sa.len() && out.len() < MAX_DIFFS {
        let (ea, _) = sa[ia];
        out.push(format!(
            "{prefix} epoch[{ea}] {}: only in {label_a}",
            fmt_epoch(&file_a.epochs[ea])
        ));
        ia += 1;
    }
    while ib < sb.len() && out.len() < MAX_DIFFS {
        let (eb, _) = sb[ib];
        out.push(format!(
            "{prefix} epoch[{eb}] {}: only in {label_b}",
            fmt_epoch(&file_b.epochs[eb])
        ));
        ib += 1;
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{EventRecord, HeaderText, ObsCode, SatId, SatSystem};

    #[test]
    fn fmt_epoch_output() {
        let e = EpochRecord {
            year: 2020,
            month: 7,
            day: 18,
            hour: 0,
            minute: 0,
            seconds_e7: 0,
            extra_ps: 0,
            flag: 0,
            clock_offset_e12: 0,
        };
        assert_eq!(fmt_epoch(&e), "2020-07-18 00:00:00.0000000");
        let e = EpochRecord {
            seconds_e7: 361_234_567,
            extra_ps: 1_234,
            ..e
        };
        assert_eq!(fmt_epoch(&e), "2020-07-18 00:00:36.123456701234");
    }

    #[test]
    fn obs_eq_lli_zero_vs_space() {
        let a = Obs {
            value: 100,
            lli: b'0',
            ssi: b' ',
        };
        let b = Obs { lli: b' ', ..a };
        assert!(obs_eq(a, a));
        assert!(!obs_eq(a, b));
    }

    fn file(events: Vec<(usize, EventRecord)>, value: Option<i64>) -> ObsFile {
        let epoch = EpochRecord {
            year: 2024,
            month: 10,
            day: 1,
            hour: 0,
            minute: 0,
            seconds_e7: 0,
            extra_ps: 0,
            flag: 0,
            clock_offset_e12: 0,
        };
        let key = (
            SatId::new(SatSystem::GPS, 1),
            ObsCode::from_v3("D1C").unwrap(),
        );
        let obs = |value| {
            vec![(
                0,
                Obs {
                    value,
                    lli: b' ',
                    ssi: b' ',
                },
            )]
        };
        ObsFile {
            header: HeaderText::from_text("END OF HEADER\n"),
            epochs: vec![epoch],
            events,
            observations: value.map(|v| (key, obs(v))).into_iter().collect(),
            // Listed even without a value: every field blank.
            satellite_epochs: [(key.0, vec![0])].into_iter().collect(),
        }
    }

    #[test]
    fn epoch_lines_compare_by_field() {
        let line = "> 2024 10 01 00 00  4.0000000  1  2       0.000000123856";
        assert!(epoch_lines_eq(
            line,
            "> 2024 10 01 00 00  4.0000000  1  2        .000000123856"
        ));
        assert!(epoch_lines_eq("  -0.5", "-.5  "));
        assert!(!epoch_lines_eq(
            line,
            "> 2024 10 01 00 00  4.0000000  1  2       0.000000123857"
        ));
    }

    #[test]
    fn satellite_epochs_are_compared() {
        let a = file(Vec::new(), None);
        let mut b = file(Vec::new(), None);
        b.satellite_epochs.clear();
        assert_eq!(
            compare(&a, &b, "a", "b"),
            ["G01 listed in epoch[0] 2024-10-01 00:00:00.0000000: only in a"]
        );
    }

    #[test]
    fn zero_differs_from_absent() {
        let zero = file(Vec::new(), Some(0));
        assert_eq!(compare(&zero, &zero, "a", "b"), [] as [String; 0]);
        let diffs = compare(&zero, &file(Vec::new(), None), "a", "b");
        assert_eq!(diffs, ["G01/D1C: only in a"]);
    }

    #[test]
    fn events_are_compared() {
        let event_at = |at, flag, epoch_line: &str, line: &str| {
            (
                at,
                EventRecord {
                    flag,
                    epoch_line: epoch_line.to_owned(),
                    lines: vec![line.to_owned()],
                },
            )
        };
        let epoch_line = "                            4  1";
        let event = |at, flag, line: &str| event_at(at, flag, epoch_line, line);
        let a = file(vec![event(1, 4, "COMMENT")], Some(1));
        let same = file(
            vec![event_at(1, 4, &format!("{epoch_line}  "), "COMMENT   ")],
            Some(1),
        );
        assert_eq!(compare(&a, &same, "a", "b"), [] as [String; 0]);
        for (b, want) in [
            (file(Vec::new(), Some(1)), "event count: a=1  b=0"),
            (
                file(vec![event(0, 4, "COMMENT")], Some(1)),
                "event[0] position: epoch[1] vs epoch[0]",
            ),
            (
                file(vec![event(1, 3, "COMMENT")], Some(1)),
                "event[0] flag: 4 vs 3",
            ),
            (
                file(vec![event_at(1, 4, "  4  2", "COMMENT")], Some(1)),
                "event[0] (flag 4) epoch line differs: \
                 a=\"                            4  1\" b=\"  4  2\"",
            ),
            (
                file(vec![event(1, 4, "COMMENTS")], Some(1)),
                "event[0] (flag 4) text differs: a=1 lines b=1 lines",
            ),
        ] {
            assert_eq!(compare(&a, &b, "a", "b"), [want]);
        }
    }
}
