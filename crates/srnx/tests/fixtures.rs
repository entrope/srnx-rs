//! Tests against the committed fixtures in `tests/data`.  They need no
//! external tools.
//!
//! The excerpts (`n302_v211_30s`, `mty2_v211_15s`, `at01_v211_1hz`,
//! `balm_v304_1hz`, `balm_v402_1hz_ps`) and every `.srnx` and `.scan` file
//! there were made by `make_fixtures.py` with the reference tools; see
//! `tests/data/README.md`.

// Fixture names are ours, so extensions are always lowercase.
#![allow(clippy::case_sensitive_file_extension_comparisons)]

mod common;

use std::collections::BTreeMap;

use srnx::file::SrnxFile;
use srnx::reader;

/// Fixture files, grouped by stem (the name up to its first `.`).  Golden
/// `crx2rnx` outputs (`*.crx.rnx`) are left out; `crinex.rs` checks those.
fn fixtures_by_stem() -> BTreeMap<String, Vec<String>> {
    let mut groups = BTreeMap::<String, Vec<String>>::new();
    for entry in std::fs::read_dir(common::data_dir()).unwrap() {
        let name = entry.unwrap().file_name().into_string().unwrap();
        let is_obs = [".rnx", ".crx", ".srnx"]
            .iter()
            .any(|ext| name.ends_with(ext));
        if !is_obs || name.ends_with(".crx.rnx") {
            continue;
        }
        let stem = name.split('.').next().unwrap().to_owned();
        groups.entry(stem).or_default().push(name);
    }
    for names in groups.values_mut() {
        names.sort();
    }
    groups
}

fn srnx_fixtures() -> Vec<String> {
    let names: Vec<String> = fixtures_by_stem()
        .into_values()
        .flatten()
        .filter(|name| name.ends_with(".srnx"))
        .collect();
    assert_eq!(names.len(), 15, "SRNX fixtures: {names:?}");
    names
}

/// Every format of a fixture (RINEX, Compact RINEX, and each SRNX variant)
/// reads to the same observations and events as its RINEX text.
#[test]
fn every_format_of_a_fixture_reads_the_same() {
    for (stem, names) in fixtures_by_stem() {
        let reference = reader::read(&common::rinex_text(&stem))
            .unwrap_or_else(|e| panic!("{stem} RINEX text: {e}"));
        for name in &names {
            let file = reader::read(&common::data(name)).unwrap_or_else(|e| panic!("{name}: {e}"));
            common::assert_same(&reference, &file, &format!("{stem} RINEX text"), name);
        }
    }
}

/// Chunk counts and SOCD histograms match the `srnx_scan` report saved
/// beside each SRNX fixture.
#[test]
fn srnx_fixtures_match_srnx_scan() {
    for name in srnx_fixtures() {
        let scan_name = format!("{}.scan", name.strip_suffix(".srnx").unwrap());
        let scan = String::from_utf8(common::data(&scan_name)).unwrap();
        common::assert_matches_scan(&common::data(&name), &scan, &name);
    }
}

/// `--socd` fixtures are revision 1.0, picosecond fixtures are revision 1.2
/// with an `EPO2` chunk, and the rest are revision 1.1.  Every digest checks.
#[test]
fn srnx_revisions_and_digests() {
    for name in srnx_fixtures() {
        let bytes = common::data(&name);
        let file = SrnxFile::parse(&bytes).unwrap_or_else(|e| panic!("{name}: {e}"));
        let picoseconds = name.contains("_ps.") || name.starts_with("rinex402_picoseconds.");
        let want_minor = if name.ends_with(".socd.srnx") {
            0
        } else if picoseconds {
            2
        } else {
            1
        };
        assert_eq!(file.header.minor, want_minor, "{name}: minor version");
        let counts = file.count_by_fourcc();
        let observation_chunk = if want_minor == 0 { b"SOCD" } else { b"SOC2" };
        assert!(
            counts.contains_key(observation_chunk),
            "{name}: no {observation_chunk:?}"
        );
        assert_eq!(counts.contains_key(b"EPO2"), picoseconds, "{name}: EPO2");
        assert_eq!(counts.contains_key(b"EPOC"), !picoseconds, "{name}: EPOC");
        assert!(file.srnx_digest_ok, "{name}: SRNX chunk digest");
        assert!(file.file_digest_ok, "{name}: file digest");
        assert!(
            file.chunks.iter().all(|c| c.digest_ok),
            "{name}: chunk digest"
        );
    }
}

#[test]
fn picosecond_epochs() {
    let file = reader::read(&common::data("rinex402_picoseconds.srnx")).unwrap();
    assert_eq!((file.epochs[0].seconds_e7, file.epochs[0].extra_ps), (1, 1));

    // `make_fixtures.py` gives epoch k (from 1) `k * 37123 % 100000` ps.
    let file = reader::read(&common::data("balm_v402_1hz_ps.srnx")).unwrap();
    assert_eq!(file.epochs.len(), 120);
    assert_eq!(
        (file.epochs[1].seconds_e7, file.epochs[1].extra_ps),
        (10_000_000, 74_246)
    );
}

#[test]
fn events() {
    // mty2_v211_15s spans one flag-4 event; every format puts it before
    // epoch 40 (`every_format_of_a_fixture_reads_the_same` checks that).
    let file = reader::read(&common::data("mty2_v211_15s.srnx")).unwrap();
    assert_eq!(file.events.len(), 1);
    assert_eq!(file.events[0].0, 40);
    assert_ne!(file.events[0].1.lines, [] as [String; 0]);

    // Flag-1 and flag-6 records are events, not epochs, in SRNX too.
    let file = reader::read(&common::data("rinex211_flag1_flag6.srnx")).unwrap();
    assert_eq!(file.epochs.len(), 2);
    assert_eq!(file.events.len(), 2);

    let file = reader::read(&common::data("n302_v211_30s.srnx")).unwrap();
    assert_eq!(file.events.len(), 0);
}

// ── Digest corruption ────────────────────────────────────────────────────────

/// Return the number of LEB128 bytes needed to encode `n`.
fn leb128_len(mut n: usize) -> usize {
    let mut len = 1;
    while n >= 128 {
        n >>= 7;
        len += 1;
    }
    len
}

/// Corrupt a byte deep inside a body chunk's payload; its CRC32C must fail.
#[test]
fn corrupt_body_chunk_detected() {
    let bytes = common::data("n302_v211_30s.srnx");
    let file = SrnxFile::parse(&bytes).expect("initial parse");
    let rhdr = file
        .chunks
        .iter()
        .find(|c| c.fourcc == *b"RHDR")
        .expect("RHDR present");

    // Payload starts at: file_offset + 4 (FOURCC) + leb128_len(payload.len()).
    let payload_start =
        usize::try_from(rhdr.file_offset).unwrap() + 4 + leb128_len(rhdr.payload.len());
    let corrupt_at = payload_start + 128; // well inside the RHDR payload

    let mut tampered = bytes.clone();
    tampered[corrupt_at] ^= 0xFF;

    let file2 = SrnxFile::parse(&tampered).expect("parse still succeeds despite corruption");
    assert!(
        file2
            .chunks
            .iter()
            .find(|c| c.fourcc == *b"RHDR")
            .is_some_and(|c| !c.digest_ok),
        "RHDR per-chunk digest should fail after payload corruption"
    );
}

/// Flip a bit in the SRNX header chunk's own per-chunk CRC32C;
/// `srnx_digest_ok` must become `false`.
#[test]
fn corrupt_srnx_chunk_digest_detected() {
    let bytes = common::data("n302_v211_30s.srnx");
    let file = SrnxFile::parse(&bytes).expect("initial parse");
    assert!(
        file.srnx_digest_ok,
        "pristine file must have valid SRNX digest"
    );

    // The SRNX chunk's CRC32C (4 bytes) sits immediately before the first
    // body chunk.
    let body_start =
        usize::try_from(file.chunks.first().expect("has body chunks").file_offset).unwrap();
    let mut tampered = bytes.clone();
    tampered[body_start - 4] ^= 0x01;

    let file2 = SrnxFile::parse(&tampered).expect("parse still succeeds");
    assert!(
        !file2.srnx_digest_ok,
        "SRNX chunk digest should fail after CRC corruption"
    );
}
