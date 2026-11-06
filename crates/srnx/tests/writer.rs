//! The SRNX writer against the committed fixtures: every fixture written
//! with each option reads back the same, and is about as small as the
//! reference encoder's output.

// Fixture names are ours, so extensions are always lowercase.
#![allow(clippy::case_sensitive_file_extension_comparisons)]

mod common;

use srnx::file::SrnxFile;
use srnx::model::{EventRecord, Obs, ObsCode, ObsFile, SatId, SatSystem};
use srnx::reader;
use srnx::writer::{self, Digest, ObsChunk, WriteError, WriteOptions};

/// Each fixture stem, with its RINEX text read.
fn fixtures() -> Vec<(String, ObsFile)> {
    let mut stems: Vec<String> = std::fs::read_dir(common::data_dir())
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .filter(|name| {
            (name.ends_with(".rnx") || name.ends_with(".crx")) && !name.ends_with(".crx.rnx")
        })
        .map(|name| name.split('.').next().unwrap().to_owned())
        .collect();
    stems.sort();
    stems.dedup();
    assert!(stems.len() >= 10, "{stems:?}");
    stems
        .into_iter()
        .map(|stem| {
            let file = reader::read(&common::rinex_text(&stem)).unwrap();
            (stem, file)
        })
        .collect()
}

fn options(
    obs_chunk: ObsChunk,
    chunk_digest: Digest,
    file_digest: Digest,
    threads: usize,
) -> WriteOptions {
    let mut options = WriteOptions::default();
    options.obs_chunk = obs_chunk;
    options.chunk_digest = chunk_digest;
    options.file_digest = file_digest;
    options.threads = threads;
    options
}

fn write(file: &ObsFile, options: &WriteOptions) -> Vec<u8> {
    let mut out = Vec::new();
    writer::write_with(file, options, &mut out).unwrap();
    out
}

#[test]
fn fixtures_round_trip() {
    for (stem, source) in fixtures() {
        for obs_chunk in [ObsChunk::Soc2, ObsChunk::Socd] {
            let bytes = write(
                &source,
                &options(obs_chunk, Digest::Crc32c, Digest::Blake3_32, 0),
            );
            let label = format!("{stem} written with {obs_chunk:?}");
            let file = reader::read(&bytes).unwrap_or_else(|e| panic!("{label}: {e}"));
            common::assert_same(&source, &file, &stem, &label);

            // The revision says which chunks the file uses.
            let parsed = SrnxFile::parse(&bytes).unwrap();
            let counts = parsed.count_by_fourcc();
            let epo2 = counts.contains_key(b"EPO2");
            let want_minor = match (epo2, obs_chunk) {
                (true, _) => 2,
                (false, ObsChunk::Soc2) => 1,
                _ => 0,
            };
            assert_eq!(parsed.header.minor, want_minor, "{label}");
            assert_eq!(
                counts.contains_key(b"SOC2"),
                obs_chunk == ObsChunk::Soc2,
                "{label}"
            );
            assert_eq!(
                counts.contains_key(b"SOCD"),
                obs_chunk == ObsChunk::Socd,
                "{label}"
            );
            assert_eq!(
                epo2,
                source.epochs.iter().any(|e| e.extra_ps != 0),
                "{label}"
            );
        }
    }
}

#[test]
fn output_does_not_depend_on_threads_and_digests_verify() {
    let (stem, source) = fixtures()
        .into_iter()
        .find(|(stem, _)| stem == "balm_v304_1hz")
        .unwrap();
    let default = write(&source, &WriteOptions::default());
    assert_eq!(
        write(
            &source,
            &options(ObsChunk::Soc2, Digest::Crc32c, Digest::Blake3_32, 1)
        ),
        default
    );
    assert_eq!(
        write(
            &source,
            &options(ObsChunk::Soc2, Digest::Crc32c, Digest::Blake3_32, 3)
        ),
        default
    );
    for chunk_digest in [
        Digest::None,
        Digest::Crc32c,
        Digest::Blake3_16,
        Digest::Blake3_32,
    ] {
        for file_digest in [
            Digest::None,
            Digest::Crc32c,
            Digest::Blake3_16,
            Digest::Blake3_32,
        ] {
            let bytes = write(
                &source,
                &options(ObsChunk::Socd, chunk_digest, file_digest, 0),
            );
            let parsed = SrnxFile::parse(&bytes).unwrap();
            assert!(parsed.file_digest_ok, "{chunk_digest:?}/{file_digest:?}");
            let file = reader::read(&bytes).unwrap();
            common::assert_same(
                &source,
                &file,
                &stem,
                &format!("{chunk_digest:?}/{file_digest:?}"),
            );
        }
    }
}

/// The writer's files are within a few percent of `rnx2srnx`'s.
#[test]
fn sizes_are_close_to_the_reference_encoder() {
    for (stem, source) in fixtures() {
        for (obs_chunk, suffix) in [(ObsChunk::Soc2, ""), (ObsChunk::Socd, ".socd")] {
            let name = format!("{stem}{suffix}.srnx");
            let Ok(reference) = std::fs::read(common::data_dir().join(&name)) else {
                continue;
            };
            let ours = write(
                &source,
                &options(obs_chunk, Digest::Crc32c, Digest::Blake3_32, 0),
            );
            let (ours, theirs) = (ours.len(), reference.len());
            eprintln!("{name}: {ours} bytes, rnx2srnx {theirs}");
            assert!(
                ours * 100 <= theirs * 102,
                "{name}: {ours} bytes, rnx2srnx {theirs}"
            );
        }
    }
}

fn minimal() -> ObsFile {
    let (_, mut file) = fixtures()
        .into_iter()
        .find(|(stem, _)| stem == "n302_v211_30s")
        .unwrap();
    file.events.clear();
    file
}

fn invalid(file: &ObsFile) -> String {
    match writer::write(file, Vec::new()) {
        Err(WriteError::Invalid(message)) => message,
        other => panic!("expected an invalid model, got {other:?}"),
    }
}

#[test]
fn unrepresentable_models_are_rejected() {
    let base = minimal();
    let g01 = SatId::new(SatSystem::GPS, 1);
    let obs = Obs {
        value: 1000,
        lli: b' ',
        ssi: b' ',
    };

    let mut file = minimal();
    file.observations
        .insert((g01, ObsCode::from_v3("L1C").unwrap()), vec![(0, obs)]);
    assert!(
        invalid(&file).contains("L1C is not in the header"),
        "{}",
        invalid(&file)
    );

    let key = *base.observations.keys().next().unwrap();
    let mut file = minimal();
    file.observations.get_mut(&key).unwrap()[0].1.lli = b'x';
    assert!(invalid(&file).contains("LLI 'x'"));

    let mut file = minimal();
    let series = file.observations.get_mut(&key).unwrap();
    series.swap(0, 1);
    assert!(invalid(&file).contains("not strictly increasing"));

    let mut file = minimal();
    let n = file.epochs.len();
    file.observations.get_mut(&key).unwrap().push((n, obs));
    assert!(invalid(&file).contains("past the last epoch"));

    let mut file = minimal();
    file.epochs[1].flag = 6;
    assert!(invalid(&file).contains("epoch 1: flag 6"));

    let mut file = minimal();
    file.events.push((
        0,
        EventRecord {
            flag: 4,
            epoch_line: "                            3  1".to_owned(),
            lines: vec!["COMMENT".to_owned()],
        },
    ));
    assert!(invalid(&file).contains("event 0: epoch line"));

    let mut file = minimal();
    file.header = srnx::model::HeaderText::from_text("not a header\n");
    assert!(invalid(&file).contains("RINEX VERSION / TYPE"));
}
