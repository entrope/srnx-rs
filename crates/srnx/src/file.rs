//! Top-level SRNX file reader.
//!
//! Owns the chicken-and-egg bootstrap that the chunk-framing layer cannot do
//! by itself: parse the `SRNX` header chunk to learn the per-chunk and
//! file-level digest IDs, verify both digests, and walk the rest of the file
//! recording every body chunk along with its absolute file offset. Per-chunk
//! payload parsing for `RHDR`, `EPOC`/`EPO2`, `SATE`, `SOCD`, `EVTF` is the job of
//! later phases — this module only frames them.

use std::collections::BTreeMap;

use thiserror::Error;

use crate::chunk::{read_chunk, read_framed, ChunkError};
use crate::digest::DigestError;
use crate::epoc::{self, Epoc, EpocError, EpochChunkKind};
use crate::evtf::{self, Evtf, EvtfError};
use crate::sate::{self, Sate, SateError};
use crate::sdir::{self, Sdir, SdirError};
use crate::soc2::{self, Soc2, Soc2Error};
use crate::socd::{self, Socd, SocdError};
use crate::srnx_header::{self, SrnxHeader, SrnxHeaderError};

/// FOURCC tag for the SRNX header chunk.
pub const SRNX_TAG: [u8; 4] = *b"SRNX";
/// FOURCC tag for the RINEX header chunk.
pub const RHDR_TAG: [u8; 4] = *b"RHDR";
/// FOURCC tag for the satellite directory chunk.
pub const SDIR_TAG: [u8; 4] = *b"SDIR";
/// FOURCC tag for the epoch records chunk.
pub const EPOC_TAG: [u8; 4] = *b"EPOC";
/// FOURCC tag for revision 1.2 epoch records with picosecond resolution.
pub const EPO2_TAG: [u8; 4] = *b"EPO2";
/// FOURCC tag for the satellite observations metadata chunk.
pub const SATE_TAG: [u8; 4] = *b"SATE";
/// FOURCC tag for the satellite observation code data chunk.
pub const SOCD_TAG: [u8; 4] = *b"SOCD";
/// FOURCC tag for revision 1.1 observation data.
pub const SOC2_TAG: [u8; 4] = *b"SOC2";
/// FOURCC tag for the event flag / special record chunk.
pub const EVTF_TAG: [u8; 4] = *b"EVTF";

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum FileError {
    #[error("chunk framing error: {0}")]
    Framing(#[from] ChunkError),
    #[error("SRNX header error: {0}")]
    Header(#[from] SrnxHeaderError),
    #[error("SDIR chunk error: {0}")]
    Sdir(#[from] SdirError),
    #[error("{} chunk error: {source}", String::from_utf8_lossy(&kind.fourcc()))]
    Epoc {
        kind: EpochChunkKind,
        source: EpocError,
    },
    #[error("file has both EPOC and EPO2 chunks")]
    MultipleEpochChunks,
    #[error("SATE chunk error: {0}")]
    Sate(#[from] SateError),
    #[error("SOCD chunk error: {0}")]
    Socd(#[from] SocdError),
    #[error("SOC2 chunk error: {0}")]
    Soc2(#[from] Soc2Error),
    #[error("observation chunk declares {declared} observations, but the satellite has {present} epochs")]
    ObservationCount { declared: u64, present: usize },
    #[error("EVTF chunk at file offset {offset}: {source}")]
    Evtf { offset: u64, source: EvtfError },
    #[error("SDIR gives EVTF offset {sdir}, but the first EVTF chunk is at {first}")]
    EvtfOffset { sdir: u64, first: u64 },
    #[error("first chunk is {got:?}, expected SRNX")]
    NotSrnxFile { got: [u8; 4] },
    #[error("file truncated reading {what}")]
    Truncated { what: &'static str },
}

/// One body chunk (i.e. anything after the `SRNX` header chunk and before
/// the trailing file-level digest), plus its absolute byte offset within
/// the input.
#[derive(Debug, Clone, Copy)]
pub struct BodyChunk<'a> {
    pub fourcc: [u8; 4],
    pub payload: &'a [u8],
    /// Result of per-chunk digest verification (always `true` under
    /// `DigestId::Null`).
    pub digest_ok: bool,
    /// Offset of the chunk's FOURCC tag from the start of the file.
    pub file_offset: u64,
}

#[derive(Debug)]
pub struct SrnxFile<'a> {
    pub header: SrnxHeader,
    /// Per-chunk digest verification of the `SRNX` chunk itself.
    pub srnx_digest_ok: bool,
    /// File-level digest verification result.
    pub file_digest_ok: bool,
    /// All body chunks in file order.
    pub chunks: Vec<BodyChunk<'a>>,
}

impl<'a> SrnxFile<'a> {
    pub fn parse(bytes: &'a [u8]) -> Result<Self, FileError> {
        // 1. Frame the SRNX chunk (digest length unknown at this point).
        let srnx_framed = read_framed(bytes)?;
        if srnx_framed.fourcc != SRNX_TAG {
            return Err(FileError::NotSrnxFile {
                got: srnx_framed.fourcc,
            });
        }
        let header = srnx_header::parse(srnx_framed.payload)?;

        // 2. Verify the SRNX chunk's per-chunk digest using the now-known ID.
        let chunk_dlen = header.chunk_digest.byte_len();
        if srnx_framed.after_payload.len() < chunk_dlen {
            return Err(FileError::Truncated {
                what: "SRNX chunk digest",
            });
        }
        let srnx_covered_len = bytes.len() - srnx_framed.after_payload.len();
        let srnx_covered = &bytes[..srnx_covered_len];
        let srnx_digest_bytes = &srnx_framed.after_payload[..chunk_dlen];
        let srnx_digest_ok = match header.chunk_digest.verify(srnx_covered, srnx_digest_bytes) {
            Ok(()) => true,
            Err(DigestError::Mismatch) => false,
            Err(DigestError::Truncated { .. } | DigestError::Unknown(_)) => unreachable!(),
        };

        // 3. Determine body bounds.
        let body_start = srnx_covered_len + chunk_dlen;
        let file_dlen = header.file_digest.byte_len();
        if bytes.len() < body_start + file_dlen {
            return Err(FileError::Truncated {
                what: "file-level digest",
            });
        }
        let body_end = bytes.len() - file_dlen;

        // 4. Verify file-level digest (covers entire file from the SRNX
        //    FOURCC until just before the file-level digest at the end).
        let file_digest_bytes = &bytes[body_end..];
        let file_digest_covered = &bytes[..body_end];
        let file_digest_ok = match header
            .file_digest
            .verify(file_digest_covered, file_digest_bytes)
        {
            Ok(()) => true,
            Err(DigestError::Mismatch) => false,
            Err(DigestError::Truncated { .. } | DigestError::Unknown(_)) => unreachable!(),
        };

        // 5. Walk body chunks, recording absolute file offsets.
        let body_data = &bytes[body_start..body_end];
        let mut chunks = Vec::new();
        let mut cursor = body_data;
        while !cursor.is_empty() {
            let consumed = body_data.len() - cursor.len();
            let file_offset = u64::try_from(body_start + consumed)
                .expect("usize fits in u64 on every supported target");
            let (chunk, rest) = read_chunk(cursor, header.chunk_digest)?;
            chunks.push(BodyChunk {
                fourcc: chunk.fourcc,
                payload: chunk.payload,
                digest_ok: chunk.digest_ok,
                file_offset,
            });
            cursor = rest;
        }

        Ok(Self {
            header,
            srnx_digest_ok,
            file_digest_ok,
            chunks,
        })
    }

    /// Locate the body chunk whose FOURCC tag begins at `offset`.
    ///
    /// Chunks are in increasing file order, so this is a binary search:
    /// the reader resolves an offset for every observation column, and a
    /// file has one chunk per (satellite, code) plus its metadata.
    #[must_use]
    pub fn chunk_at_offset(&self, offset: u64) -> Option<&BodyChunk<'a>> {
        let index = self.chunks.partition_point(|c| c.file_offset < offset);
        self.chunks.get(index).filter(|c| c.file_offset == offset)
    }

    /// Convenience: the payload of the first `RHDR` chunk, if any. The spec
    /// requires `RHDR` to be the second chunk of the file, but we don't
    /// enforce that here — callers can check chunk ordering themselves if
    /// they care.
    #[must_use]
    pub fn rhdr_payload(&self) -> Option<&'a [u8]> {
        self.chunks
            .iter()
            .find(|c| c.fourcc == RHDR_TAG)
            .map(|c| c.payload)
    }

    /// Parse the `SDIR` chunk if present. Looks it up by the offset stored
    /// in the SRNX header (rather than scanning by FOURCC) so that
    /// `header.sdir_offset == 0` cleanly maps to `Ok(None)`.
    pub fn sdir(&self) -> Result<Option<Sdir>, FileError> {
        if self.header.sdir_offset == 0 {
            return Ok(None);
        }
        let chunk = self
            .chunk_at_offset(self.header.sdir_offset)
            .ok_or(FileError::Truncated {
                what: "SDIR chunk at advertised offset",
            })?;
        if chunk.fourcc != SDIR_TAG {
            return Err(FileError::Truncated {
                what: "SDIR chunk at advertised offset (FOURCC mismatch)",
            });
        }
        Ok(Some(sdir::parse(chunk.payload)?))
    }

    /// Parse the epoch chunk (`EPOC` or `EPO2`) if present. Located via the
    /// `SDIR` chunk's `epoc_offset`, mirroring the SRNX-spec discovery
    /// mechanism. Returns `Ok(None)` when there is no SDIR or no epochs.
    pub fn epoc(&self) -> Result<Option<Epoc>, FileError> {
        let Some(sdir) = self.sdir()? else {
            return Ok(None);
        };
        if sdir.epoc_offset == 0 {
            return Ok(None);
        }
        if self.chunks.iter().any(|c| c.fourcc == EPOC_TAG)
            && self.chunks.iter().any(|c| c.fourcc == EPO2_TAG)
        {
            return Err(FileError::MultipleEpochChunks);
        }
        let chunk = self
            .chunk_at_offset(sdir.epoc_offset)
            .ok_or(FileError::Truncated {
                what: "epoch chunk at advertised offset",
            })?;
        let kind = match chunk.fourcc {
            EPOC_TAG => EpochChunkKind::Epoc,
            EPO2_TAG => EpochChunkKind::Epo2,
            _ => {
                return Err(FileError::Truncated {
                    what: "epoch chunk at advertised offset (FOURCC mismatch)",
                })
            }
        };
        let epoc =
            epoc::parse(chunk.payload, kind).map_err(|source| FileError::Epoc { kind, source })?;
        Ok(Some(epoc))
    }

    /// Parse a `SATE` chunk at the given absolute file offset.
    ///
    /// `num_obs_codes` must match the number of observation codes declared in
    /// the `RHDR` header; `total_epochs` is the count from the epoch chunk.
    pub fn sate_at(
        &self,
        offset: u64,
        num_obs_codes: usize,
        total_epochs: usize,
    ) -> Result<Sate, FileError> {
        let payload = self.sate_payload(offset)?;
        Ok(sate::parse(payload, num_obs_codes, total_epochs)?)
    }

    /// Parse only the observation chunk offsets of the `SATE` chunk at the
    /// given absolute file offset, skipping its epoch presence data.
    pub fn sate_offsets_at(
        &self,
        offset: u64,
        num_obs_codes: usize,
    ) -> Result<Vec<i64>, FileError> {
        let payload = self.sate_payload(offset)?;
        Ok(sate::parse_offsets(payload, num_obs_codes)?)
    }

    fn sate_payload(&self, offset: u64) -> Result<&'a [u8], FileError> {
        let chunk = self.chunk_at_offset(offset).ok_or(FileError::Truncated {
            what: "SATE chunk at given offset",
        })?;
        if chunk.fourcc != SATE_TAG {
            return Err(FileError::Truncated {
                what: "SATE chunk at given offset (FOURCC mismatch)",
            });
        }
        Ok(chunk.payload)
    }

    /// Parse a `SOCD` chunk at the given absolute file offset.
    pub fn socd_at(&self, offset: u64) -> Result<Socd, FileError> {
        let chunk = self.chunk_at_offset(offset).ok_or(FileError::Truncated {
            what: "SOCD chunk at given offset",
        })?;
        if chunk.fourcc != SOCD_TAG {
            return Err(FileError::Truncated {
                what: "SOCD chunk at given offset (FOURCC mismatch)",
            });
        }
        Ok(socd::parse(chunk.payload)?)
    }

    /// Parse SOC2 residuals and metadata at an absolute offset. Cross-signal
    /// references must subsequently be resolved with `soc2::SatelliteColumns`.
    pub fn soc2_at(&self, offset: u64) -> Result<Soc2, FileError> {
        let chunk = self.chunk_at_offset(offset).ok_or(FileError::Truncated {
            what: "SOC2 chunk at given offset",
        })?;
        if chunk.fourcc != SOC2_TAG {
            return Err(FileError::Truncated {
                what: "SOC2 FOURCC mismatch",
            });
        }
        Ok(soc2::parse(chunk.payload)?)
    }

    /// Parse every `EVTF` chunk, in file order, which keeps events at the
    /// same epoch index in their original order.
    ///
    /// The spec only says `EVTF` chunks SHOULD be stored consecutively, so
    /// this reads all of them rather than walking from the SDIR's
    /// `evtf_offset`, which must name the first one (or be zero if there
    /// are none).
    pub fn evtf_chunks(&self, evtf_offset: u64) -> Result<Vec<Evtf>, FileError> {
        let chunks: Vec<&BodyChunk<'a>> = self
            .chunks
            .iter()
            .filter(|c| c.fourcc == EVTF_TAG)
            .collect();
        let first = chunks.first().map_or(0, |c| c.file_offset);
        if first != evtf_offset {
            return Err(FileError::EvtfOffset {
                sdir: evtf_offset,
                first,
            });
        }
        chunks
            .into_iter()
            .map(|c| {
                evtf::parse(c.payload).map_err(|source| FileError::Evtf {
                    offset: c.file_offset,
                    source,
                })
            })
            .collect()
    }

    /// Histogram of body-chunk FOURCCs (does not include the SRNX chunk).
    #[must_use]
    pub fn count_by_fourcc(&self) -> BTreeMap<[u8; 4], usize> {
        let mut counts = BTreeMap::new();
        for c in &self.chunks {
            *counts.entry(c.fourcc).or_insert(0) += 1;
        }
        counts
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::digest::DigestId;

    /// Expected per-FOURCC body-chunk counts from `srnx_scan` reference output.
    struct Expected {
        rhdr: usize,
        epoc: usize,
        epo2: usize,
        sdir: usize,
        evtf: usize,
        sate: usize,
        socd: usize,
        soc2: usize,
    }

    fn check_counts(file: &SrnxFile<'_>, want: &Expected, name: &str) {
        let counts = file.count_by_fourcc();
        let g = |tag: &[u8; 4]| counts.get(tag).copied().unwrap_or(0);
        assert_eq!(g(b"RHDR"), want.rhdr, "{name}: RHDR count");
        assert_eq!(g(b"EPOC"), want.epoc, "{name}: EPOC count");
        assert_eq!(g(b"EPO2"), want.epo2, "{name}: EPO2 count");
        assert_eq!(g(b"SDIR"), want.sdir, "{name}: SDIR count");
        assert_eq!(g(b"EVTF"), want.evtf, "{name}: EVTF count");
        assert_eq!(g(b"SATE"), want.sate, "{name}: SATE count");
        assert_eq!(g(b"SOCD"), want.socd, "{name}: SOCD count");
        assert_eq!(g(b"SOC2"), want.soc2, "{name}: SOC2 count");
    }

    fn check_full(name: &str, want: &Expected) {
        let bytes = crate::test_data(name);
        let file = SrnxFile::parse(&bytes).expect("parse should succeed");
        assert!(file.srnx_digest_ok, "{name}: SRNX-chunk digest mismatch");
        assert!(file.file_digest_ok, "{name}: file-level digest mismatch");
        for chunk in &file.chunks {
            assert!(
                chunk.digest_ok,
                "{name}: digest mismatch on body chunk at file_offset={}",
                chunk.file_offset
            );
        }
        check_counts(&file, want, name);

        // SDIR lookup via offset must succeed and contain the expected number
        // of satellites (= count of distinct SATE chunks).
        let sdir = file
            .sdir()
            .expect("sdir() should not error")
            .expect("SDIR present");
        assert_eq!(
            sdir.satellites.len(),
            want.sate,
            "{name}: SDIR satellite count should match SATE chunk count"
        );

        // The epoch offset from SDIR must point at the EPOC or EPO2 chunk.
        let epoc_chunk = file
            .chunk_at_offset(sdir.epoc_offset)
            .expect("SDIR.epoc_offset should resolve");
        let want_fourcc = if want.epo2 == 1 { b"EPO2" } else { b"EPOC" };
        assert_eq!(&epoc_chunk.fourcc, want_fourcc);

        // Every SATE offset in SDIR must resolve to a SATE chunk.
        for entry in &sdir.satellites {
            let c = file.chunk_at_offset(entry.sate_offset).unwrap_or_else(|| {
                panic!(
                    "{name}: SATE offset {} for sat {:?} did not resolve",
                    entry.sate_offset, entry.sat_id
                )
            });
            assert_eq!(&c.fourcc, b"SATE");
        }

        // RHDR payload should be reachable.
        assert!(file.rhdr_payload().is_some(), "{name}: RHDR missing");
    }

    #[test]
    fn soc2_fixture() {
        check_full(
            "at01_v211_1hz.srnx",
            &Expected {
                rhdr: 1,
                epoc: 1,
                epo2: 0,
                sdir: 1,
                evtf: 0,
                sate: 29,
                socd: 0,
                soc2: 284,
            },
        );
    }

    #[test]
    fn socd_fixture_with_event() {
        check_full(
            "mty2_v211_15s.socd.srnx",
            &Expected {
                rhdr: 1,
                epoc: 1,
                epo2: 0,
                sdir: 1,
                evtf: 1,
                sate: 21,
                socd: 246,
                soc2: 0,
            },
        );
    }

    #[test]
    fn digest_variant_fixtures() {
        let want = Expected {
            rhdr: 1,
            epoc: 1,
            epo2: 0,
            sdir: 1,
            evtf: 0,
            sate: 22,
            socd: 0,
            soc2: 168,
        };
        for (name, chunk_digest, file_digest) in [
            ("n302_v211_30s.srnx", DigestId::Crc32c, DigestId::Blake3_32),
            ("n302_v211_30s.null.srnx", DigestId::Null, DigestId::Null),
            (
                "n302_v211_30s.blake3.srnx",
                DigestId::Blake3_32,
                DigestId::Blake3_16,
            ),
        ] {
            check_full(name, &want);
            let bytes = crate::test_data(name);
            let header = SrnxFile::parse(&bytes).unwrap().header;
            assert_eq!(
                (header.chunk_digest, header.file_digest),
                (chunk_digest, file_digest),
                "{name}"
            );
        }
    }

    #[test]
    fn epo2_fixture() {
        check_full(
            "balm_v402_1hz_ps.srnx",
            &Expected {
                rhdr: 1,
                epoc: 0,
                epo2: 1,
                sdir: 1,
                evtf: 0,
                sate: 47,
                socd: 0,
                soc2: 282,
            },
        );
    }

    #[test]
    fn rejects_non_srnx_file() {
        // 4-byte FOURCC + 1-byte length + 0-byte payload + 0 padding.
        let bytes = [b'F', b'A', b'K', b'E', 0x00];
        let err = SrnxFile::parse(&bytes).unwrap_err();
        assert!(matches!(err, FileError::NotSrnxFile { got } if &got == b"FAKE"));
    }

    #[test]
    fn chunk_at_offset_requires_an_exact_offset() {
        let bytes = crate::test_data("at01_v211_1hz.srnx");
        let file = SrnxFile::parse(&bytes).unwrap();
        // The SRNX chunk at offset 0 is framing, not a body chunk.
        assert!(file.chunk_at_offset(0).is_none());
        for chunk in &file.chunks {
            assert_eq!(
                file.chunk_at_offset(chunk.file_offset)
                    .map(|found| found.fourcc),
                Some(chunk.fourcc)
            );
            assert!(file.chunk_at_offset(chunk.file_offset + 1).is_none());
        }
        assert!(file.chunk_at_offset(u64::MAX).is_none());
    }

    #[test]
    fn detects_tampered_file_digest() {
        let mut bytes = crate::test_data("at01_v211_1hz.srnx");
        // Flip a bit in the last byte (file-level digest).
        let last = bytes.len() - 1;
        bytes[last] ^= 0x01;
        let file = SrnxFile::parse(&bytes).expect("parsing still succeeds");
        assert!(!file.file_digest_ok, "tampered file digest should fail");
    }
}
