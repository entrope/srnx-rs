//! Parser for the `SDIR` (Satellite Directory) chunk's payload.
//!
//! Per `docs/Succinct_RINEX.md` §SDIR, the payload is:
//!
//! 1. ULEB128 file offset of the `EPOC` chunk.
//! 2. ULEB128 file offset of the first `EVTF` chunk (0 = no event chunks).
//! 3. A sequence of satellite directory entries, each consisting of:
//!    * three-byte satellite identifier (RINEX-style, e.g. `b"G01"`),
//!    * ULEB128 file offset of that satellite's `SATE` chunk.
//!
//! The number of entries is implicit: read until the payload is exhausted.

use thiserror::Error;

use crate::compat::split_first_chunk;
use crate::leb128::{read_uleb128, LebError};

#[derive(Debug, Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum SdirError {
    #[error("invalid LEB128 in SDIR chunk: {0}")]
    Leb(#[from] LebError),
    #[error("SDIR truncated reading {what}")]
    Truncated { what: &'static str },
    #[error("SDIR lists satellite {0:?} more than once")]
    DuplicateSatellite(String),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SatelliteEntry {
    pub sat_id: [u8; 3],
    pub sate_offset: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sdir {
    pub epoc_offset: u64,
    /// File offset of the first `EVTF` chunk; 0 means no event chunks exist.
    pub evtf_offset: u64,
    pub satellites: Vec<SatelliteEntry>,
}

pub fn parse(payload: &[u8]) -> Result<Sdir, SdirError> {
    let mut cursor: &[u8] = payload;
    let epoc_offset = read_uleb128(&mut cursor)?;
    let evtf_offset = read_uleb128(&mut cursor)?;
    let mut satellites = Vec::new();
    let mut seen = std::collections::HashSet::new();
    while !cursor.is_empty() {
        let (sat_id_arr, rest) = split_first_chunk::<3>(cursor).ok_or(SdirError::Truncated {
            what: "satellite id",
        })?;
        let sat_id = *sat_id_arr;
        cursor = rest;
        let sate_offset = read_uleb128(&mut cursor)?;
        // "A given satellite name MUST be used in at most one SATE chunk",
        // so a second entry would repeat its observations.
        if !seen.insert(sat_id) {
            return Err(SdirError::DuplicateSatellite(
                String::from_utf8_lossy(&sat_id).into_owned(),
            ));
        }
        satellites.push(SatelliteEntry {
            sat_id,
            sate_offset,
        });
    }
    Ok(Sdir {
        epoc_offset,
        evtf_offset,
        satellites,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn encode_uleb(mut v: u64) -> Vec<u8> {
        let mut out = Vec::new();
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

    fn build(epoc: u64, evtf: u64, sats: &[(&[u8; 3], u64)]) -> Vec<u8> {
        let mut out = encode_uleb(epoc);
        out.extend(encode_uleb(evtf));
        for &(sat_id, sate_off) in sats {
            out.extend_from_slice(sat_id);
            out.extend(encode_uleb(sate_off));
        }
        out
    }

    #[test]
    fn parse_empty_sdir() {
        // No satellites and no events.
        let bytes = build(0, 0, &[]);
        let sdir = parse(&bytes).unwrap();
        assert_eq!(sdir.epoc_offset, 0);
        assert_eq!(sdir.evtf_offset, 0);
        assert_eq!(sdir.satellites, []);
    }

    #[test]
    fn parse_with_satellites() {
        let bytes = build(
            42,
            0,
            &[
                (b"G01", 100),
                (b"G02", 200),
                (b"R10", 300),
                (b"E12", 25_668_000),
            ],
        );
        let sdir = parse(&bytes).unwrap();
        assert_eq!(sdir.epoc_offset, 42);
        assert_eq!(sdir.evtf_offset, 0);
        assert_eq!(sdir.satellites.len(), 4);
        assert_eq!(&sdir.satellites[0].sat_id, b"G01");
        assert_eq!(sdir.satellites[0].sate_offset, 100);
        assert_eq!(&sdir.satellites[3].sat_id, b"E12");
        assert_eq!(sdir.satellites[3].sate_offset, 25_668_000);
    }

    #[test]
    fn parse_rejects_duplicate_satellites() {
        let bytes = build(0, 0, &[(b"G01", 100), (b"G02", 200), (b"G01", 300)]);
        assert_eq!(
            parse(&bytes),
            Err(SdirError::DuplicateSatellite("G01".to_string()))
        );
    }

    #[test]
    fn parse_with_events_offset() {
        let bytes = build(10, 999, &[(b"G01", 50)]);
        let sdir = parse(&bytes).unwrap();
        assert_eq!(sdir.evtf_offset, 999);
    }

    #[test]
    fn parse_truncated_satellite_id() {
        // epoc + evtf, then "G0" — cut off mid-id.
        let mut bytes = encode_uleb(0);
        bytes.extend(encode_uleb(0));
        bytes.extend_from_slice(b"G0");
        assert_eq!(
            parse(&bytes),
            Err(SdirError::Truncated {
                what: "satellite id"
            })
        );
    }

    #[test]
    fn parse_truncated_offset() {
        // satellite id present but offset ULEB cut off (continuation bit set).
        let mut bytes = encode_uleb(0);
        bytes.extend(encode_uleb(0));
        bytes.extend_from_slice(b"G01");
        bytes.push(0x80); // says "more bytes follow" but none do
        assert!(matches!(
            parse(&bytes),
            Err(SdirError::Leb(LebError::Truncated))
        ));
    }
}
