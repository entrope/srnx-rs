//! Per-chunk and file-level digest verification, per `docs/Succinct_RINEX.md`
//! §"Digest Identifiers".
//!
//! Supported identifiers:
//!
//! | id | bytes | algorithm |
//! |----|-------|-----------|
//! |  0 |   0   | null      |
//! |  2 |   4   | CRC32C    |
//! |  4 |  16   | BLAKE3 truncated to 16 bytes |
//! |  5 |  32   | BLAKE3 (32-byte output)      |
//!
//! BLAKE3-16 is the first 16 bytes of the standard 32-byte output, which by
//! the BLAKE3 XOF construction equals `finalize_xof().fill(16)`.
//!
//! CRC32C bytes appear in little-endian order on disk (RFC 3720 §12.1
//! convention; verified against the per-chunk digest of the SRNX header in
//! `rnx-data/at012000.srnx`).

use thiserror::Error;

#[derive(Debug, Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum DigestError {
    #[error("unknown digest identifier {0}")]
    Unknown(u64),
    #[error("digest buffer too short: need {expected} bytes, have {actual}")]
    Truncated { expected: usize, actual: usize },
    #[error("digest mismatch")]
    Mismatch,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DigestId {
    Null,
    Crc32c,
    Blake3_16,
    Blake3_32,
}

impl DigestId {
    pub fn from_raw(id: u64) -> Result<Self, DigestError> {
        match id {
            0 => Ok(Self::Null),
            2 => Ok(Self::Crc32c),
            4 => Ok(Self::Blake3_16),
            5 => Ok(Self::Blake3_32),
            _ => Err(DigestError::Unknown(id)),
        }
    }

    #[must_use]
    pub fn byte_len(self) -> usize {
        match self {
            Self::Null => 0,
            Self::Crc32c => 4,
            Self::Blake3_16 => 16,
            Self::Blake3_32 => 32,
        }
    }

    /// The identifier written in the `SRNX` chunk.
    #[must_use]
    pub fn raw(self) -> u64 {
        match self {
            Self::Null => 0,
            Self::Crc32c => 2,
            Self::Blake3_16 => 4,
            Self::Blake3_32 => 5,
        }
    }

    /// Append the digest of `data` (`byte_len()` bytes) to `out`.
    pub fn append(self, data: &[u8], out: &mut Vec<u8>) {
        match self {
            Self::Null => {}
            Self::Crc32c => out.extend_from_slice(&crc32c::crc32c(data).to_le_bytes()),
            Self::Blake3_16 => out.extend_from_slice(&blake3::hash(data).as_bytes()[..16]),
            Self::Blake3_32 => out.extend_from_slice(blake3::hash(data).as_bytes()),
        }
    }

    /// Compute the digest of `data` and compare it to the first `byte_len()`
    /// bytes of `expected`.
    pub fn verify(self, data: &[u8], expected: &[u8]) -> Result<(), DigestError> {
        let want_len = self.byte_len();
        if expected.len() < want_len {
            return Err(DigestError::Truncated {
                expected: want_len,
                actual: expected.len(),
            });
        }
        let expected = &expected[..want_len];
        let ok = match self {
            Self::Null => true,
            Self::Crc32c => crc32c::crc32c(data).to_le_bytes() == *expected,
            Self::Blake3_16 => &blake3::hash(data).as_bytes()[..16] == expected,
            Self::Blake3_32 => blake3::hash(data).as_bytes() == expected,
        };
        if ok {
            Ok(())
        } else {
            Err(DigestError::Mismatch)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_raw_known_and_unknown() {
        assert_eq!(DigestId::from_raw(0).unwrap(), DigestId::Null);
        assert_eq!(DigestId::from_raw(2).unwrap(), DigestId::Crc32c);
        assert_eq!(DigestId::from_raw(4).unwrap(), DigestId::Blake3_16);
        assert_eq!(DigestId::from_raw(5).unwrap(), DigestId::Blake3_32);
        assert_eq!(DigestId::from_raw(1), Err(DigestError::Unknown(1)));
        assert_eq!(DigestId::from_raw(3), Err(DigestError::Unknown(3)));
        assert_eq!(DigestId::from_raw(99), Err(DigestError::Unknown(99)));
    }

    #[test]
    fn append_verifies_and_round_trips_raw_ids() {
        for id in [
            DigestId::Null,
            DigestId::Crc32c,
            DigestId::Blake3_16,
            DigestId::Blake3_32,
        ] {
            let mut out = Vec::new();
            id.append(b"chunk", &mut out);
            assert_eq!(out.len(), id.byte_len());
            assert_eq!(id.verify(b"chunk", &out), Ok(()));
            assert_eq!(DigestId::from_raw(id.raw()), Ok(id));
        }
    }

    #[test]
    fn byte_lengths() {
        assert_eq!(DigestId::Null.byte_len(), 0);
        assert_eq!(DigestId::Crc32c.byte_len(), 4);
        assert_eq!(DigestId::Blake3_16.byte_len(), 16);
        assert_eq!(DigestId::Blake3_32.byte_len(), 32);
    }

    #[test]
    fn null_digest_always_verifies() {
        assert!(DigestId::Null.verify(b"", b"").is_ok());
        assert!(DigestId::Null.verify(b"anything", b"").is_ok());
        // Excess "expected" bytes are ignored (length 0 means take 0 bytes).
        assert!(DigestId::Null.verify(b"x", b"junk").is_ok());
    }

    #[test]
    fn crc32c_standard_test_vector() {
        // CRC32C("123456789") = 0xE3069283 (Castagnoli check value).
        let want = 0xE306_9283_u32.to_le_bytes();
        assert!(DigestId::Crc32c.verify(b"123456789", &want).is_ok());

        let mut wrong = want;
        wrong[0] ^= 1;
        assert_eq!(
            DigestId::Crc32c.verify(b"123456789", &wrong),
            Err(DigestError::Mismatch)
        );
    }

    #[test]
    fn crc32c_truncated_buffer() {
        assert_eq!(
            DigestId::Crc32c.verify(b"x", b"\x00\x00\x00"),
            Err(DigestError::Truncated {
                expected: 4,
                actual: 3
            })
        );
    }

    /// Well-known BLAKE3 hash of the empty string, from the BLAKE3 reference
    /// test vectors.
    const BLAKE3_EMPTY: [u8; 32] = [
        0xaf, 0x13, 0x49, 0xb9, 0xf5, 0xf9, 0xa1, 0xa6, 0xa0, 0x40, 0x4d, 0xea, 0x36, 0xdc, 0xc9,
        0x49, 0x9b, 0xcb, 0x25, 0xc9, 0xad, 0xc1, 0x12, 0xb7, 0xcc, 0x9a, 0x93, 0xca, 0xe4, 0x1f,
        0x32, 0x62,
    ];

    #[test]
    fn blake3_32_empty_input() {
        assert!(DigestId::Blake3_32.verify(b"", &BLAKE3_EMPTY).is_ok());
        assert_eq!(
            DigestId::Blake3_32.verify(b"x", &BLAKE3_EMPTY),
            Err(DigestError::Mismatch)
        );
    }

    #[test]
    fn blake3_16_is_prefix_of_blake3_32() {
        assert!(DigestId::Blake3_16.verify(b"", &BLAKE3_EMPTY[..16]).is_ok());
        // Full 32-byte buffer also works (we ignore the trailing bytes).
        assert!(DigestId::Blake3_16.verify(b"", &BLAKE3_EMPTY).is_ok());
    }

    /// End-to-end check that our CRC32C wrapper agrees with a real SRNX file
    /// on both polynomial and on-disk byte order.
    ///
    /// The first chunk of `at01_v211_1hz.srnx` is `SRNX` with a 16-byte payload
    /// (single-byte ULEB128 length 0x10), followed by a 4-byte CRC32C of
    /// `FOURCC || length || payload`.
    #[test]
    fn crc32c_matches_real_srnx_chunk() {
        let bytes = crate::test_data("at01_v211_1hz.srnx");

        // 4-byte FOURCC + 1-byte ULEB128 length (0x10) + 16-byte payload = 21.
        let covered = &bytes[..21];
        let on_disk_digest = &bytes[21..25];
        assert_eq!(
            DigestId::Crc32c.verify(covered, on_disk_digest),
            Ok(()),
            "SRNX-chunk CRC32C mismatch — likely byte-order or polynomial bug"
        );
    }
}
