//! LEB128 codecs as specified in `docs/Succinct_RINEX.md` §"General Overview".
//!
//! Two flavours:
//!
//! * **ULEB128** — little-endian base-128, MSB of each non-terminal byte set.
//! * **SLEB128** — ULEB128 with **`ZigZag`** (Protocol Buffers) sign mapping
//!   *not* the DWARF sign-extending SLEB128. Odd encoded values are negative.
//!
//! Decoding policy:
//!
//! * Reject inputs that overflow `u64` / `i64`. The top of the value space is
//!   used by SRNX (e.g. `INT64_MIN` is the absent-observation sentinel), so
//!   silent wraparound would be unsafe.
//! * Accept overlong-but-in-range encodings (e.g. `0x80 0x00` for zero). The
//!   spec does not forbid them and well-formed encoders simply do not emit
//!   them; failing on them would be needlessly fragile.

use thiserror::Error;

#[derive(Debug, Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum LebError {
    #[error("input truncated while decoding LEB128")]
    Truncated,
    #[error("LEB128 value does not fit in 64 bits")]
    Overflow,
}

/// Decode a ULEB128 value, advancing `input` past the consumed bytes.
pub fn read_uleb128(input: &mut &[u8]) -> Result<u64, LebError> {
    let mut result: u64 = 0;
    let mut shift: u32 = 0;
    loop {
        let (&byte, rest) = input.split_first().ok_or(LebError::Truncated)?;
        *input = rest;
        let payload = u64::from(byte & 0x7F);
        let cont = byte & 0x80 != 0;

        if shift < 57 {
            // 7 payload bits always fit (max bit position = shift + 6 <= 62).
            result |= payload << shift;
        } else if shift < 64 {
            // 1..=6 bits of room left in `result`.
            let room = 64 - shift;
            if payload >> room != 0 {
                return Err(LebError::Overflow);
            }
            result |= payload << shift;
        } else {
            // No room left; only zero-payload continuation bytes are tolerated.
            if payload != 0 {
                return Err(LebError::Overflow);
            }
        }

        if !cont {
            return Ok(result);
        }
        shift = shift.saturating_add(7);
    }
}

/// Decode an SLEB128 (`ZigZag`) value, advancing `input` past the consumed bytes.
pub fn read_sleb128(input: &mut &[u8]) -> Result<i64, LebError> {
    let u = read_uleb128(input)?;
    // ZigZag decode: even -> non-negative, odd -> negative.
    #[allow(clippy::cast_possible_wrap)]
    Ok(((u >> 1) as i64) ^ -((u & 1) as i64))
}

/// Append the ULEB128 encoding of `value` to `out`.
pub(crate) fn write_uleb128(out: &mut Vec<u8>, mut value: u64) {
    while value >= 0x80 {
        out.push(u8::try_from(value & 0x7F).unwrap() | 0x80);
        value >>= 7;
    }
    out.push(u8::try_from(value).unwrap());
}

/// Append the SLEB128 (`ZigZag`) encoding of `value` to `out`.
pub(crate) fn write_sleb128(out: &mut Vec<u8>, value: i64) {
    write_uleb128(out, zigzag(value));
}

/// The `ZigZag` mapping of `value`: even for non-negative, odd for negative.
#[allow(clippy::cast_sign_loss)]
pub(crate) const fn zigzag(value: i64) -> u64 {
    ((value << 1) ^ (value >> 63)) as u64
}

/// Length in bytes of the ULEB128 encoding of `value`.
pub(crate) const fn uleb128_len(value: u64) -> usize {
    let bits = 64 - (value | 1).leading_zeros();
    bits.div_ceil(7) as usize
}

/// Length in bytes of the SLEB128 encoding of `value`.
pub(crate) const fn sleb128_len(value: i64) -> usize {
    uleb128_len(zigzag(value))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn decode_u(bytes: &[u8]) -> Result<(u64, usize), LebError> {
        let mut cursor = bytes;
        let v = read_uleb128(&mut cursor)?;
        Ok((v, bytes.len() - cursor.len()))
    }

    fn decode_s(bytes: &[u8]) -> Result<(i64, usize), LebError> {
        let mut cursor = bytes;
        let v = read_sleb128(&mut cursor)?;
        Ok((v, bytes.len() - cursor.len()))
    }

    #[test]
    fn writers_round_trip() {
        for v in [0, 1, 127, 128, 300, 16_383, 16_384, u64::MAX / 3, u64::MAX] {
            let mut out = Vec::new();
            write_uleb128(&mut out, v);
            assert_eq!(out.len(), uleb128_len(v), "{v}");
            assert_eq!(decode_u(&out), Ok((v, out.len())));
        }
        for v in [0, 1, -1, 63, -64, 64, -65, i64::MAX, i64::MIN] {
            let mut out = Vec::new();
            write_sleb128(&mut out, v);
            assert_eq!(out.len(), sleb128_len(v), "{v}");
            assert_eq!(decode_s(&out), Ok((v, out.len())));
        }
    }

    #[test]
    fn uleb128_canonical_encodings() {
        assert_eq!(decode_u(&[0x00]).unwrap(), (0, 1));
        assert_eq!(decode_u(&[0x01]).unwrap(), (1, 1));
        assert_eq!(decode_u(&[0x7F]).unwrap(), (127, 1));
        assert_eq!(decode_u(&[0x80, 0x01]).unwrap(), (128, 2));
        assert_eq!(decode_u(&[0xFF, 0x7F]).unwrap(), (16_383, 2));
        assert_eq!(decode_u(&[0x80, 0x80, 0x01]).unwrap(), (16_384, 3));
        // u64::MAX = 9 bytes of 0xFF + 0x01 (10 bytes total).
        assert_eq!(
            decode_u(&[0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x01]).unwrap(),
            (u64::MAX, 10)
        );
    }

    #[test]
    fn uleb128_stops_at_terminal_byte() {
        // Trailing bytes must be left in the cursor.
        assert_eq!(decode_u(&[0x05, 0xAA, 0xBB]).unwrap(), (5, 1));
    }

    #[test]
    fn uleb128_overlong_zero_is_accepted() {
        // 0x80 0x80 0x00 = overlong but in-range encoding of 0.
        assert_eq!(decode_u(&[0x80, 0x80, 0x00]).unwrap(), (0, 3));
    }

    #[test]
    fn uleb128_truncated() {
        assert_eq!(decode_u(&[]), Err(LebError::Truncated));
        assert_eq!(decode_u(&[0x80]), Err(LebError::Truncated));
        assert_eq!(
            decode_u(&[0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF]),
            Err(LebError::Truncated)
        );
    }

    #[test]
    fn uleb128_overflow_in_tenth_byte() {
        // 9 bytes of 0xFF (continuation), then 0x02 in the 10th: bit 1 of the
        // 10th byte would land at bit 64, which is out of range.
        assert_eq!(
            decode_u(&[0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x02]),
            Err(LebError::Overflow)
        );
    }

    #[test]
    fn uleb128_overflow_eleventh_byte_with_payload() {
        // 10 zero-payload continuation bytes then a nonzero terminal byte:
        // shift has advanced past 64, payload must be zero.
        assert_eq!(
            decode_u(&[0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x01]),
            Err(LebError::Overflow)
        );
    }

    #[test]
    fn sleb128_zigzag_small_values() {
        assert_eq!(decode_s(&[0x00]).unwrap(), (0, 1));
        assert_eq!(decode_s(&[0x01]).unwrap(), (-1, 1));
        assert_eq!(decode_s(&[0x02]).unwrap(), (1, 1));
        assert_eq!(decode_s(&[0x03]).unwrap(), (-2, 1));
        assert_eq!(decode_s(&[0x04]).unwrap(), (2, 1));
    }

    #[test]
    fn sleb128_extremes() {
        // i64::MAX zigzags to 2 * i64::MAX = u64::MAX - 1.
        assert_eq!(
            decode_s(&[0xFE, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x01]).unwrap(),
            (i64::MAX, 10)
        );
        // i64::MIN zigzags to u64::MAX.
        assert_eq!(
            decode_s(&[0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x01]).unwrap(),
            (i64::MIN, 10)
        );
    }

    #[test]
    fn sleb128_propagates_uleb_errors() {
        assert_eq!(decode_s(&[]), Err(LebError::Truncated));
        assert_eq!(
            decode_s(&[0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x02]),
            Err(LebError::Overflow)
        );
    }
}
