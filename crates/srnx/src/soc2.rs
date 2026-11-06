//! Revision 1.1 observation decoding, including split Rice and cross-signal prediction.
//!
//! [`parse`] decodes residuals and metadata without references. [`SatelliteColumns`]
//! resolves value dependencies, then SSI dependencies, in RHDR observation-code order.
use std::collections::BTreeMap;

use thiserror::Error;

use crate::compat::{cast_signed, split_first_chunk};
use crate::leb128::{read_sleb128, read_uleb128, LebError};
use crate::rinex_header::ObsCode;
use crate::socd::{DeltaIntegrator, Socd, ABSENT};

#[derive(Debug, Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum Soc2Error {
    #[error("invalid LEB128 in SOC2: {0}")]
    Leb(#[from] LebError),
    #[error("truncated SOC2 {0}")]
    Truncated(&'static str),
    #[error("invalid SOC2: {0}")]
    Invalid(&'static str),
    #[error("SOC2 arithmetic overflow")]
    Overflow,
}

/// Value predictor; reference indices use RHDR/SATE order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Predictor {
    None,
    Value {
        reference: usize,
        num: i64,
        den: u64,
    },
    Rate {
        reference: usize,
        num: i64,
        den: u64,
        stencil: u8,
        max_span: u64,
    },
}

impl Predictor {
    fn reference(self) -> Option<usize> {
        match self {
            Self::None => None,
            Self::Value { reference, .. } | Self::Rate { reference, .. } => Some(reference),
        }
    }
}

/// Parsed residual column. `data.ssi` uses byte 11 for a predicted SSI until resolved.
#[derive(Debug, Clone)]
pub struct Soc2 {
    pub data: Socd,
    pub predictor: Predictor,
    pub ssi_reference: Option<usize>,
    /// Slots containing direct scaled observations instead of residuals.
    pub unpredicted: Vec<bool>,
}

fn index(v: u64) -> Result<usize, Soc2Error> {
    usize::try_from(v).map_err(|_| Soc2Error::Overflow)
}

fn section<'a>(cursor: &mut &'a [u8], what: &'static str) -> Result<&'a [u8], Soc2Error> {
    let len = index(read_uleb128(cursor)?)?;
    if len > cursor.len() {
        return Err(Soc2Error::Truncated(what));
    }
    let (data, rest) = cursor.split_at(len);
    *cursor = rest;
    Ok(data)
}

fn rle(cursor: &mut &[u8], count: usize, alphabet: u64) -> Result<Vec<u8>, Soc2Error> {
    let mut bytes = section(cursor, "indicator RLE")?;
    let mut out = Vec::new();
    while !bytes.is_empty() {
        let packed = read_uleb128(&mut bytes)?;
        let run = index(packed / alphabet)?
            .checked_add(1)
            .ok_or(Soc2Error::Overflow)?;
        if run > count - out.len() {
            return Err(Soc2Error::Invalid("indicator run exceeds count"));
        }
        let symbol = index(packed % alphabet)?;
        let ch = if symbol == 11 {
            11
        } else {
            b" 0123456789"[symbol]
        };
        out.resize(out.len() + run, ch);
    }
    out.resize(count, b' ');
    Ok(out)
}

fn initial(cursor: &mut &[u8], order: usize) -> Result<DeltaIntegrator, Soc2Error> {
    if order > 7 {
        return Err(Soc2Error::Invalid("delta order exceeds 7"));
    }
    let mut state = [0; 7];
    for item in &mut state[..order] {
        *item = read_sleb128(cursor)?;
    }
    Ok(DeltaIntegrator::new(order, &state[..order]))
}

/// Decode a SOC2 payload to scaled residuals. References are resolved separately.
#[allow(clippy::too_many_lines)]
pub fn parse(mut cursor: &[u8]) -> Result<Soc2, Soc2Error> {
    let (name, rest) =
        split_first_chunk::<8>(cursor).ok_or(Soc2Error::Truncated("observation name"))?;
    let obs_name = *name;
    cursor = rest;
    let count = index(
        read_uleb128(&mut cursor)?
            .checked_add(1)
            .ok_or(Soc2Error::Overflow)?,
    )?;
    let lli = rle(&mut cursor, count, 11)?;
    let mode = read_uleb128(&mut cursor)?;
    let ssi_reference = if mode >> 2 == 0 {
        None
    } else {
        Some(index((mode >> 2) - 1)?)
    };
    let kind = mode & 3;
    if kind == 3 {
        let _reserved_kind = read_uleb128(&mut cursor)?;
        return Err(Soc2Error::Invalid("reserved predictor kind"));
    }
    let ssi = rle(
        &mut cursor,
        count,
        if ssi_reference.is_some() { 12 } else { 11 },
    )?;
    let predictor = if kind == 0 {
        Predictor::None
    } else {
        let reference = index(read_uleb128(&mut cursor)?)?;
        let num = read_sleb128(&mut cursor)?;
        let den = read_uleb128(&mut cursor)?;
        if den == 0 {
            return Err(Soc2Error::Invalid("zero predictor denominator"));
        }
        if kind == 1 {
            Predictor::Value {
                reference,
                num,
                den,
            }
        } else {
            let stencil = read_uleb128(&mut cursor)?;
            if stencil > 2 {
                return Err(Soc2Error::Invalid("reserved rate stencil"));
            }
            Predictor::Rate {
                reference,
                num,
                den,
                stencil: u8::try_from(stencil).unwrap(),
                max_span: read_uleb128(&mut cursor)?,
            }
        }
    };
    let mut packed = section(&mut cursor, "packed data")?;
    if !cursor.is_empty() {
        return Err(Soc2Error::Invalid("trailing payload bytes"));
    }
    let schema = read_uleb128(&mut packed)?;
    if schema > 15 {
        return Err(Soc2Error::Invalid("reserved encoding schema"));
    }
    let delta_order = u8::try_from(schema & 7).unwrap();
    let scale = if schema >= 8 {
        let v = read_uleb128(&mut packed)?;
        if v == 0 || v > 1_000_000_000 {
            return Err(Soc2Error::Invalid("scale out of range"));
        }
        u32::try_from(v).unwrap()
    } else {
        1
    };
    let mut integrator = initial(&mut packed, usize::from(delta_order))?;
    let mut values = Vec::new();
    let mut unpredicted = Vec::new();
    let mut block_header_counts = BTreeMap::new();
    while let Some((&header, rest)) = packed.split_first() {
        packed = rest;
        *block_header_counts.entry(header).or_default() += 1;
        if header == 0xFB {
            let order = index(read_uleb128(&mut packed)?)?;
            integrator = initial(&mut packed, order)?;
            continue;
        }
        let run = match header {
            0x00..=0xDF => [8, 24, 56, 144, 240, 640, 2048][usize::from(header >> 5)],
            0xFC..=0xFF => index(
                read_uleb128(&mut packed)?
                    .checked_add(1)
                    .ok_or(Soc2Error::Overflow)?,
            )?,
            _ => return Err(Soc2Error::Invalid("reserved block header")),
        };
        if run > count - values.len() {
            return Err(Soc2Error::Invalid("block exceeds observation count"));
        }
        if header <= 0xDF {
            let r = usize::from(header & 31);
            let row_bytes = run / 8;
            let matrix_len = r * row_bytes;
            if packed.len() < matrix_len {
                return Err(Soc2Error::Truncated("Rice remainders"));
            }
            let (matrix, unary) = packed.split_at(matrix_len);
            let mut bit = 0usize;
            for i in 0..run {
                let mut remainder = 0u64;
                for row in 0..r {
                    remainder = (remainder << 1)
                        | u64::from((matrix[row * row_bytes + i / 8] >> (7 - i % 8)) & 1);
                }
                let mut quotient = 0u64;
                loop {
                    let byte = *unary
                        .get(bit / 8)
                        .ok_or(Soc2Error::Truncated("Rice quotients"))?;
                    let one = (byte >> (7 - bit % 8)) & 1;
                    bit += 1;
                    if one == 0 {
                        break;
                    }
                    quotient = quotient.checked_add(1).ok_or(Soc2Error::Overflow)?;
                    if quotient > (u64::MAX >> r) {
                        return Err(Soc2Error::Overflow);
                    }
                }
                let u = (quotient << r) | remainder;
                let delta = cast_signed(u >> 1) ^ -cast_signed(u & 1);
                values.push(integrator.integrate(delta));
            }
            if bit % 8 != 0 && unary[bit / 8] & ((1 << (8 - bit % 8)) - 1) != 0 {
                return Err(Soc2Error::Invalid("nonzero Rice padding"));
            }
            packed = &unary[bit.div_ceil(8)..];
        } else {
            if header == 0xFC && predictor == Predictor::None {
                return Err(Soc2Error::Invalid("unpredicted block without predictor"));
            }
            for _ in 0..run {
                values.push(match header {
                    0xFC => read_sleb128(&mut packed)?,
                    0xFD => ABSENT,
                    0xFE => integrator.integrate(0),
                    0xFF => integrator.integrate(read_sleb128(&mut packed)?),
                    _ => unreachable!(),
                });
            }
        }
        unpredicted.resize(values.len(), header == 0xFC);
    }
    if values.len() != count {
        return Err(Soc2Error::Invalid("observation count mismatch"));
    }
    Ok(Soc2 {
        data: Socd {
            obs_name,
            values,
            lli,
            ssi,
            scale,
            delta_order,
            block_header_counts,
        },
        predictor,
        ssi_reference,
        unpredicted,
    })
}

/// `num / den` rounded to the nearest integer, ties away from zero.
pub(crate) fn rounded(num: i128, den: i128) -> Result<i128, Soc2Error> {
    if den <= 0 {
        return Err(Soc2Error::Invalid("nonpositive predictor denominator"));
    }
    let q = num / den;
    let rem = num % den;
    if rem.abs() >= den / 2 + den % 2 {
        q.checked_add(num.signum()).ok_or(Soc2Error::Overflow)
    } else {
        Ok(q)
    }
}

/// The prediction `p[k]` in milli-units, or `None` if it is not available.
/// `milli(i)` is the reference value at satellite epoch `i`, if present.
pub(crate) fn prediction(
    p: Predictor,
    k: usize,
    milli: impl Fn(usize) -> Option<i128>,
    times: &[i128],
    units_per_second: u64,
) -> Result<Option<i128>, Soc2Error> {
    let (numerator, den) = match p {
        Predictor::None => return Ok(None),
        Predictor::Value { num, den, .. } => {
            let Some(value) = milli(k) else {
                return Ok(None);
            };
            (
                value
                    .checked_mul(i128::from(num))
                    .ok_or(Soc2Error::Overflow)?,
                i128::from(den),
            )
        }
        Predictor::Rate {
            num,
            den,
            stencil,
            max_span,
            ..
        } => {
            let pair = match stencil {
                0 => k.checked_sub(1).map(|a| (a, k)),
                1 => k.checked_sub(1).zip(k.checked_add(1)),
                2 => k.checked_add(1).map(|b| (k, b)),
                _ => return Err(Soc2Error::Invalid("reserved rate stencil")),
            };
            let Some((a, b)) = pair.filter(|&(_, b)| b < times.len()) else {
                return Ok(None);
            };
            let span = times[b].checked_sub(times[a]).ok_or(Soc2Error::Overflow)?;
            if span <= 0 || (max_span != 0 && span > i128::from(max_span)) {
                return Ok(None);
            }
            let Some((a, b)) = milli(a).zip(milli(b)) else {
                return Ok(None);
            };
            // A numerator that does not fit in i128 makes the prediction
            // unavailable rather than the chunk invalid.
            let Some(numerator) = (b - a)
                .checked_mul(i128::from(units_per_second))
                .and_then(|d| d.checked_mul(i128::from(num)))
            else {
                return Ok(None);
            };
            (
                numerator,
                span.checked_mul(i128::from(den))
                    .ok_or(Soc2Error::Overflow)?,
            )
        }
    };
    rounded(numerator, den).map(Some)
}

/// A satellite's SOCD/SOC2 columns in RHDR order; `None` denotes a zero SATE offset.
#[derive(Debug, Clone)]
pub enum Column {
    Socd(Socd),
    Soc2(Soc2),
}

impl Column {
    fn predictor(&self) -> Predictor {
        match self {
            Self::Socd(_) => Predictor::None,
            Self::Soc2(s) => s.predictor,
        }
    }
}

/// Resolve all of a satellite's observation values and indicators. Times
/// are continuous timestamps in units of `1 / units_per_second` seconds (the
/// epoch chunk's time unit), one per satellite epoch.
#[cfg(test)]
pub fn decode_satellite(
    columns: &[Option<Column>],
    codes: &[ObsCode],
    times: &[i128],
    units_per_second: u64,
) -> Result<Vec<Option<Socd>>, Soc2Error> {
    if columns.len() != codes.len() {
        return Err(Soc2Error::Invalid("observation code count mismatch"));
    }
    let mut cache = SatelliteColumns::new(codes, times, units_per_second);
    let mut load = |i: usize| Ok::<_, Soc2Error>(columns[i].clone());
    for i in 0..columns.len() {
        cache.decode(i, &mut load)?;
    }
    Ok(cache.into_columns())
}

/// One observation code's column in a [`SatelliteColumns`] cache.
enum Slot {
    Unloaded,
    /// Zero SATE offset: never observed.
    Missing,
    /// Parsed, with values not yet resolved.
    Parsed(Column),
    /// Values resolved.  `ssi_reference` is cleared once predicted SSIs are
    /// filled in; `requested` marks columns passed to `decode`.
    Decoded {
        data: Socd,
        predictor: Predictor,
        ssi_reference: Option<usize>,
        requested: bool,
    },
}

/// Lazily decoded observation columns of one satellite, indexed by the
/// position of their observation code in the RHDR list (the SATE order,
/// which predictor references also use).
///
/// Decoding a column also decodes the columns it references for value or
/// SSI prediction.  Every decoded column is kept, so a reference shared by
/// several columns is decoded once.  After an error the cache is
/// inconsistent and should be dropped.
pub struct SatelliteColumns<'a> {
    codes: &'a [ObsCode],
    times: &'a [i128],
    units_per_second: u64,
    slots: Vec<Slot>,
}

impl<'a> SatelliteColumns<'a> {
    /// `times` are continuous timestamps in units of `1 / units_per_second`
    /// seconds (the epoch chunk's time unit), one per satellite epoch.
    #[must_use]
    pub fn new(codes: &'a [ObsCode], times: &'a [i128], units_per_second: u64) -> Self {
        Self {
            codes,
            times,
            units_per_second,
            slots: codes.iter().map(|_| Slot::Unloaded).collect(),
        }
    }

    /// Decode column `i` with its values and indicators resolved, using
    /// `load(j)` to parse column `j` (`None` for a zero SATE offset) the
    /// first time it is needed.  Returns `None` for a zero SATE offset.
    ///
    /// # Panics
    ///
    /// If `i` is not less than the number of observation codes.
    pub fn decode<E: From<Soc2Error>>(
        &mut self,
        i: usize,
        load: &mut impl FnMut(usize) -> Result<Option<Column>, E>,
    ) -> Result<Option<&Socd>, E> {
        self.values(i, load)?;
        self.ssi(i, load)?;
        Ok(match &mut self.slots[i] {
            Slot::Decoded {
                data, requested, ..
            } => {
                *requested = true;
                Some(data)
            }
            _ => None,
        })
    }

    /// The columns passed to [`decode`](Self::decode), in RHDR order;
    /// `None` for the others and for zero SATE offsets.
    #[must_use]
    pub fn into_columns(self) -> Vec<Option<Socd>> {
        self.slots
            .into_iter()
            .map(|slot| match slot {
                Slot::Decoded {
                    data,
                    requested: true,
                    ..
                } => Some(data),
                _ => None,
            })
            .collect()
    }

    fn load<E>(
        &mut self,
        i: usize,
        load: &mut impl FnMut(usize) -> Result<Option<Column>, E>,
    ) -> Result<(), E> {
        if let Slot::Unloaded = self.slots[i] {
            self.slots[i] = load(i)?.map_or(Slot::Missing, Slot::Parsed);
        }
        Ok(())
    }

    /// Predictor of a loaded column.
    fn predictor(&self, i: usize) -> Predictor {
        match &self.slots[i] {
            Slot::Parsed(column) => column.predictor(),
            Slot::Decoded { predictor, .. } => *predictor,
            Slot::Unloaded | Slot::Missing => Predictor::None,
        }
    }

    fn decoded(&self, i: usize) -> Option<&Socd> {
        match &self.slots[i] {
            Slot::Decoded { data, .. } => Some(data),
            _ => None,
        }
    }

    /// Resolve column `i`'s values, after those of its predictor reference.
    fn values<E: From<Soc2Error>>(
        &mut self,
        i: usize,
        load: &mut impl FnMut(usize) -> Result<Option<Column>, E>,
    ) -> Result<(), E> {
        self.load(i, load)?;
        let Slot::Parsed(column) = &self.slots[i] else {
            return Ok(());
        };
        let n = self.times.len();
        let data = match column {
            Column::Socd(s) => s,
            Column::Soc2(s) => &s.data,
        };
        if data.values.len() != n {
            return Err(Soc2Error::Invalid("column count differs from satellite epochs").into());
        }
        if let Column::Soc2(s) = column {
            if s.unpredicted.len() != n || data.ssi.len() != n || data.scale == 0 {
                return Err(Soc2Error::Invalid("invalid residual column metadata").into());
            }
        }
        let predictor = column.predictor();
        if let Some(j) = predictor.reference() {
            if j == i {
                return Err(Soc2Error::Invalid("self predictor reference").into());
            }
            if j >= self.slots.len() {
                return Err(Soc2Error::Invalid("predictor index out of range").into());
            }
            self.load(j, load)?;
            // Checking the chain before recursing bounds the recursion.
            if let Some(base) = self.predictor(j).reference() {
                if base != i && base < self.slots.len() {
                    self.load(base, load)?;
                }
                if base == i || base >= self.slots.len() || self.predictor(base) != Predictor::None
                {
                    return Err(Soc2Error::Invalid(
                        "prediction chain exceeds two steps or is cyclic",
                    )
                    .into());
                }
            }
            self.values(j, load)?;
        }
        let Slot::Parsed(column) = std::mem::replace(&mut self.slots[i], Slot::Missing) else {
            unreachable!()
        };
        let reference = predictor.reference().and_then(|j| self.decoded(j));
        let (data, ssi_reference) = match column {
            Column::Socd(data) => (data, None),
            Column::Soc2(s) => {
                let ssi_reference = s.ssi_reference;
                (
                    apply_prediction(s, reference, self.times, self.units_per_second)?,
                    ssi_reference,
                )
            }
        };
        self.slots[i] = Slot::Decoded {
            data,
            predictor,
            ssi_reference,
            requested: false,
        };
        Ok(())
    }

    /// Fill in column `i`'s predicted SSIs.  SSI depends only on values, so
    /// SSI references may be cyclic.
    fn ssi<E: From<Soc2Error>>(
        &mut self,
        i: usize,
        load: &mut impl FnMut(usize) -> Result<Option<Column>, E>,
    ) -> Result<(), E> {
        let Slot::Decoded {
            data,
            ssi_reference: Some(j),
            ..
        } = &self.slots[i]
        else {
            return Ok(());
        };
        let j = *j;
        if j == i {
            return Err(Soc2Error::Invalid("self SSI reference").into());
        }
        if !self
            .codes
            .get(j)
            .is_some_and(|c| c.as_str().starts_with('S'))
        {
            return Err(Soc2Error::Invalid("SSI reference is not a signal-strength code").into());
        }
        if data.ssi.contains(&11) {
            self.values(j, load)?;
        }
        let Slot::Decoded { data, .. } = &self.slots[i] else {
            unreachable!()
        };
        let reference = self.decoded(j);
        let mut ssi = data.ssi.clone();
        for (k, ch) in ssi.iter_mut().enumerate() {
            if *ch == 11 {
                let reference =
                    reference.ok_or(Soc2Error::Invalid("predicted SSI reference absent"))?;
                let v = reference.values[k];
                if v == ABSENT {
                    return Err(Soc2Error::Invalid("predicted SSI reference absent").into());
                }
                let strength = i128::from(v) * i128::from(reference.scale);
                *ch = b'0' + u8::try_from(strength.div_euclid(6000).clamp(1, 9)).unwrap();
            }
        }
        if let Slot::Decoded {
            data,
            ssi_reference,
            ..
        } = &mut self.slots[i]
        {
            data.ssi = ssi;
            *ssi_reference = None;
        }
        Ok(())
    }
}

/// Add the predictions from `reference` (already resolved) to a column's
/// residuals.
fn apply_prediction(
    column: Soc2,
    reference: Option<&Socd>,
    times: &[i128],
    units_per_second: u64,
) -> Result<Socd, Soc2Error> {
    let Soc2 {
        mut data,
        predictor,
        unpredicted,
        ..
    } = column;
    if predictor.reference().is_none() {
        return Ok(data);
    }
    for (k, value) in data.values.iter_mut().enumerate() {
        if *value == ABSENT {
            continue;
        }
        let p = reference
            .map(|r| {
                prediction(
                    predictor,
                    k,
                    |i| reference_milli(r, i),
                    times,
                    units_per_second,
                )
            })
            .transpose()?
            .flatten();
        if unpredicted[k] != p.is_none() {
            return Err(Soc2Error::Invalid(
                "unpredicted block disagrees with prediction availability",
            ));
        }
        if let Some(p) = p {
            let q = rounded(p, i128::from(data.scale))?;
            *value = i64::try_from(
                i128::from(*value)
                    .checked_add(q)
                    .ok_or(Soc2Error::Overflow)?,
            )
            .map_err(|_| Soc2Error::Overflow)?;
            if *value == ABSENT {
                return Err(Soc2Error::Invalid(
                    "present observation equals absent sentinel",
                ));
            }
        }
    }
    Ok(data)
}

/// Value `i` of a decoded reference column in milli-units, if present.
fn reference_milli(reference: &Socd, i: usize) -> Option<i128> {
    reference
        .values
        .get(i)
        .filter(|&&v| v != ABSENT)
        .map(|&v| i128::from(v) * i128::from(reference.scale))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compat::cast_unsigned;

    /// `EPOC` time unit; `EPO2` is [`U12`].
    const U7: u64 = 10_000_000;
    const U12: u64 = 1_000_000_000_000;

    fn uleb(mut v: u64) -> Vec<u8> {
        let mut out = Vec::new();
        loop {
            let byte = u8::try_from(v & 127).unwrap();
            v >>= 7;
            out.push(byte | if v == 0 { 0 } else { 128 });
            if v == 0 {
                return out;
            }
        }
    }
    fn sleb(v: i64) -> Vec<u8> {
        uleb(cast_unsigned((v << 1) ^ (v >> 63)))
    }
    fn payload(count: u64, mode: u64, ssi: &[u8], fields: &[u8], packed: &[u8]) -> Vec<u8> {
        let mut out = b"G01\0L1\0\0".to_vec();
        out.extend(uleb(count - 1));
        out.push(0); // LLI
        out.extend(uleb(mode));
        out.extend(uleb(ssi.len().try_into().unwrap()));
        out.extend(ssi);
        out.extend(fields);
        out.extend(uleb(packed.len().try_into().unwrap()));
        out.extend(packed);
        out
    }
    fn data(values: &[i64], scale: u32) -> Socd {
        Socd {
            obs_name: *b"G01\0L1\0\0",
            values: values.to_vec(),
            lli: vec![b' '; values.len()],
            ssi: vec![b' '; values.len()],
            scale,
            delta_order: 0,
            block_header_counts: BTreeMap::new(),
        }
    }
    #[allow(clippy::unnecessary_wraps)] // Matches the optional SATE columns used by each test.
    fn column(
        values: &[i64],
        scale: u32,
        predictor: Predictor,
        unpredicted: &[bool],
    ) -> Option<Column> {
        Some(Column::Soc2(Soc2 {
            data: data(values, scale),
            predictor,
            ssi_reference: None,
            unpredicted: unpredicted.to_vec(),
        }))
    }
    fn codes(names: &[&str]) -> Vec<ObsCode> {
        names.iter().map(|s| ObsCode::from_v2(s).unwrap()).collect()
    }

    #[test]
    fn split_rice_all_sizes_and_parameters() {
        for (size_index, m) in [8usize, 24, 56, 144, 240, 640, 2048]
            .into_iter()
            .enumerate()
        {
            for r in [0usize, 1, 7, 31] {
                let unsigned: Vec<u64> = (0..m)
                    .map(|i| {
                        if r == 31 {
                            (1u64 << 31) + u64::try_from(i % 7).unwrap()
                        } else {
                            u64::try_from(i % 7).unwrap()
                        }
                    })
                    .collect();
                let mut packed = vec![0, u8::try_from((size_index << 5) | r).unwrap()];
                let mut matrix = vec![0u8; r * (m / 8)];
                for (i, &u) in unsigned.iter().enumerate() {
                    for row in 0..r {
                        matrix[row * (m / 8) + i / 8] |=
                            u8::try_from((u >> (r - row - 1)) & 1).unwrap() << (7 - i % 8);
                    }
                }
                packed.extend(matrix);
                let mut bits = Vec::new();
                for &u in &unsigned {
                    bits.extend(std::iter::repeat(true).take(usize::try_from(u >> r).unwrap()));
                    bits.push(false);
                }
                for byte_bits in bits.chunks(8) {
                    let mut byte = 0;
                    for (i, &one) in byte_bits.iter().enumerate() {
                        byte |= u8::from(one) << (7 - i);
                    }
                    packed.push(byte);
                }
                let s = parse(&payload(m.try_into().unwrap(), 0, &[], &[], &packed)).unwrap();
                let expected: Vec<i64> = unsigned
                    .iter()
                    .map(|&u| cast_signed(u >> 1) ^ -cast_signed(u & 1))
                    .collect();
                assert_eq!(s.data.values, expected, "m={m}, r={r}");
            }
        }
    }

    #[test]
    fn resets_and_absence_preserve_delta_state() {
        // Order 1 at 10, zero delta ->10, absent, reset to order 2 [20,2],
        // two zero deltas ->22,24; reset to order 0, direct delta -3.
        let mut packed = vec![1];
        packed.extend(sleb(10));
        packed.extend([0xFE, 0, 0xFD, 0, 0xFB, 2]);
        packed.extend(sleb(20));
        packed.extend(sleb(2));
        packed.extend([0xFE, 1, 0xFB, 0, 0xFF, 0]);
        packed.extend(sleb(-3));
        let s = parse(&payload(5, 0, &[], &[], &packed)).unwrap();
        assert_eq!(s.data.values, [10, ABSENT, 22, 24, -3]);
        assert_eq!(s.data.block_header_counts[&0xFB], 2);
    }

    #[test]
    fn value_prediction_rounding_scale_and_two_step_chain() {
        let columns = [
            column(
                &[1, -1, 0],
                2,
                Predictor::Value {
                    reference: 1,
                    num: 1,
                    den: 2,
                },
                &[false; 3],
            ),
            column(
                &[0, 0, 0],
                1,
                Predictor::Value {
                    reference: 2,
                    num: 1,
                    den: 2,
                },
                &[false; 3],
            ),
            Some(Column::Socd(data(&[6, -6, 2], 1))),
        ];
        let out = decode_satellite(&columns, &codes(&["C1", "L1", "L2"]), &[0, 1, 2], U7).unwrap();
        assert_eq!(out[1].as_ref().unwrap().values, [3, -3, 1]);
        // First rounding: +/-3/2 -> +/-2; scale rounding: +/-2/2 -> +/-1.
        assert_eq!(out[0].as_ref().unwrap().values, [2, -2, 1]);
    }

    #[test]
    fn lazy_decode_loads_only_requested_columns_and_references() {
        let columns = [
            column(
                &[3, 4],
                1,
                Predictor::Value {
                    reference: 2,
                    num: 1,
                    den: 1,
                },
                &[false; 2],
            ),
            Some(Column::Socd(data(&[7, 8], 1))),
            Some(Column::Socd(data(&[10, 20], 1))),
        ];
        let codes = codes(&["C1", "L2", "L1"]);
        let mut loaded = Vec::new();
        let mut cache = SatelliteColumns::new(&codes, &[0, 1], U7);
        let mut load = |i: usize| {
            loaded.push(i);
            Ok::<_, Soc2Error>(columns[i].clone())
        };
        let c1 = cache.decode(0, &mut load).unwrap().unwrap();
        assert_eq!(c1.values, [13, 24]);
        // The shared reference is decoded once.
        cache.decode(2, &mut load).unwrap();
        assert_eq!(loaded, [0, 2]);
        // Only requested columns are returned, not the unrequested L2.
        let out = cache.into_columns();
        assert_eq!(out[0].as_ref().unwrap().values, [13, 24]);
        assert!(out[1].is_none());
        assert_eq!(out[2].as_ref().unwrap().values, [10, 20]);
    }

    #[test]
    fn unpredicted_values_do_not_advance_integrator() {
        let mut fields = vec![1];
        fields.extend(sleb(1));
        fields.push(1);
        let mut packed = vec![1];
        packed.extend(sleb(5));
        packed.extend([0xFC, 0]);
        packed.extend(sleb(100));
        packed.extend([0xFE, 0, 0xFC, 0]);
        packed.extend(sleb(200));
        packed.extend([0xFE, 0]);
        let s = parse(&payload(4, 1, &[], &fields, &packed)).unwrap();
        let columns = [
            Some(Column::Soc2(s)),
            Some(Column::Socd(data(&[ABSENT, 10, ABSENT, 20], 1))),
        ];
        let out = decode_satellite(&columns, &codes(&["C1", "L1"]), &[0, 1, 2, 3], U7).unwrap();
        assert_eq!(out[0].as_ref().unwrap().values, [100, 15, 200, 25]);
    }

    #[test]
    fn rate_stencils_missing_values_and_max_span() {
        let reference = data(&[0, 30, 90, 180], 1);
        let times = [0, 10_000_000, 30_000_000, 60_000_000];
        for (stencil, expected) in [
            (0, [None, Some(-30), Some(-30), Some(-30)]),
            (1, [None, Some(-30), Some(-30), None]),
            (2, [Some(-30), Some(-30), Some(-30), None]),
        ] {
            let p = Predictor::Rate {
                reference: 0,
                num: -1,
                den: 1,
                stencil,
                max_span: 0,
            };
            for (k, exp) in expected.into_iter().enumerate() {
                assert_eq!(
                    prediction(p, k, |i| reference_milli(&reference, i), &times, U7).unwrap(),
                    exp
                );
            }
        }
        let p = Predictor::Rate {
            reference: 0,
            num: 1,
            den: 1,
            stencil: 0,
            max_span: 10_000_000,
        };
        assert_eq!(
            prediction(p, 1, |i| reference_milli(&reference, i), &times, U7).unwrap(),
            Some(30)
        );
        assert_eq!(
            prediction(p, 2, |i| reference_milli(&reference, i), &times, U7).unwrap(),
            None
        );
        let mut missing = reference.clone();
        missing.values[1] = ABSENT;
        assert_eq!(
            prediction(p, 1, |i| reference_milli(&missing, i), &times, U7).unwrap(),
            None
        );
        assert_eq!(
            prediction(p, 1, |i| reference_milli(&reference, i), &[0, 0, 1, 2], U7).unwrap(),
            None
        );
    }

    #[test]
    fn rate_prediction_uses_epoch_chunk_time_unit() {
        // Same 1 s and 2 s spans as `rate_stencils_…`, in EPO2 picoseconds.
        let reference = data(&[0, 30, 90], 1);
        let p = Predictor::Rate {
            reference: 0,
            num: 1,
            den: 1,
            stencil: 0,
            max_span: 0,
        };
        let times = [0, 1_000_000_000_000, 3_000_000_000_000];
        assert_eq!(
            prediction(p, 1, |i| reference_milli(&reference, i), &times, U12).unwrap(),
            Some(30)
        );
        assert_eq!(
            prediction(p, 2, |i| reference_milli(&reference, i), &times, U12).unwrap(),
            Some(30)
        );
    }

    #[test]
    fn rate_numerator_overflow_makes_prediction_unavailable() {
        // (R[b] - R[a]) ≈ 2^93 milli-units: × 10^7 fits in i128, × 10^12 does not.
        let reference = data(&[0, i64::MAX], 1_000_000_000);
        let p = Predictor::Rate {
            reference: 0,
            num: 1,
            den: 1,
            stencil: 0,
            max_span: 0,
        };
        assert!(prediction(
            p,
            1,
            |i| reference_milli(&reference, i),
            &[0, U7.into()],
            U7
        )
        .unwrap()
        .is_some());
        assert_eq!(
            prediction(
                p,
                1,
                |i| reference_milli(&reference, i),
                &[0, U12.into()],
                U12
            )
            .unwrap(),
            None
        );
    }

    #[test]
    fn predicted_ssi_uses_milliunits_and_clamps() {
        let reference = Some(Column::Socd(data(
            &[-1, 0, 5999, 6000, 12000, 54000, 60000],
            1,
        )));
        let mut s = Soc2 {
            data: data(&[0; 7], 1),
            predictor: Predictor::None,
            ssi_reference: Some(1),
            unpredicted: vec![false; 7],
        };
        s.data.ssi = vec![11; 7];
        let out = decode_satellite(
            &[Some(Column::Soc2(s.clone())), reference],
            &codes(&["L1", "S1"]),
            &[0; 7],
            U7,
        )
        .unwrap();
        assert_eq!(out[0].as_ref().unwrap().ssi, b"1111299");
        assert!(decode_satellite(
            &[
                Some(Column::Soc2(s)),
                Some(Column::Socd(data(&[ABSENT; 7], 1)))
            ],
            &codes(&["L1", "S1"]),
            &[0; 7],
            U7,
        )
        .is_err());
        let s = parse(&payload(3, 8, &[35], &[], &[0, 0xFE, 2])).unwrap();
        assert_eq!(s.data.ssi, [11; 3]); // twelve-symbol RLE
    }

    #[test]
    fn rejects_bad_dependencies_and_prediction_availability() {
        let p = |j| Predictor::Value {
            reference: j,
            num: 1,
            den: 1,
        };
        for columns in [
            vec![column(&[0], 1, p(0), &[false])],
            vec![
                column(&[0], 1, p(1), &[false]),
                column(&[0], 1, p(0), &[false]),
            ],
            vec![column(&[0], 1, p(1), &[false]), None],
            vec![
                column(&[0], 1, p(1), &[false]),
                Some(Column::Socd(data(&[ABSENT], 1))),
            ],
            vec![
                column(&[0], 1, p(1), &[true]),
                Some(Column::Socd(data(&[1], 1))),
            ],
            vec![
                column(&[0], 1, p(1), &[false]),
                column(&[0], 1, p(2), &[false]),
                column(&[0], 1, p(3), &[false]),
                Some(Column::Socd(data(&[1], 1))),
            ],
        ] {
            let c = vec![ObsCode::from_v2("L1").unwrap(); columns.len()];
            assert!(decode_satellite(&columns, &c, &[0], U7).is_err());
        }
    }

    #[test]
    fn absent_reference_columns_allow_only_unpredicted_values() {
        let p = Predictor::Value {
            reference: 1,
            num: 1,
            den: 1,
        };
        let out = decode_satellite(
            &[column(&[50, 60], 1, p, &[true; 2]), None],
            &codes(&["C1", "L1"]),
            &[0, 1],
            U7,
        )
        .unwrap();
        assert_eq!(out[0].as_ref().unwrap().values, [50, 60]);
    }

    #[test]
    fn central_rate_prediction_resolves_scaled_reference_and_endpoints() {
        let p = Predictor::Rate {
            reference: 1,
            num: -1,
            den: 2,
            stencil: 1,
            max_span: 0,
        };
        let out = decode_satellite(
            &[
                column(&[100, 2, 3, 200], 1, p, &[true, false, false, true]),
                Some(Column::Socd(data(&[0, 30, 90, 180], 2))),
            ],
            &codes(&["D1", "L1"]),
            &[0, 10_000_000, 30_000_000, 60_000_000],
            U7,
        )
        .unwrap();
        assert_eq!(out[0].as_ref().unwrap().values, [100, -28, -27, 200]);
    }

    #[test]
    fn arithmetic_uses_i128_and_reports_overflow() {
        let reference = data(&[i64::MAX], 1);
        let p = Predictor::Value {
            reference: 0,
            num: i64::MAX,
            den: cast_unsigned(i64::MAX),
        };
        assert_eq!(
            prediction(p, 0, |i| reference_milli(&reference, i), &[0], U7).unwrap(),
            Some(i128::from(i64::MAX))
        );
        let reference = data(&[i64::MAX], 1_000_000_000);
        assert_eq!(
            prediction(p, 0, |i| reference_milli(&reference, i), &[0], U7),
            Err(Soc2Error::Overflow)
        );
        let columns = [
            column(
                &[1],
                1,
                Predictor::Value {
                    reference: 1,
                    num: 1,
                    den: 1,
                },
                &[false],
            ),
            Some(Column::Socd(data(&[i64::MAX], 1))),
        ];
        assert!(matches!(
            decode_satellite(&columns, &codes(&["C1", "L1"]), &[0], U7),
            Err(Soc2Error::Overflow)
        ));
    }

    #[test]
    fn ssi_dependencies_do_not_create_value_cycles() {
        let mut a = Soc2 {
            data: data(&[12000], 1),
            predictor: Predictor::None,
            ssi_reference: Some(1),
            unpredicted: vec![false],
        };
        a.data.ssi = vec![11];
        let mut b = a.clone();
        b.ssi_reference = Some(0);
        b.data.values = vec![18000];
        let out = decode_satellite(
            &[Some(Column::Soc2(a)), Some(Column::Soc2(b))],
            &codes(&["S1", "S2"]),
            &[0],
            U7,
        )
        .unwrap();
        assert_eq!(out[0].as_ref().unwrap().ssi, b"3");
        assert_eq!(out[1].as_ref().unwrap().ssi, b"2");
    }

    #[test]
    fn malformed_payloads_return_errors() {
        for packed in [
            &[0, 0xE0][..],
            &[0, 0xFA],
            &[0, 0xFB, 8],
            &[0, 0xFC, 0, 0],
            &[0, 0xFD, 1],
            &[0, 0],
            &[16],
            &[8, 0],
            &[0, 0xFF, 0],
        ] {
            assert!(
                parse(&payload(1, 0, &[], &[], packed)).is_err(),
                "{packed:?}"
            );
        }
        // r=0 with eight terminating bits followed by one nonzero padding bit.
        assert!(parse(&payload(8, 0, &[], &[], &[0, 0, 0x80, 0x01])).is_err());
        assert!(parse(&payload(1, 3, &[], &[], &[0])).is_err());
        assert!(parse(&payload(1, 1, &[], &[0, 2, 0], &[0])).is_err());
        assert!(parse(&payload(1, 0, &[11], &[], &[0])).is_err());
    }
}
