//! `SOC2` predictor planning for the SRNX writer, following the spec's
//! encoder guidance as the reference encoder (`rnx2srnx`) does.
//!
//! Each satellite's reference carrier phase (the one present at the most
//! epochs) is coded without prediction.  Other carrier phases are predicted
//! from it, and pseudoranges and Doppler from the carrier phase on their
//! own band (or the reference).  The writer uses a planned predictor, or an
//! SSI reference, only where it makes the chunk smaller.

use crate::model::{ObsCode, SatId, SatSystem};
use crate::soc2::{rounded, Predictor};
use crate::socd::ABSENT;

/// Doppler predictor `max_span`, in seconds.
const DOPPLER_MAX_SPAN: u64 = 300;
const SPEED_OF_LIGHT: i64 = 299_792_458;

/// The planned coding of one observation code of a satellite.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Plan {
    /// The predictor to try, besides none.
    pub predictor: Predictor,
    /// A signal-strength code on the same band to try as SSI reference.
    pub ssi_reference: Option<usize>,
}

/// GLONASS frequency channels by slot, from the header's
/// `GLONASS SLOT / FRQ #` records.
pub(crate) fn glonass_channels(header: &str) -> [Option<i8>; 100] {
    let mut channels = [None; 100];
    for line in header.lines() {
        if line.get(60..80) != Some("GLONASS SLOT / FRQ #") {
            continue;
        }
        let content = line.as_bytes();
        let mut i = 0;
        while i + 3 <= 60 {
            if content[i] == b'R'
                && content[i + 1].is_ascii_digit()
                && content[i + 2].is_ascii_digit()
            {
                let slot = usize::from((content[i + 1] - b'0') * 10 + content[i + 2] - b'0');
                let channel = line
                    .get(i + 3..i + 7)
                    .and_then(leading_int)
                    .filter(|c| (-7..=13).contains(c));
                if channel.is_some() {
                    channels[slot] = channel;
                }
                i += 7;
            } else {
                i += 1;
            }
        }
    }
    channels
}

/// The integer at the start of `field`, after any blanks.
fn leading_int(field: &str) -> Option<i8> {
    let field = field.trim_start();
    let sign = usize::from(field.starts_with(['-', '+']));
    let end = sign + field[sign..].bytes().take_while(u8::is_ascii_digit).count();
    field[..end].parse().ok()
}

/// Carrier frequency in Hz of `band` (`b'1'`, …) of `system`, if known.
/// `channel` is the GLONASS FDMA frequency channel.
fn frequency(system: SatSystem, band: u8, channel: i64) -> Option<i64> {
    let hz = match (system.0, band) {
        (b'G' | b'J' | b'S' | b'I' | b'E' | b'C', b'1') => 1_575_420_000,
        (b'G' | b'J', b'2') => 1_227_600_000,
        (b'G' | b'J' | b'S' | b'I' | b'E' | b'C', b'5') => 1_176_450_000,
        (b'J' | b'E', b'6') => 1_278_750_000,
        (b'I', b'9') => 2_492_028_000,
        (b'E' | b'C', b'7') => 1_207_140_000,
        (b'E' | b'C', b'8') => 1_191_795_000,
        (b'C', b'2') => 1_561_098_000,
        (b'C', b'6') => 1_268_520_000,
        (b'R', b'1') => 1_602_000_000 + channel * 562_500,
        (b'R', b'2') => 1_246_000_000 + channel * 437_500,
        (b'R', b'3') => 1_202_025_000,
        (b'R', b'4') => 1_600_995_000,
        (b'R', b'6') => 1_248_060_000,
        _ => return None,
    };
    Some(hz)
}

/// Byte `i` of an observation code: type, band, attribute (0 for RINEX 2).
fn byte(code: ObsCode, i: usize) -> u8 {
    code.as_str().as_bytes().get(i).copied().unwrap_or(0)
}

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

/// `num / den` in lowest terms, with `den` positive.
fn ratio(num: i64, den: i64) -> (i64, u64) {
    let div = gcd(num.unsigned_abs(), den.unsigned_abs()).max(1);
    let div_i = i64::try_from(div).unwrap();
    (num / div_i, (den / div_i).unsigned_abs())
}

/// One satellite's columns, as the planner sees them.
struct Columns<'a> {
    codes: &'a [ObsCode],
    /// Values in milli-units, `ABSENT` where missing; `None` if never
    /// observed.
    values: &'a [Option<&'a [i64]>],
    /// Present values in each column.
    present: Vec<usize>,
}

impl Columns<'_> {
    /// The observed code of `kind` (`b'L'`, `b'S'`, …) on the band of
    /// `code`, other than `exclude`: one with the same attribute if any,
    /// else the one present at the most epochs.
    fn same_band(&self, kind: u8, code: ObsCode, exclude: Option<usize>) -> Option<usize> {
        let mut best: Option<usize> = None;
        for (j, &c) in self.codes.iter().enumerate() {
            if self.values[j].is_none()
                || Some(j) == exclude
                || byte(c, 0) != kind
                || byte(c, 1) != byte(code, 1)
            {
                continue;
            }
            if byte(c, 2) == byte(code, 2) {
                return Some(j);
            }
            if best.map_or(true, |b| self.present[j] > self.present[b]) {
                best = Some(j);
            }
        }
        best
    }

    /// The first observed code of `kinds` on band 1.
    fn first_band1(&self, kinds: &[u8]) -> Option<usize> {
        (0..self.codes.len()).find(|&j| {
            let c = self.codes[j];
            self.values[j].is_some() && kinds.contains(&byte(c, 0)) && byte(c, 1) == b'1'
        })
    }
}

/// Estimate a GLONASS frequency channel as the one making the band-1
/// code-minus-carrier sequence smoothest (least total |second difference|).
fn estimate_glonass_channel(code: &[i64], phase: &[i64]) -> i64 {
    let mut best = (0, f64::INFINITY);
    for channel in -7..=6 {
        let wavelength_den = i128::from(frequency(SatSystem::GLONASS, b'1', channel).unwrap());
        let (mut prev1, mut prev2): (Option<i128>, Option<i128>) = (None, None);
        let mut sum = 0.0;
        for (&c, &l) in code.iter().zip(phase) {
            if c == ABSENT || l == ABSENT {
                (prev1, prev2) = (None, None);
                continue;
            }
            let range =
                rounded(i128::from(l) * i128::from(SPEED_OF_LIGHT), wavelength_den).unwrap_or(0);
            let cmc = i128::from(c) - range;
            if let (Some(p1), Some(p2)) = (prev1, prev2) {
                #[allow(clippy::cast_precision_loss)]
                let d2 = (cmc - 2 * p1 + p2).unsigned_abs() as f64;
                sum += d2;
            }
            (prev2, prev1) = (prev1, Some(cmc));
        }
        if sum < best.1 {
            best = (channel, sum);
        }
    }
    best.0
}

/// Plan the `SOC2` coding of each observation code of `sat`.  `values`
/// are its columns over its satellite epochs; `units_per_second` is the
/// epoch chunk's time unit.
pub(crate) fn plan_satellite(
    sat: SatId,
    codes: &[ObsCode],
    values: &[Option<&[i64]>],
    glonass: &[Option<i8>; 100],
    units_per_second: u64,
) -> Vec<Plan> {
    let columns = Columns {
        codes,
        values,
        present: values
            .iter()
            .map(|v| v.map_or(0, |v| v.iter().filter(|&&x| x != ABSENT).count()))
            .collect(),
    };
    let system = sat.system;
    let mut channel = 0;
    if system == SatSystem::GLONASS {
        if let Some(c) = glonass
            .get(usize::from(sat.prn))
            .copied()
            .flatten()
            .filter(|_| sat.prn > 0)
        {
            channel = i64::from(c);
        } else if let (Some(c), Some(l)) = (columns.first_band1(b"CP"), columns.first_band1(b"L")) {
            channel = estimate_glonass_channel(values[c].unwrap(), values[l].unwrap());
        }
    }
    let freq = |code: ObsCode| frequency(system, byte(code, 1), channel);

    // Reference carrier phase: the one present at the most epochs.
    let mut reference: Option<usize> = None;
    for (j, &code) in codes.iter().enumerate() {
        if values[j].is_some()
            && byte(code, 0) == b'L'
            && freq(code).is_some()
            && reference.map_or(true, |r| columns.present[j] > columns.present[r])
        {
            reference = Some(j);
        }
    }

    codes
        .iter()
        .enumerate()
        .map(|(j, &code)| {
            let mut plan = Plan {
                predictor: Predictor::None,
                ssi_reference: None,
            };
            if values[j].is_none() {
                return plan;
            }
            plan.ssi_reference = columns.same_band(b'S', code, Some(j));
            let (Some(reference), Some(f_target)) = (reference, freq(code)) else {
                return plan;
            };
            let source = columns
                .same_band(b'L', code, None)
                .filter(|&s| freq(codes[s]).is_some())
                .unwrap_or(reference);
            let f_source = freq(codes[source]).unwrap();
            plan.predictor = match byte(code, 0) {
                b'L' if j != reference => {
                    let (num, den) = ratio(f_target, freq(codes[reference]).unwrap());
                    Predictor::Value {
                        reference,
                        num,
                        den,
                    }
                }
                b'C' | b'P' => {
                    let (num, den) = ratio(SPEED_OF_LIGHT, f_source);
                    Predictor::Value {
                        reference: source,
                        num,
                        den,
                    }
                }
                b'D' => {
                    let (num, den) = ratio(-f_target, f_source);
                    Predictor::Rate {
                        reference: source,
                        num,
                        den,
                        stencil: 1,
                        max_span: DOPPLER_MAX_SPAN * units_per_second,
                    }
                }
                _ => Predictor::None,
            };
            plan
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn codes(names: &[&str]) -> Vec<ObsCode> {
        names.iter().map(|n| ObsCode::from_v3(n).unwrap()).collect()
    }

    #[test]
    fn plans_follow_the_encoder_guidance() {
        let codes = codes(&["C1C", "L1C", "D1C", "S1C", "C2W", "L2W", "S2W", "X5X"]);
        let full = [1i64, 2, 3];
        let sparse = [1i64, ABSENT, 3];
        let values = [
            Some(&full[..]),
            Some(&full[..]),
            Some(&full[..]),
            Some(&full[..]),
            Some(&full[..]),
            Some(&sparse[..]),
            None,
            Some(&full[..]),
        ];
        let sat = SatId::new(SatSystem::GPS, 3);
        let plans = plan_satellite(sat, &codes, &values, &[None; 100], 10_000_000);
        let value = |reference, num, den| Predictor::Value {
            reference,
            num,
            den,
        };
        // L1C is the reference phase: present at the most epochs.
        assert_eq!(plans[1].predictor, Predictor::None);
        let (num, den) = ratio(299_792_458, 1_575_420_000);
        assert_eq!(plans[0].predictor, value(1, num, den));
        assert_eq!(
            plans[2].predictor,
            Predictor::Rate {
                reference: 1,
                num: -1,
                den: 1,
                stencil: 1,
                max_span: 3_000_000_000,
            }
        );
        // An S code's SSIs may use another S code on its band, not itself.
        assert_eq!(plans[3].ssi_reference, None);
        assert_eq!(plans[0].ssi_reference, Some(3));
        // L2W from L1C; C2W from L2W.  S2W was never observed.
        assert_eq!(plans[5].predictor, value(1, 60, 77));
        let (num, den) = ratio(299_792_458, 1_227_600_000);
        assert_eq!(plans[4].predictor, value(5, num, den));
        assert_eq!(plans[4].ssi_reference, None);
        assert_eq!(plans[6].predictor, Predictor::None);
        // An unknown observation type gets no predictor.
        assert_eq!(plans[7].predictor, Predictor::None);
    }

    #[test]
    fn glonass_channels_from_header_or_data() {
        let header = format!(
            "{:<60}GLONASS SLOT / FRQ #\n",
            "  3 R01  1 R02 -4 R24  2 R25 1x"
        );
        let channels = glonass_channels(&header);
        assert_eq!(channels[1], Some(1));
        assert_eq!(channels[2], Some(-4));
        assert_eq!(channels[24], Some(2));
        assert_eq!(channels[25], Some(1));
        assert_eq!(channels[3], None);

        // Without a header record, the channel comes from code minus carrier.
        let f = frequency(SatSystem::GLONASS, b'1', -3).unwrap();
        let phase: Vec<i64> = (0..50).map(|i| 1_000_000_000 + i * 3_456_789).collect();
        let code: Vec<i64> = phase
            .iter()
            .map(|&l| {
                i64::try_from(rounded(i128::from(l) * 299_792_458, i128::from(f)).unwrap()).unwrap()
                    + 7
            })
            .collect();
        assert_eq!(estimate_glonass_channel(&code, &phase), -3);
    }
}
