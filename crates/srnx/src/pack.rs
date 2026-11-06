//! Packed observation data for the SRNX writer: the delta order and the
//! blocks of a `SOCD` or `SOC2` chunk's values.
//!
//! As the spec's encoder guidance suggests, and as the reference encoder
//! (`rnx2srnx`) does, blocks are chosen by a shortest path over positions
//! `0..=n`, walked from the end: each edge is a block, a zero run or an
//! SLEB128 run, weighted by its size in bytes.  Bit-matrix widths come from
//! a table of window maxima and split-Rice sizes from per-parameter prefix
//! sums, so each delta order costs O(n).  `SOCD` delta orders are scanned
//! upward until two fail to improve; `SOC2` orders start from an estimate
//! and move to neighbouring orders while that helps.

use crate::leb128::{sleb128_len, uleb128_len, write_sleb128, write_uleb128, zigzag};

/// `SOCD` bit-matrix block sizes, by the top three bits of the header byte.
const MATRIX_SIZES: [usize; 7] = [8, 16, 24, 32, 40, 48, 112];
/// `SOC2` split-Rice block sizes, by the top three bits of the header byte.
const RICE_SIZES: [usize; 7] = [8, 24, 56, 144, 240, 640, 2048];
/// Widest bit-matrix value: five bits of the header byte hold `bits - 1`.
const MAX_MATRIX_BITS: usize = 32;
/// Largest Rice parameter: five bits of the header byte.
const MAX_RICE_R: usize = 31;
/// Longest SLEB128 run the search considers, so a run's header and count
/// take two bytes.
const MAX_SLEB_RUN: usize = 16;
/// Rice quotients above this count as this in block sizes.  A block with
/// one is estimated at 2^20 bits or more, larger than SLEB128 runs of its
/// values (at most 2048 × 10 + 128 × 2 bytes), so the search never chooses
/// it; and a block's sum of capped quotients fits in a `u32`.
const RICE_Q_CAP: u64 = 1 << 20;
/// `ZigZag` values above this count as this when estimating Rice
/// parameters.
const RICE_U_CAP: u64 = 1 << 40;
/// A `SOC2` delta order next to the estimated best is tried only if its
/// estimate is within this fraction of the best's.
const ORDER_MARGIN: f64 = 0.05;

/// Which block coding the packed data uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Coding {
    /// `SOCD` bit matrices.
    Matrix,
    /// `SOC2` split-Rice blocks.
    Rice,
}

/// One slot of a column: one per satellite epoch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Entry {
    Absent,
    /// A scaled observation in a `SOC2` unpredicted-values block.
    Unpredicted(i64),
    /// A scaled observation, or with a `SOC2` predictor its residual, coded
    /// through the delta coder.
    Coded(i64),
}

/// An edge of the block search: what codes the values from a position.
#[derive(Debug, Clone, Copy, Default)]
enum Block {
    #[default]
    Unset,
    /// An SLEB128-value block of this many values.
    Sleb(usize),
    /// A zero-value block of this many values.
    Zeros(usize),
    /// A bit matrix of `MATRIX_SIZES[size]` values, `bits` wide.
    Matrix { size: usize, bits: usize },
    /// A split-Rice block of `RICE_SIZES[size]` values, parameter `r`.
    Rice { size: usize, r: usize },
}

impl Block {
    fn len(self) -> usize {
        match self {
            Self::Unset => unreachable!("block search left a position unset"),
            Self::Sleb(n) | Self::Zeros(n) => n,
            Self::Matrix { size, .. } => MATRIX_SIZES[size],
            Self::Rice { size, .. } => RICE_SIZES[size],
        }
    }
}

/// Reusable buffers for packing, so a file's many columns need few
/// allocations.
#[derive(Default)]
pub(crate) struct Scratch {
    coded: Vec<i64>,
    /// Row `d` (`n` values from `d * n`) holds the order-`d` deltas.
    deltas: Vec<i64>,
    /// Bit widths of `deltas`, for the rows filled so far.
    widths: Vec<u8>,
    /// Block choices of the best delta order so far, and of a candidate.
    best: Vec<Block>,
    cand: Vec<Block>,
    search: Search,
}

/// The block search's working arrays.
#[derive(Default)]
struct Search {
    /// `cost[i]`: fewest bytes coding the values from position `i`.
    cost: Vec<u64>,
    /// Prefix sums of SLEB128 lengths.
    sleb_sum: Vec<u32>,
    /// Row `k` holds the widest value of each window of `MATRIX_SIZES[k]`.
    window: Vec<u8>,
    /// Prefix sums of capped `ZigZag` values.
    u_sum: Vec<u64>,
    /// Position-major prefix sums of capped Rice quotients plus one, one
    /// per Rice parameter.  Sums wrap; differences over a block are exact.
    q_sum: Vec<u32>,
}

/// Two's-complement bit width of `v`, at least 1.
#[allow(clippy::cast_possible_truncation)]
fn width(v: i64) -> u8 {
    let magnitude = if v < 0 { !v } else { v };
    (65 - magnitude.leading_zeros()) as u8
}

/// Bits needed for `u`, or 1 for zero.
fn bit_len(u: u64) -> usize {
    (64 - u.leading_zeros()).max(1) as usize
}

/// Fill `deltas` with rows of order-0 to -7 deltas of `values`, each
/// starting from a zero decoder state.  The decoder adds with wrapping, so
/// deltas wrap too.
fn fill_deltas(values: &[i64], deltas: &mut Vec<i64>) {
    let n = values.len();
    deltas.clear();
    deltas.resize(8 * n, 0);
    deltas[..n].copy_from_slice(values);
    for order in 1..8 {
        let (prev, cur) = deltas[(order - 1) * n..(order + 1) * n].split_at_mut(n);
        if n > 0 {
            cur[0] = prev[0];
        }
        for i in 1..n {
            cur[i] = prev[i].wrapping_sub(prev[i - 1]);
        }
    }
}

impl Search {
    /// Choose `SOCD` blocks for `vals`, whose bit widths are `widths`,
    /// into `choices`; returns their size in bytes.
    fn matrix(&mut self, vals: &[i64], widths: &[u8], choices: &mut Vec<Block>) -> u64 {
        let n = vals.len();
        self.fill_windows(widths);
        self.sleb_sum.clear();
        self.sleb_sum.push(0);
        let mut sum = 0u32;
        for &w in widths {
            sum = sum.wrapping_add(u32::from(w).div_ceil(7));
            self.sleb_sum.push(sum);
        }
        self.walk(vals, choices, |search, i, best| {
            let remaining = n - i;
            for (size, &m) in MATRIX_SIZES.iter().enumerate() {
                if remaining < m {
                    break;
                }
                let bits = usize::from(search.window[size * n + i]);
                if bits > MAX_MATRIX_BITS {
                    break;
                }
                let cost = 1 + (bits * m / 8) as u64 + search.cost[i + m];
                if cost < best.0 {
                    *best = (cost, Block::Matrix { size, bits });
                }
            }
        })
    }

    /// Fill `window` with the widest value in each window of each
    /// bit-matrix size.  Windows of 8 to 48 extend the next smaller one by
    /// 8; 112 is 48 + 48 + 16.
    fn fill_windows(&mut self, widths: &[u8]) {
        let n = widths.len();
        self.window.clear();
        self.window.resize(MATRIX_SIZES.len() * n, 0);
        if n < 8 {
            return;
        }
        for i in 0..=n - 8 {
            self.window[i] = *widths[i..i + 8].iter().max().unwrap();
        }
        for k in 1..6 {
            let m = MATRIX_SIZES[k];
            if n < m {
                return;
            }
            let (lower, upper) = self.window.split_at_mut(k * n);
            let (row0, prev) = (&lower[..n], &lower[(k - 1) * n..]);
            for i in 0..=n - m {
                upper[i] = row0[i].max(prev[i + 8]);
            }
        }
        if n < 112 {
            return;
        }
        let (lower, upper) = self.window.split_at_mut(6 * n);
        let (row16, row48) = (&lower[n..2 * n], &lower[5 * n..]);
        for i in 0..=n - 112 {
            upper[i] = row48[i].max(row48[i + 48]).max(row16[i + 96]);
        }
    }

    /// Choose `SOC2` blocks for `vals` into `choices`; returns their size
    /// in bytes.
    fn rice(&mut self, vals: &[i64], choices: &mut Vec<Block>) -> u64 {
        let n = vals.len();
        let u_max = vals.iter().map(|&v| zigzag(v)).max().unwrap_or(0);
        let r_max = bit_len(u_max.min(RICE_U_CAP)).min(MAX_RICE_R);
        let width = r_max + 1;
        self.u_sum.clear();
        self.u_sum.push(0);
        self.sleb_sum.clear();
        self.sleb_sum.push(0);
        self.q_sum.clear();
        self.q_sum.resize(width, 0);
        let (mut u_sum, mut sleb_sum) = (0u64, 0u32);
        for (i, &v) in vals.iter().enumerate() {
            let u = zigzag(v);
            u_sum = u_sum.wrapping_add(u.min(RICE_U_CAP));
            self.u_sum.push(u_sum);
            sleb_sum = sleb_sum.wrapping_add(u32::try_from(sleb128_len(v)).unwrap());
            self.sleb_sum.push(sleb_sum);
            for r in 0..width {
                let q = u32::try_from((u >> r).min(RICE_Q_CAP)).unwrap();
                let next = self.q_sum[i * width + r].wrapping_add(q + 1);
                self.q_sum.push(next);
            }
        }
        self.walk(vals, choices, |search, i, best| {
            let remaining = n - i;
            // The best parameter for one size is a good start for the next.
            let mut r_start = None;
            for (size, &m) in RICE_SIZES.iter().enumerate() {
                if remaining < m {
                    break;
                }
                let mean = search.u_sum[i + m].wrapping_sub(search.u_sum[i]) / m as u64;
                let r_hi = bit_len(mean).min(r_max);
                let (qa, qb) = (i * width, (i + m) * width);
                let bits = |r: usize| {
                    (r * m) as u64
                        + u64::from(search.q_sum[qb + r].wrapping_sub(search.q_sum[qa + r]))
                };
                // A block's bit count is convex in r, so walk downhill.
                let start = r_start.map_or(r_hi, |s: usize| s.min(r_hi));
                let (mut r, mut best_bits) = (start, bits(start));
                while r > 0 && bits(r - 1) < best_bits {
                    r -= 1;
                    best_bits = bits(r);
                }
                if r == start {
                    while r < r_hi && bits(r + 1) < best_bits {
                        r += 1;
                        best_bits = bits(r);
                    }
                }
                r_start = Some(r);
                let cost = 1 + best_bits.div_ceil(8) + search.cost[i + m];
                if cost < best.0 {
                    *best = (cost, Block::Rice { size, r });
                }
            }
        })
    }

    /// The shortest-path walk shared by both codings: zero runs and SLEB128
    /// runs here, coding-specific blocks from `blocks(self, i, best)`.
    /// `sleb_sum` must be filled.
    fn walk(
        &mut self,
        vals: &[i64],
        choices: &mut Vec<Block>,
        mut blocks: impl FnMut(&Self, usize, &mut (u64, Block)),
    ) -> u64 {
        let n = vals.len();
        choices.clear();
        choices.resize(n, Block::Unset);
        self.cost.clear();
        self.cost.resize(n + 1, 0);
        let mut zeros = 0;
        for i in (0..n).rev() {
            let mut best = (u64::MAX, Block::Unset);
            zeros = if vals[i] == 0 { zeros + 1 } else { 0 };
            // Only the longest zero run: splitting one never helps.
            if zeros > 0 {
                let cost = 1 + uleb128_len(zeros as u64 - 1) as u64 + self.cost[i + zeros];
                best = (cost, Block::Zeros(zeros));
            }
            blocks(self, i, &mut best);
            let base = self.sleb_sum[i];
            for len in 1..=MAX_SLEB_RUN.min(n - i) {
                let bytes = self.sleb_sum[i + len].wrapping_sub(base);
                let cost = 2 + u64::from(bytes) + self.cost[i + len];
                if cost < best.0 {
                    best = (cost, Block::Sleb(len));
                }
            }
            self.cost[i] = best.0;
            choices[i] = best.1;
        }
        self.cost[0]
    }
}

impl Scratch {
    /// Bit widths of delta row `order`, computed on first use.
    fn widths(&mut self, order: usize, n: usize) -> &[u8] {
        let row = order * n..(order + 1) * n;
        if self.widths.len() < row.end {
            self.widths.resize(row.end, 0);
            for (w, &d) in self.widths[row.clone()]
                .iter_mut()
                .zip(&self.deltas[row.clone()])
            {
                *w = width(d);
            }
        }
        &self.widths[row]
    }

    /// Choose the bit-matrix delta order for the `n` coded values, leaving
    /// its block choices in `best`.  Each order's initial state costs a
    /// byte per level.
    fn matrix_order(&mut self, n: usize) -> usize {
        self.widths.clear();
        self.widths(0, n);
        let mut best_cost =
            self.search
                .matrix(&self.deltas[..n], &self.widths[..n], &mut self.best);
        let (mut best_order, mut misses) = (0, 0);
        for order in 1..8.min(n) {
            self.widths(order, n);
            let row = order * n..(order + 1) * n;
            let cost = order as u64
                + self
                    .search
                    .matrix(&self.deltas[row.clone()], &self.widths[row], &mut self.cand);
            if cost < best_cost {
                std::mem::swap(&mut self.best, &mut self.cand);
                (best_cost, best_order, misses) = (cost, order, 0);
            } else {
                misses += 1;
                if misses == 2 {
                    break;
                }
            }
        }
        best_order
    }

    /// Choose the split-Rice delta order for the `n` coded values, leaving
    /// its block choices in `best`.
    fn rice_order(&mut self, n: usize) -> usize {
        let estimates = order_estimates(&self.deltas, n);
        let guess = best_estimate(&estimates);
        let row = |order: usize| order * n..(order + 1) * n;
        let mut best_total =
            guess as u64 + self.search.rice(&self.deltas[row(guess)], &mut self.best);
        let mut best_order = guess;
        for step in [-1isize, 1] {
            let mut order = guess;
            while let Some(next) = order.checked_add_signed(step).filter(|&o| o < 8.min(n)) {
                order = next;
                #[allow(clippy::cast_precision_loss)]
                if estimates[order] as f64 > estimates[guess] as f64 * (1.0 + ORDER_MARGIN) {
                    break;
                }
                let total =
                    order as u64 + self.search.rice(&self.deltas[row(order)], &mut self.cand);
                if total >= best_total {
                    break;
                }
                std::mem::swap(&mut self.best, &mut self.cand);
                (best_total, best_order) = (total, order);
            }
        }
        best_order
    }
}

/// Estimated bits of each delta order's residuals: their `ZigZag` bit
/// lengths plus a byte of initial state per level.  Orders `n` and above
/// do not apply.
fn order_estimates(deltas: &[i64], n: usize) -> [u64; 8] {
    let mut estimates = [u64::MAX; 8];
    for (order, estimate) in estimates.iter_mut().enumerate().take(8.min(n)) {
        let bits: u64 = deltas[order * n..(order + 1) * n]
            .iter()
            .map(|&d| u64::from(64 - zigzag(d).leading_zeros()))
            .sum();
        *estimate = 8 * order as u64 + bits;
    }
    estimates
}

fn best_estimate(estimates: &[u64; 8]) -> usize {
    (0..8).min_by_key(|&order| estimates[order]).unwrap()
}

/// Bytes of the absent-value and unpredicted-value blocks of `entries`.
fn marker_bytes(entries: &[Entry]) -> usize {
    let mut bytes = 0;
    let mut i = 0;
    while i < entries.len() {
        let run = run_len(entries, i);
        match entries[i] {
            Entry::Coded(_) => {}
            Entry::Absent => bytes += 1 + uleb128_len(run as u64 - 1),
            Entry::Unpredicted(_) => {
                bytes += 1 + uleb128_len(run as u64 - 1);
                bytes += entries[i..i + run]
                    .iter()
                    .map(|e| match *e {
                        Entry::Unpredicted(v) => sleb128_len(v),
                        _ => unreachable!(),
                    })
                    .sum::<usize>();
            }
        }
        i += run;
    }
    bytes
}

/// Length of the run of entries of the same kind as `entries[i]`.
fn run_len(entries: &[Entry], i: usize) -> usize {
    let kind = std::mem::discriminant(&entries[i]);
    entries[i..]
        .iter()
        .take_while(|e| std::mem::discriminant(*e) == kind)
        .count()
}

/// Gather the coded values of `entries` into `coded`.
fn gather(entries: &[Entry], coded: &mut Vec<i64>) {
    coded.clear();
    coded.extend(entries.iter().filter_map(|e| match *e {
        Entry::Coded(v) => Some(v),
        _ => None,
    }));
}

/// A cheap estimate of the split-Rice packed size of `entries`, in bytes,
/// to skip coding a `SOC2` candidate that is clearly larger than another.
pub(crate) fn estimate(entries: &[Entry], scratch: &mut Scratch) -> f64 {
    gather(entries, &mut scratch.coded);
    let n = scratch.coded.len();
    #[allow(clippy::cast_precision_loss)]
    let markers = marker_bytes(entries) as f64;
    if n == 0 {
        return markers;
    }
    fill_deltas(&scratch.coded, &mut scratch.deltas);
    let estimates = order_estimates(&scratch.deltas, n);
    #[allow(clippy::cast_precision_loss)]
    let bits = estimates[best_estimate(&estimates)] as f64;
    bits / 8.0 + markers
}

/// Append packed observation data for `entries` at explicit scale `scale`
/// (1 for none) to `out`, starting with its ULEB128 length.
pub(crate) fn encode(
    out: &mut Vec<u8>,
    entries: &[Entry],
    scale: u64,
    coding: Coding,
    scratch: &mut Scratch,
) {
    gather(entries, &mut scratch.coded);
    let n = scratch.coded.len();
    let mut order = 0;
    if n > 0 {
        fill_deltas(&scratch.coded, &mut scratch.deltas);
        order = match coding {
            Coding::Matrix => scratch.matrix_order(n),
            Coding::Rice => scratch.rice_order(n),
        };
    }

    let mut packed = Vec::new();
    let explicit = if scale > 1 { 8 } else { 0 };
    write_uleb128(&mut packed, order as u64 + explicit);
    if scale > 1 {
        write_uleb128(&mut packed, scale);
    }
    packed.resize(packed.len() + order, 0); // initial state: zeros

    let row = order * n;
    let (mut i, mut coded) = (0, 0);
    while i < entries.len() {
        let run = run_len(entries, i);
        match entries[i] {
            Entry::Absent => {
                packed.push(0xFD);
                write_uleb128(&mut packed, run as u64 - 1);
            }
            Entry::Unpredicted(_) => {
                packed.push(0xFC);
                write_uleb128(&mut packed, run as u64 - 1);
                for e in &entries[i..i + run] {
                    if let Entry::Unpredicted(v) = *e {
                        write_sleb128(&mut packed, v);
                    }
                }
            }
            Entry::Coded(_) => {
                let vals = &scratch.deltas[row + coded..row + coded + run];
                // The order search chose blocks for all coded values; a
                // run that is only some of them needs its own.
                let choices = if run == n {
                    &scratch.best
                } else {
                    match coding {
                        Coding::Matrix => {
                            let widths = &scratch.widths[row + coded..row + coded + run];
                            scratch.search.matrix(vals, widths, &mut scratch.cand);
                        }
                        Coding::Rice => {
                            scratch.search.rice(vals, &mut scratch.cand);
                        }
                    }
                    &scratch.cand
                };
                emit(&mut packed, vals, choices);
                coded += run;
            }
        }
        i += run;
    }

    write_uleb128(out, packed.len() as u64);
    out.extend_from_slice(&packed);
}

/// Append the blocks `choices` picked for `vals`.
fn emit(out: &mut Vec<u8>, vals: &[i64], choices: &[Block]) {
    let mut i = 0;
    while i < vals.len() {
        let block = choices[i];
        let len = block.len();
        let vals = &vals[i..i + len];
        match block {
            Block::Unset => unreachable!(),
            Block::Sleb(_) => {
                out.push(0xFF);
                write_uleb128(out, len as u64 - 1);
                for &v in vals {
                    write_sleb128(out, v);
                }
            }
            Block::Zeros(_) => {
                out.push(0xFE);
                write_uleb128(out, len as u64 - 1);
            }
            Block::Matrix { size, bits } => {
                out.push(header(size, bits - 1));
                #[allow(clippy::cast_sign_loss)]
                write_matrix(out, vals.iter().map(|&v| v as u64), len, bits);
            }
            Block::Rice { size, r } => {
                out.push(header(size, r));
                let mask = (1u64 << r) - 1;
                write_matrix(out, vals.iter().map(|&v| zigzag(v) & mask), len, r);
                let mut bits = BitWriter::new(out);
                for &v in vals {
                    bits.unary(zigzag(v) >> r);
                }
                bits.finish();
            }
        }
        i += len;
    }
}

fn header(size: usize, low: usize) -> u8 {
    u8::try_from(size << 5 | low).unwrap()
}

/// Append a bit matrix of `m` values `bits` wide: their bits from most
/// significant, transposed, each row `m / 8` bytes.
fn write_matrix(out: &mut Vec<u8>, vals: impl Iterator<Item = u64>, m: usize, bits: usize) {
    let start = out.len();
    out.resize(start + bits * m / 8, 0);
    let rows = &mut out[start..];
    for (i, v) in vals.enumerate() {
        let (col, bit) = (i / 8, 7 - i % 8);
        for row in 0..bits {
            let b = (v >> (bits - 1 - row)) & 1;
            rows[row * m / 8 + col] |= u8::try_from(b).unwrap() << bit;
        }
    }
}

/// Appends bits most significant first, padding the last byte with zeros.
struct BitWriter<'a> {
    out: &'a mut Vec<u8>,
    acc: u8,
    n: u32,
}

impl<'a> BitWriter<'a> {
    fn new(out: &'a mut Vec<u8>) -> Self {
        Self { out, acc: 0, n: 0 }
    }

    fn bit(&mut self, b: bool) {
        self.acc = self.acc << 1 | u8::from(b);
        self.n += 1;
        if self.n == 8 {
            self.out.push(self.acc);
            (self.acc, self.n) = (0, 0);
        }
    }

    /// `q` one bits, then a zero bit.
    fn unary(&mut self, mut q: u64) {
        while q > 0 && self.n != 0 {
            self.bit(true);
            q -= 1;
        }
        while q >= 8 {
            self.out.push(0xFF);
            q -= 8;
        }
        for _ in 0..q {
            self.bit(true);
        }
        self.bit(false);
    }

    fn finish(self) {
        if self.n > 0 {
            self.out.push(self.acc << (8 - self.n));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{soc2, socd};

    /// A `SOCD` payload around packed data, as `socd::parse` expects.
    fn socd_payload(count: usize, packed: &[u8]) -> Vec<u8> {
        let mut p = b"G01\0L1C\0".to_vec();
        write_uleb128(&mut p, count as u64 - 1);
        p.extend_from_slice(&[0, 0]); // empty LLI and SSI
        p.extend_from_slice(packed);
        p
    }

    /// A predictor-free `SOC2` payload around packed data.
    fn soc2_payload(count: usize, packed: &[u8]) -> Vec<u8> {
        let mut p = b"G01\0L1C\0".to_vec();
        write_uleb128(&mut p, count as u64 - 1);
        p.extend_from_slice(&[0, 0, 0]); // LLI, mode, SSI
        p.extend_from_slice(packed);
        p
    }

    fn round_trip(values: &[Option<i64>], scale: u64) {
        let entries: Vec<Entry> = values
            .iter()
            .map(|v| v.map_or(Entry::Absent, Entry::Coded))
            .collect();
        let want: Vec<i64> = values.iter().map(|v| v.unwrap_or(socd::ABSENT)).collect();
        let mut scratch = Scratch::default();
        for coding in [Coding::Matrix, Coding::Rice] {
            let mut packed = Vec::new();
            encode(&mut packed, &entries, scale, coding, &mut scratch);
            let got = match coding {
                Coding::Matrix => socd::parse(&socd_payload(values.len(), &packed)).unwrap(),
                Coding::Rice => {
                    soc2::parse(&soc2_payload(values.len(), &packed))
                        .unwrap()
                        .data
                }
            };
            assert_eq!(got.values, want, "{coding:?}");
            assert_eq!(u64::from(got.scale), scale, "{coding:?}");
        }
    }

    /// A deterministic pseudo-random sequence (xorshift).
    fn noise(n: usize, seed: u64) -> impl Iterator<Item = i64> {
        let mut x = seed;
        (0..n).map(move |_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            #[allow(clippy::cast_possible_wrap)]
            let v = x as i64;
            v
        })
    }

    #[test]
    fn smooth_noisy_sparse_and_extreme_columns_round_trip() {
        // A carrier-phase-like ramp: high delta orders win.
        let ramp: Vec<_> = (0..3000i64)
            .map(|i| Some(20_000_000_000 + i * 1_234_567 + i * i * 37))
            .collect();
        round_trip(&ramp, 1);
        // Small noise around a level, with gaps and zeros.
        let noisy: Vec<_> = noise(2500, 1)
            .enumerate()
            .map(|(i, v)| match i % 97 {
                0..=4 => None,
                5..=40 => Some(0),
                _ => Some(v % 5000),
            })
            .collect();
        round_trip(&noisy, 1);
        // Full-range values and the extremes: wrapping deltas, bit widths
        // over 32, quotients over the cap.
        let wild: Vec<_> = noise(700, 7)
            .map(Some)
            .chain([Some(i64::MAX), Some(i64::MIN + 1), Some(0), Some(-1)])
            .collect();
        round_trip(&wild, 1);
        // Mostly small values with a few huge outliers.
        let outliers: Vec<_> = noise(4200, 3)
            .enumerate()
            .map(|(i, v)| Some(if i % 501 == 0 { v } else { v % 9 }))
            .collect();
        round_trip(&outliers, 1);
    }

    #[test]
    fn short_columns_round_trip() {
        for n in 1..20i64 {
            let column: Vec<_> = (0..n).map(|i| Some(i * i - 7)).collect();
            round_trip(&column, 1);
        }
        round_trip(&[None, None], 1);
        round_trip(&[Some(5)], 20);
    }

    /// Every bit-matrix and split-Rice size and width decodes, alone and
    /// between other blocks.
    #[test]
    fn every_block_kind_decodes() {
        for size in 0..7 {
            for (bits, r) in [(1, 0), (2, 1), (13, 9), (32, 31)] {
                for coding in [Coding::Matrix, Coding::Rice] {
                    let (block, m) = match coding {
                        Coding::Matrix => (Block::Matrix { size, bits }, MATRIX_SIZES[size]),
                        Coding::Rice => (Block::Rice { size, r }, RICE_SIZES[size]),
                    };
                    let limit = 1i64 << (bits - 1);
                    let mut vals: Vec<i64> = noise(m, 5)
                        .map(|v| match coding {
                            Coding::Matrix => v.rem_euclid(2 * limit) - limit,
                            Coding::Rice => v % 1000,
                        })
                        .collect();
                    vals.extend([0, 0, 7, -3]);
                    let mut choices = vec![Block::Unset; vals.len()];
                    choices[0] = block;
                    choices[m] = Block::Zeros(2);
                    choices[m + 2] = Block::Sleb(2);
                    let mut packed = vec![0]; // schema: order 0, no scale
                    emit(&mut packed, &vals, &choices);
                    let mut framed = Vec::new();
                    write_uleb128(&mut framed, packed.len() as u64);
                    framed.extend_from_slice(&packed);
                    let got = match coding {
                        Coding::Matrix => socd::parse(&socd_payload(vals.len(), &framed)).unwrap(),
                        Coding::Rice => {
                            soc2::parse(&soc2_payload(vals.len(), &framed))
                                .unwrap()
                                .data
                        }
                    };
                    assert_eq!(got.values, vals, "{coding:?} size {size}");
                }
            }
        }
    }

    #[test]
    fn unpredicted_values_need_a_predictor_to_decode() {
        // The packed data itself is coding-agnostic for markers.
        let entries = [
            Entry::Coded(3),
            Entry::Unpredicted(-40),
            Entry::Unpredicted(41),
            Entry::Absent,
            Entry::Coded(4),
        ];
        let mut packed = Vec::new();
        encode(
            &mut packed,
            &entries,
            1,
            Coding::Rice,
            &mut Scratch::default(),
        );
        assert_eq!(marker_bytes(&entries), 2 + 1 + 1 + 2);
        let mut p = b"G01\0L1C\0".to_vec();
        p.extend_from_slice(&[4, 0, 1, 0, 0, 2, 1]); // count, LLI, mode 1, SSI, ref 0, num 1, den 1
        p.extend_from_slice(&packed);
        let soc2 = soc2::parse(&p).unwrap();
        assert_eq!(soc2.data.values, [3, -40, 41, socd::ABSENT, 4]);
        assert_eq!(soc2.unpredicted, [false, true, true, false, false]);
    }
}
