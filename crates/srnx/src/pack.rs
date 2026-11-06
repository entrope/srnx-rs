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

use std::collections::VecDeque;

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

/// A [`Block`] packed into a `u64`: the kind in the top byte, then the
/// size index and the bit width or Rice parameter, or a run length.  The
/// search keeps its best edge in one, so choosing between edges needs no
/// branches.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct Choice(u64);

impl Choice {
    const SLEB: u64 = 1 << 56;
    const ZEROS: u64 = 2 << 56;
    const MATRIX: u64 = 3 << 56;
    const RICE: u64 = 4 << 56;

    fn block(self) -> Block {
        let low = self.0 & ((1 << 56) - 1);
        #[allow(clippy::cast_possible_truncation)]
        let (size, param) = ((low >> 8) as usize, (low & 0xFF) as usize);
        match self.0 & !((1 << 56) - 1) {
            Self::SLEB => Block::Sleb(low as usize),
            Self::ZEROS => Block::Zeros(low as usize),
            Self::MATRIX => Block::Matrix { size, bits: param },
            Self::RICE => Block::Rice { size, r: param },
            _ => Block::Unset,
        }
    }
}

impl From<Block> for Choice {
    fn from(block: Block) -> Self {
        Self(match block {
            Block::Unset => 0,
            Block::Sleb(n) => Self::SLEB | n as u64,
            Block::Zeros(n) => Self::ZEROS | n as u64,
            Block::Matrix { size, bits } => Self::MATRIX | (size as u64) << 8 | bits as u64,
            Block::Rice { size, r } => Self::RICE | (size as u64) << 8 | r as u64,
        })
    }
}

/// Make `(cost, choice)` the best edge if it is cheaper, without a branch.
#[inline]
fn offer(best: &mut (u64, Choice), cost: u64, choice: Block) {
    let choice = Choice::from(choice);
    let better = cost < best.0;
    best.0 = if better { cost } else { best.0 };
    best.1 = Choice(if better { choice.0 } else { best.1 .0 });
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
    /// Sum of the bit widths of the nonzero values of each row in `widths`.
    width_sums: [u64; 8],
    /// Sum of the SLEB128 lengths of each row in `widths`.
    sleb_sums: [u64; 8],
    /// Runs of consecutive coded values, as (first, count) in `coded`.
    runs: Vec<(usize, usize)>,
    /// Block choices of the best delta order so far, and of a candidate,
    /// one per coded value.
    best: Vec<Choice>,
    cand: Vec<Choice>,
    search: Search,
}

/// The block search's working arrays.
#[derive(Default)]
struct Search {
    /// `cost[i]`: fewest bytes coding the values from position `i`.
    cost: Vec<u64>,
    /// Prefix sums of SLEB128 lengths.
    sleb_sum: Vec<u64>,
    /// Candidate ends `j` of an SLEB128 run, with `sleb_sum[j] + cost[j]`,
    /// for a sliding-window minimum: ends decrease and values increase
    /// from front to back.
    sleb_ends: VecDeque<(usize, u64)>,
    /// Row `k` holds the widest value of each window of `MATRIX_SIZES[k]`.
    window: Vec<u8>,
    /// Prefix sums of capped `ZigZag` values.
    u_sum: Vec<u64>,
    /// Position-major prefix sums of capped Rice quotients plus one, one
    /// per Rice parameter.  Sums wrap; differences over a block are exact.
    q_sum: Vec<u32>,
}

/// The split-Rice blocks from one position of [`Search::rice`].
struct RiceStep<'a> {
    search: &'a Search,
    i: usize,
    remaining: usize,
    /// Rice parameters with prefix sums, and the largest one.
    width: usize,
    r_max: usize,
}

impl RiceStep<'_> {
    /// Offer the best split-Rice block of size `RICE_SIZES[SIZE]`, starting
    /// the parameter search from `r_start`.  Returns false if the block
    /// does not fit, nor any larger one.
    #[inline]
    fn size<const SIZE: usize>(
        &self,
        r_start: &mut Option<usize>,
        best: &mut (u64, Choice),
    ) -> bool {
        let m = RICE_SIZES[SIZE];
        if self.remaining < m {
            return false;
        }
        let (search, i, width) = (self.search, self.i, self.width);
        let mean = search.u_sum[i + m].wrapping_sub(search.u_sum[i]) / m as u64;
        let r_hi = bit_len(mean).min(self.r_max);
        let qa = &search.q_sum[i * width..(i + 1) * width];
        let qb = &search.q_sum[(i + m) * width..(i + m + 1) * width];
        let bits = |r: usize| (r * m) as u64 + u64::from(qb[r].wrapping_sub(qa[r]));
        // A block's bit count is convex in r, so walk downhill.
        let start = r_start.map_or(r_hi, |s| s.min(r_hi));
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
        *r_start = Some(r);
        let cost = 1 + best_bits.div_ceil(8) + search.cost[i + m];
        offer(best, cost, Block::Rice { size: SIZE, r });
        true
    }
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
    fn matrix(&mut self, vals: &[i64], widths: &[u8], choices: &mut [Choice]) -> u64 {
        let n = vals.len();
        self.fill_windows(widths);
        self.sleb_sum.clear();
        self.sleb_sum.push(0);
        let mut sum = 0;
        for &w in widths {
            sum += u64::from(w).div_ceil(7);
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
                offer(best, cost, Block::Matrix { size, bits });
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
    fn rice(&mut self, vals: &[i64], choices: &mut [Choice]) -> u64 {
        let n = vals.len();
        let u_max = vals.iter().map(|&v| zigzag(v)).max().unwrap_or(0);
        let r_max = bit_len(u_max.min(RICE_U_CAP)).min(MAX_RICE_R);
        let width = r_max + 1;
        self.u_sum.clear();
        self.u_sum.push(0);
        self.sleb_sum.clear();
        self.sleb_sum.push(0);
        self.q_sum.clear();
        self.q_sum.resize((n + 1) * width, 0);
        let (mut u_sum, mut sleb_sum) = (0u64, 0u64);
        for (i, &v) in vals.iter().enumerate() {
            let u = zigzag(v);
            u_sum = u_sum.wrapping_add(u.min(RICE_U_CAP));
            self.u_sum.push(u_sum);
            sleb_sum += sleb128_len(v) as u64;
            self.sleb_sum.push(sleb_sum);
            let (prev, next) = self.q_sum[i * width..(i + 2) * width].split_at_mut(width);
            for (r, (next, &prev)) in next.iter_mut().zip(&*prev).enumerate() {
                #[allow(clippy::cast_possible_truncation)]
                let q = (u >> r).min(RICE_Q_CAP) as u32;
                *next = prev.wrapping_add(q + 1);
            }
        }
        self.walk(vals, choices, |search, i, best| {
            let rice = RiceStep {
                search,
                i,
                remaining: n - i,
                width,
                r_max,
            };
            // The best parameter for one size is a good start for the next.
            // Each size is its own instance, so the mean divides by a
            // constant.
            let mut r = None;
            let _ = rice.size::<0>(&mut r, best)
                && rice.size::<1>(&mut r, best)
                && rice.size::<2>(&mut r, best)
                && rice.size::<3>(&mut r, best)
                && rice.size::<4>(&mut r, best)
                && rice.size::<5>(&mut r, best)
                && rice.size::<6>(&mut r, best);
        })
    }

    /// The shortest-path walk shared by both codings: zero runs and SLEB128
    /// runs here, coding-specific blocks from `blocks(self, i, best)`.
    /// `sleb_sum` must be filled.  Writes one choice per value.
    fn walk(
        &mut self,
        vals: &[i64],
        choices: &mut [Choice],
        mut blocks: impl FnMut(&Self, usize, &mut (u64, Choice)),
    ) -> u64 {
        let n = vals.len();
        assert_eq!(choices.len(), n);
        self.cost.clear();
        self.cost.resize(n + 1, 0);
        let mut zeros = 0;
        self.sleb_ends.clear();
        for i in (0..n).rev() {
            let mut best = (u64::MAX, Choice::default());
            zeros = if vals[i] == 0 { zeros + 1 } else { 0 };
            // Only the longest zero run: splitting one never helps.
            if zeros > 0 {
                let cost = 1 + uleb128_len(zeros as u64 - 1) as u64 + self.cost[i + zeros];
                best = (cost, Choice::from(Block::Zeros(zeros)));
            }
            blocks(self, i, &mut best);
            // An SLEB128 run to `j` costs `2 + sleb_sum[j] - sleb_sum[i] +
            // cost[j]`, so the best is the least `sleb_sum[j] + cost[j]`
            // over the window `i + 1..=i + MAX_SLEB_RUN`: a sliding-window
            // minimum.  On ties the shortest run wins.
            let j = i + 1;
            let g = self.sleb_sum[j] + self.cost[j];
            while self.sleb_ends.back().is_some_and(|&(_, h)| h >= g) {
                self.sleb_ends.pop_back();
            }
            self.sleb_ends.push_back((j, g));
            if self.sleb_ends[0].0 > i + MAX_SLEB_RUN {
                self.sleb_ends.pop_front();
            }
            let (j, g) = self.sleb_ends[0];
            offer(&mut best, 2 + g - self.sleb_sum[i], Block::Sleb(j - i));
            self.cost[i] = best.0;
            choices[i] = best.1;
        }
        self.cost[0]
    }
}

impl Scratch {
    /// Bit widths of delta row `order`, the sum of those of nonzero values,
    /// and the sum of their SLEB128 lengths, computed on first use.  Rows
    /// are computed in order.
    fn widths(&mut self, order: usize, n: usize) {
        let row = order * n..(order + 1) * n;
        if self.widths.len() < row.end {
            self.widths.resize(row.end, 0);
            let (mut bits, mut bytes) = (0, 0);
            for (w, &d) in self.widths[row.clone()].iter_mut().zip(&self.deltas[row]) {
                *w = width(d);
                bits += if d == 0 { 0 } else { u64::from(*w) };
                bytes += u64::from(*w).div_ceil(7);
            }
            (self.width_sums[order], self.sleb_sums[order]) = (bits, bytes);
        }
    }

    /// Choose blocks for each run of delta row `order` into `cand`, and
    /// return their size in bytes plus a byte of initial state per level.
    fn search_order(&mut self, coding: Coding, order: usize, n: usize) -> u64 {
        self.cand.clear();
        self.cand.resize(n, Choice::default());
        let row = order * n;
        let mut total = order as u64;
        for &(start, len) in &self.runs {
            let range = row + start..row + start + len;
            let choices = &mut self.cand[start..start + len];
            total += match coding {
                Coding::Matrix => {
                    self.search
                        .matrix(&self.deltas[range.clone()], &self.widths[range], choices)
                }
                Coding::Rice => self.search.rice(&self.deltas[range], choices),
            };
        }
        total
    }

    /// Choose the bit-matrix delta order for the `n` coded values, leaving
    /// its block choices in `best`: orders 0, 1, … until two fail to
    /// improve.  Bounds avoid searching orders whose outcome is known, so
    /// the choice is the same as searching them all.  Every nonzero value
    /// costs at least its bit width in any block, which bounds an order's
    /// cost from below, and coding every value in SLEB128 runs bounds it
    /// from above.
    fn matrix_order(&mut self, n: usize) -> usize {
        self.widths.clear();
        let orders = 8.min(n);
        let runs: u64 = self
            .runs
            .iter()
            .map(|&(_, len)| 2 * len.div_ceil(MAX_SLEB_RUN) as u64)
            .sum();
        let lower = |s: &Self, order: usize| order as u64 + s.width_sums[order].div_ceil(8);
        let upper = |s: &Self, order: usize| order as u64 + s.sleb_sums[order] + runs;
        let (mut best_cost, mut best_order, mut misses) = (u64::MAX, 0, 0);
        for order in 0..orders {
            self.widths(order, n);
            // If the next order surely beats this one and the best so far,
            // this one cannot affect the outcome, unless it would be the
            // second miss in a row that ends the scan.
            if order + 1 < orders && misses == 0 {
                self.widths(order + 1, n);
                let next = upper(self, order + 1);
                if next < lower(self, order) && next < best_cost {
                    continue;
                }
            }
            let cost = if lower(self, order) >= best_cost {
                u64::MAX
            } else {
                self.search_order(Coding::Matrix, order, n)
            };
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
        let order_bits = order_bits(&self.coded);
        let estimates = order_bits.map(|(bits, _)| bits);
        let guess = best_estimate(&estimates);
        let mut best_total = self.search_order(Coding::Rice, guess, n);
        std::mem::swap(&mut self.best, &mut self.cand);
        let mut best_order = guess;
        for step in [-1isize, 1] {
            let mut order = guess;
            while let Some(next) = order.checked_add_signed(step).filter(|&o| o < 8.min(n)) {
                order = next;
                #[allow(clippy::cast_precision_loss)]
                if estimates[order] as f64 > estimates[guess] as f64 * (1.0 + ORDER_MARGIN)
                    || order_bits[order].1 >= best_total
                {
                    break;
                }
                let total = self.search_order(Coding::Rice, order, n);
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

/// Estimated bits of each delta order's residuals for `values`: their
/// `ZigZag` bit lengths plus a byte of initial state per level.  Orders
/// `values.len()` and above do not apply.  The deltas are computed as they
/// are summed, in one pass.
fn order_estimates(values: &[i64]) -> [u64; 8] {
    order_bits(values).map(|(bits, _)| bits)
}

/// For each delta order of `values`, as for [`order_estimates`], and a
/// lower bound in bytes on its split-Rice coding: a nonzero value costs at
/// least one bit more than its `ZigZag` bit length in any block.
fn order_bits(values: &[i64]) -> [(u64, u64); 8] {
    let mut previous = [0i64; 8];
    let mut bits = [0u64; 8];
    let mut nonzero = [0u64; 8];
    for &v in values {
        let mut d = v;
        for ((bits, nonzero), previous) in bits.iter_mut().zip(&mut nonzero).zip(&mut previous) {
            *bits += u64::from(64 - zigzag(d).leading_zeros());
            *nonzero += u64::from(d != 0);
            (d, *previous) = (d.wrapping_sub(*previous), d);
        }
    }
    let mut out = [(u64::MAX, u64::MAX); 8];
    for (order, out) in out.iter_mut().enumerate().take(values.len()) {
        let order_bytes = order as u64;
        *out = (
            8 * order_bytes + bits[order],
            order_bytes + (bits[order] + nonzero[order]).div_ceil(8),
        );
    }
    out
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

/// Gather the coded values of `entries` into `coded`, and their runs of
/// consecutive entries into `runs`.
fn gather(entries: &[Entry], coded: &mut Vec<i64>, runs: &mut Vec<(usize, usize)>) {
    coded.clear();
    runs.clear();
    let mut previous_coded = false;
    for e in entries {
        if let Entry::Coded(v) = *e {
            if previous_coded {
                runs.last_mut().unwrap().1 += 1;
            } else {
                runs.push((coded.len(), 1));
            }
            coded.push(v);
        }
        previous_coded = matches!(e, Entry::Coded(_));
    }
}

/// A cheap estimate of the split-Rice packed size of `entries`, in bytes,
/// to skip coding a `SOC2` candidate that is clearly larger than another.
pub(crate) fn estimate(entries: &[Entry], scratch: &mut Scratch) -> f64 {
    gather(entries, &mut scratch.coded, &mut scratch.runs);
    let n = scratch.coded.len();
    #[allow(clippy::cast_precision_loss)]
    let markers = marker_bytes(entries) as f64;
    if n == 0 {
        return markers;
    }
    let estimates = order_estimates(&scratch.coded);
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
    gather(entries, &mut scratch.coded, &mut scratch.runs);
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
                emit(&mut packed, vals, &scratch.best[coded..coded + run]);
                coded += run;
            }
        }
        i += run;
    }

    write_uleb128(out, packed.len() as u64);
    out.extend_from_slice(&packed);
}

/// Append the blocks `choices` picked for `vals`.
fn emit(out: &mut Vec<u8>, vals: &[i64], choices: &[Choice]) {
    let mut i = 0;
    while i < vals.len() {
        let block = choices[i].block();
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
    if bits == 0 {
        return;
    }
    let start = out.len();
    let row_bytes = m / 8;
    out.resize(start + bits * row_bytes, 0);
    let rows = &mut out[start..];
    let mut group = [0u64; 8];
    let mut vals = vals;
    for col in 0..row_bytes {
        for v in &mut group {
            *v = vals.next().unwrap_or(0);
        }
        // Each row's byte holds bit `shift` of the eight values.
        for (row, byte) in rows[col..]
            .iter_mut()
            .step_by(row_bytes)
            .take(bits)
            .enumerate()
        {
            let shift = bits - 1 - row;
            let mut b = 0u64;
            for &v in &group {
                b = b << 1 | (v >> shift) & 1;
            }
            *byte = u8::try_from(b).unwrap();
        }
    }
}

/// Appends bits most significant first, padding the last byte with zeros.
struct BitWriter<'a> {
    out: &'a mut Vec<u8>,
    /// The `n` (fewer than 8) bits not yet written, in the low bits.
    acc: u64,
    n: u32,
}

impl<'a> BitWriter<'a> {
    fn new(out: &'a mut Vec<u8>) -> Self {
        Self { out, acc: 0, n: 0 }
    }

    /// The low `len` (at most 56) bits of `value`.
    fn put(&mut self, value: u64, len: u32) {
        self.acc = self.acc << len | value;
        self.n += len;
        while self.n >= 8 {
            self.n -= 8;
            #[allow(clippy::cast_possible_truncation)]
            self.out.push((self.acc >> self.n) as u8);
        }
        self.acc &= (1 << self.n) - 1;
    }

    /// `q` one bits, then a zero bit.
    fn unary(&mut self, mut q: u64) {
        while q >= 48 {
            self.put((1 << 48) - 1, 48);
            q -= 48;
        }
        #[allow(clippy::cast_possible_truncation)]
        self.put(((1 << q) - 1) << 1, q as u32 + 1);
    }

    fn finish(self) {
        if self.n > 0 {
            #[allow(clippy::cast_possible_truncation)]
            self.out.push((self.acc << (8 - self.n)) as u8);
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

    /// The bounds in `matrix_order` choose as a scan that searches every
    /// order would: orders 0, 1, … until two in a row fail to improve.
    #[test]
    fn matrix_order_bounds_do_not_change_the_choice() {
        let mut scratch = Scratch::default();
        let mut columns: Vec<Vec<Entry>> = Vec::new();
        for seed in 1..40u64 {
            let degree = seed % 5;
            let noise_bits = seed % 23;
            let column = noise(3000, seed)
                .enumerate()
                .map(|(i, v)| {
                    let t = i64::try_from(i).unwrap();
                    let trend =
                        (0..=degree).fold(0i64, |acc, _| acc.wrapping_mul(t).wrapping_add(977));
                    if (i as u64 + seed) % 211 < seed % 7 {
                        Entry::Absent
                    } else {
                        Entry::Coded(trend.wrapping_add(v >> (63 - noise_bits)))
                    }
                })
                .collect();
            columns.push(column);
        }
        for entries in &columns {
            gather(entries, &mut scratch.coded, &mut scratch.runs);
            let n = scratch.coded.len();
            fill_deltas(&scratch.coded, &mut scratch.deltas);
            let chosen = scratch.matrix_order(n);
            let (mut best_cost, mut best_order, mut misses) = (u64::MAX, 0, 0);
            for order in 0..8.min(n) {
                scratch.widths(order, n);
                let cost = scratch.search_order(Coding::Matrix, order, n);
                if cost < best_cost {
                    (best_cost, best_order, misses) = (cost, order, 0);
                } else {
                    misses += 1;
                    if misses == 2 {
                        break;
                    }
                }
            }
            assert_eq!(chosen, best_order);
        }
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
                    let mut choices = vec![Choice::default(); vals.len()];
                    choices[0] = block.into();
                    choices[m] = Block::Zeros(2).into();
                    choices[m + 2] = Block::Sleb(2).into();
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
