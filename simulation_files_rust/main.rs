use rand::prelude::*;
use rand_distr::{Exp, Distribution};
use csv::WriterBuilder;
use serde::Serialize;
use std::fs::OpenOptions;
use std::io::BufWriter;
use std::path::Path;
use std::time::Instant;

// =============================================================================
// FenwickTree
//
// A Binary Indexed Tree (BIT) supporting:
//   update(i, delta)  — add delta to position i         — O(log n)
//   query(i)          — prefix sum [0..=i]               — O(log n)
//   find_kth(k)       — smallest i with prefix_sum > k   — O(log n)
//
// All indices are 0-based externally; converted to 1-based internally.
//
// `find_kth` is the key operation for weighted sampling: draw a uniform
// integer r in [0, total), then find_kth(r) returns the sampled index.
// =============================================================================
struct FenwickTree {
    n: usize,
    tree: Vec<i64>,
}

impl FenwickTree {
    fn new(size: usize) -> Self {
        FenwickTree {
            n: size,
            tree: vec![0; size + 1],
        }
    }

    /// Add `delta` to position `i` (0-indexed).
    fn update(&mut self, i: usize, delta: i64) {
        let mut i = i + 1; // convert to 1-indexed
        while i <= self.n {
            self.tree[i] += delta;
            i += i & i.wrapping_neg();
        }
    }

    /// Return prefix sum [0..=i] (0-indexed).
    fn query(&self, i: usize) -> i64 {
        let mut i = i + 1; // convert to 1-indexed
        let mut s = 0i64;
        while i > 0 {
            s += self.tree[i];
            i -= i & i.wrapping_neg();
        }
        s
    }

    /// Total sum of all elements.
    fn total(&self) -> i64 {
        if self.n == 0 { 0 } else { self.query(self.n - 1) }
    }

    /// Return the smallest 0-indexed position whose prefix sum exceeds `k`.
    /// Used for O(log n) weighted sampling.
    fn find_kth(&self, mut k: i64) -> usize {
        let mut pos = 0usize;
        // Find highest power of 2 that is <= n
        let mut step = 1usize << (usize::BITS - self.n.leading_zeros()) as usize >> 1;
        if step == 0 { step = 1; }
        while step > 0 {
            let next = pos + step;
            if next <= self.n && self.tree[next] <= k {
                k -= self.tree[next];
                pos = next;
            }
            step >>= 1;
        }
        pos // 0-indexed result
    }
}

// =============================================================================
// LevelStore
//
// Maintains the set of occupied levels and their miner counts.
//
// Internal representation:
//   levels      — sorted Vec<i64> of distinct occupied level values
//   counts      — parallel Vec<i64>: counts[i] = miners at levels[i]
//   sum_sq_val  — running sum of counts[i]^2, updated incrementally — O(1)
//   bit         — FenwickTree over counts[], rebuilt on structural changes
//                 (level insertions / deletions). Used for O(log n) sampling.
//
// Why a sorted Vec rather than a BTreeMap?
//   The number of distinct levels is tiny (empirically O(log N) for this
//   model), so cache-friendly Vec operations beat pointer-chasing trees.
//   Binary search (partition_point) gives O(log(distinct)) for lookups.
//
// The FenwickTree is sized to len(levels) and fully rebuilt on every
// structural change (insertion or deletion of a level). Because structural
// changes happen at most once per step and distinct levels stays tiny,
// rebuilds are cheap — and this design avoids the index-mismatch bug that
// arises when trying to maintain a fixed-capacity BIT alongside a
// shrinking/growing Vec.
// =============================================================================
struct LevelStore {
    levels:      Vec<i64>,    // sorted, distinct level values
    counts:      Vec<i64>,    // counts[i] = miners at levels[i]
    sum_sq_val:  i64,         // sum of counts[i]^2, maintained incrementally
    bit:         FenwickTree, // BIT over counts[], in sync with levels/counts
}

impl LevelStore {
    fn new(initial_level: i64, initial_count: i64) -> Self {
        let mut bit = FenwickTree::new(1);
        bit.update(0, initial_count);
        LevelStore {
            levels:     vec![initial_level],
            counts:     vec![initial_count],
            sum_sq_val: initial_count * initial_count,
            bit,
        }
    }

    // ── internal ──────────────────────────────────────────────────────────────

    /// Rebuild the BIT from scratch to match current counts[]. O(distinct log distinct).
    fn rebuild_bit(&mut self) {
        let n = self.levels.len();
        self.bit = FenwickTree::new(n);
        for (i, &c) in self.counts.iter().enumerate() {
            self.bit.update(i, c);
        }
    }

    // ── public API ────────────────────────────────────────────────────────────

    fn get_level(&self, idx: usize) -> i64 { self.levels[idx] }

    /// Sum of counts[i]^2 — needed for the interaction rate Ri. O(1).
    fn sum_sq(&self) -> i64 { self.sum_sq_val }

    /// Add `delta` miners to `level_val`.
    /// Creates the level if new; removes it if count drops to zero.
    fn add_miners(&mut self, level_val: i64, delta: i64) {
        let ins = self.levels.partition_point(|&x| x < level_val);

        if ins < self.levels.len() && self.levels[ins] == level_val {
            // ── existing level ────────────────────────────────────────────────
            let old_c = self.counts[ins];
            let new_c = old_c + delta;
            self.sum_sq_val += new_c * new_c - old_c * old_c;
            self.counts[ins] = new_c;

            if new_c == 0 {
                // Level vanishes — structural change, must rebuild BIT
                self.levels.remove(ins);
                self.counts.remove(ins);
                self.rebuild_bit();
            } else {
                // Count changed only — cheap in-place BIT update
                self.bit.update(ins, delta);
            }
        } else {
            // ── new level ─────────────────────────────────────────────────────
            // Insertion shifts all BIT indices to the right of `ins`,
            // so we must rebuild the BIT from scratch.
            self.levels.insert(ins, level_val);
            self.counts.insert(ins, delta);
            self.sum_sq_val += delta * delta;
            self.rebuild_bit();
        }
    }

    // ── sampling ──────────────────────────────────────────────────────────────

    /// Sample a level index proportional to counts[i].
    /// O(log(distinct levels)) via BIT find_kth.
    /// The BIT is always sized to match the current levels Vec,
    /// so the returned index is guaranteed to be in range.
    fn sample_idx(&self, rng: &mut ThreadRng) -> usize {
        let total = self.bit.total();
        debug_assert!(total > 0);
        let k = rng.gen_range(0..total);
        self.bit.find_kth(k)
    }

    /// Sample j proportional to counts[j] * prefix_count[j-1].
    /// This is the weight for the "leading" level in an interaction pair.
    /// O(distinct levels) linear scan — fast since distinct levels is tiny.
    /// Returns None if all miners share the same level (no valid pair exists).
    fn sample_j(&self, rng: &mut ThreadRng) -> Option<usize> {
        let n = self.levels.len();
        let mut weights: Vec<i64> = Vec::with_capacity(n);
        let mut cum: i64 = 0;
        for i in 0..n {
            weights.push(self.counts[i] * cum);
            cum += self.counts[i];
        }
        let total: i64 = weights.iter().sum();
        if total == 0 {
            return None; // all miners at same level
        }
        let mut r = rng.gen_range(0..total);
        for (i, &w) in weights.iter().enumerate() {
            if r < w {
                return Some(i);
            }
            r -= w;
        }
        Some(n - 1) // numerical safety
    }

    /// Sample i < j proportional to counts[i].
    /// O(log(distinct levels)) via BIT prefix-sum query + find_kth.
    fn sample_i_below_j(&self, j: usize, rng: &mut ThreadRng) -> usize {
        debug_assert!(j > 0);
        let prefix_total = self.bit.query(j - 1);
        debug_assert!(prefix_total > 0);
        let k = rng.gen_range(0..prefix_total);
        self.bit.find_kth(k)
    }
}

// =============================================================================
// simulate_speed  —  O(N · T · log N)
//
// Same CTMC model as the original O(N²T) version; LevelStore + FenwickTree
// replace all O(N) list operations with O(log N) equivalents.
//
//   Operation                       Original        This version
//   ──────────────────────────────────────────────────────────────
//   Compute sum_sq for rate         O(distinct)     O(1)  incremental
//   Weighted sample — birth         O(N)            O(log N)
//   Weighted sample j — interact    O(N)            O(distinct) ≈ O(log N)
//   Weighted sample i < j           O(N)            O(log N)
//   Insert / remove a level         O(N) Vec shift  O(distinct) + BIT rebuild
//   Steps per unit time             O(N)            O(N)   (rate unchanged)
//   ──────────────────────────────────────────────────────────────
//   Total                           O(N²T)          O(N·T·log N)
// =============================================================================
fn simulate_speed(n: usize, lam: f64, omega: f64, t_final: f64) -> f64 {
    let mut rng = thread_rng();
    let n_f = n as f64;

    // All N miners start at level 0
    let mut store = LevelStore::new(0, n as i64);

    let mut t       = 0.0f64;
    let mut mean    = 0.0f64;
    let t0          = 0.5 * t_final;
    let mut mean_t0 = 0.0f64;
    let mut t0_recorded = false;

    while t < t_final {

        // ── rates ─────────────────────────────────────────────────────────────
        //
        //   Rb = λ · N                              (total birth rate)
        //   Ri = (ω / 2N) · (N² − Σᵢ nᵢ²)          (total interaction rate)
        //   R  = Rb + Ri                             (total event rate)
        //
        let sq  = store.sum_sq();                    // O(1)
        let rb  = lam * n_f;
        let ri  = (omega / (2.0 * n_f)) * (n_f * n_f - sq as f64);
        let r   = rb + ri;

        // ── exponential holding time ──────────────────────────────────────────
        t += Exp::new(r).expect("rate must be positive").sample(&mut rng);

        // ── choose event ──────────────────────────────────────────────────────
        if rng.gen::<f64>() < rb / r {

            // ── BIRTH ─────────────────────────────────────────────────────────
            //
            // Pick a miner uniformly (proportional to counts).
            // That miner advances from level k to k+1.
            //
            let i = store.sample_idx(&mut rng);      // O(log distinct)
            let k = store.get_level(i);

            store.add_miners(k,     -1);
            store.add_miners(k + 1,  1);

            mean += 1.0 / n_f;

        } else {

            // ── INTERACTION ───────────────────────────────────────────────────
            //
            // Pick an ordered pair of levels (k < ℓ) where:
            //   j = index of ℓ,  sampled ∝ counts[j] · (miners below j)
            //   i = index of k,  sampled ∝ counts[i]  for i < j
            //
            // One miner at level k jumps to level ℓ.
            //
            let Some(j) = store.sample_j(&mut rng) else {
                continue; // all miners at same level — skip
            };
            let ell = store.get_level(j);

            let i = store.sample_i_below_j(j, &mut rng); // O(log distinct)
            let k = store.get_level(i);

            store.add_miners(k,   -1);
            store.add_miners(ell,  1);

            mean += (ell - k) as f64 / n_f;
        }

        // ── record mean at half-time ──────────────────────────────────────────
        if !t0_recorded && t > t0 {
            mean_t0 = mean;
            t0_recorded = true;
        }
    }

    (mean - mean_t0) / (t - t0)
}

// =============================================================================
// Record — one row written to the output CSV
// =============================================================================
#[derive(Serialize)]
struct Record {
    lambda: f64,
    omega:  f64,
    #[serde(rename = "T")]
    t_final: f64,
    #[serde(rename = "N")]
    n: usize,
    v_n: f64,
}

// =============================================================================
// main
// =============================================================================
// To confirm version compiled in terminal run [grep -n "FenwickTree" src/main.rs] in IDE terminal
// When running run [cargo build --release] in IDE terminal
fn main() {
    let p_samples: usize = 100;
    let t_final:   f64   = 1000.0;
    let lam:       f64   = 1.0;
    let omegas: &[f64]   = &[0.5];
    let ns: &[usize]     = &[100, 1_000, 10_000, 100_000, 1_000_000];

    let file = "vN_speed_simulations_2.csv";
    let file_exists = Path::new(file).exists();

    let f = OpenOptions::new()
        .create(true)
        .append(true)
        .open(file)
        .expect("Cannot open output CSV");

    let mut wtr = WriterBuilder::new()
        .has_headers(!file_exists)
        .from_writer(BufWriter::new(f));

    for &omega in omegas {
        for &n in ns {
            for p in 0..p_samples {
                println!(
                    "Simulating sample {} of {} samples of N={}",
                    p + 1, p_samples, n
                );
                let start = Instant::now();

                let v_n = simulate_speed(n, lam, omega, t_final);

                wtr.serialize(Record { lambda: lam, omega, t_final, n, v_n })
                   .expect("CSV write failed");
                wtr.flush().expect("CSV flush failed");

                println!(
                    "{}% of {} samples of N={} complete",
                    (p + 1) * 100 / p_samples, p_samples, n
                );
                println!("Runtime: {:.1} seconds", start.elapsed().as_secs_f64());
                println!();
            }
        }
    }

    println!("Results written to {}", file);
}
