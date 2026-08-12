// Outlier detection for model parameters, matching sctransform's `is_outlier`.
//
// Copyright (C) 2026 ORIC Labs.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, version 3 only.
//
// Derived from `satijalab/sctransform` commit
// 49e35b5aeb76a602910207cbfde1561093340be3, R/utils.R (`is_outlier`,
// `robust_scale_binned`, `robust_scale`), GPL-3. Also reproduces R's
// `seq.default` and `cut.default` binning conventions, from the R statistical
// computing environment (GPL-2 or later, conveyed here under GPL-3).
//
//! Why an outlier rule needs to be this literal.
//!
//! `reg_model_pars` drops outlying genes before smoothing, and outliers are by
//! construction the extreme values -- so which genes the rule catches moves the
//! regularisation curve where they sit. The rule scores each gene against the
//! other genes in its expression bin, twice, on two bin grids offset by half a
//! bin width, and keeps the *smaller* of the two absolute scores. A gene has to
//! look extreme under both griddings to count, which is what stops a bin
//! boundary falling in an awkward place from inventing outliers.

use crate::bandwidth::bw_sj;

/// `.Machine$double.eps * 10`, as `is_outlier` spells it.
const EPS_10: f64 = f64::EPSILON * 10.0;

/// R's `mad` constant.
const MAD_CONSTANT: f64 = 1.4826;

/// Median with R's even-length convention: the mean of the two middle values.
fn median_sorted(sorted: &[f64]) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return f64::NAN;
    }
    if n % 2 == 0 {
        (sorted[n / 2 - 1] + sorted[n / 2]) * 0.5
    } else {
        sorted[n / 2]
    }
}

fn median(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    median_sorted(&sorted)
}

/// R's `mad(x)`: `median(|x - median(x)|) * 1.4826`.
fn mad(values: &[f64]) -> f64 {
    let centre = median(values);
    let deviations: Vec<f64> = values.iter().map(|v| (v - centre).abs()).collect();
    median(&deviations) * MAD_CONSTANT
}

/// R's `seq(from, to, by)`.
///
/// Not `while (x <= to)`. R computes the count once as
/// `as.integer((to - from) / by + 1e-10)` and then clamps the last element to
/// `to`, so the final break can be shorter than `by`. Generating breaks by
/// accumulation instead would drift and could add or drop a whole bin.
fn r_seq(from: f64, to: f64, by: f64) -> Vec<f64> {
    if by == 0.0 || !by.is_finite() {
        return vec![from];
    }
    let span = (to - from) / by;
    if span < 0.0 {
        return vec![from];
    }
    let n = (span + 1e-10) as usize;
    (0..=n)
        .map(|index| {
            let value = from + index as f64 * by;
            if by > 0.0 {
                value.min(to)
            } else {
                value.max(to)
            }
        })
        .collect()
}

/// R's `cut(x, breaks)` with the default `right = TRUE`: the intervals are
/// `(breaks[i], breaks[i + 1]]`, and anything outside is `NA`.
fn bin_index(value: f64, breaks: &[f64]) -> Option<usize> {
    if breaks.len() < 2 || !(value > breaks[0]) || value > breaks[breaks.len() - 1] {
        return None;
    }
    // Largest i with breaks[i] < value.
    let mut low = 0usize;
    let mut high = breaks.len() - 1;
    while low + 1 < high {
        let middle = (low + high) / 2;
        if breaks[middle] < value {
            low = middle;
        } else {
            high = middle;
        }
    }
    Some(low)
}

/// `robust_scale_binned`: `(y - median) / (mad + eps)` within each bin of `x`.
///
/// Upstream reaches this through `aggregate` and `order`, which is R's idiom
/// for a grouped transform; the result is the same as scoring each element
/// against its own bin, and genes that fall in no bin keep a score of zero.
fn robust_scale_binned(y: &[f64], x: &[f64], breaks: &[f64]) -> Vec<f64> {
    let mut members: std::collections::HashMap<usize, Vec<usize>> =
        std::collections::HashMap::new();
    for (index, &value) in x.iter().enumerate() {
        if let Some(bin) = bin_index(value, breaks) {
            members.entry(bin).or_default().push(index);
        }
    }
    let mut score = vec![0.0f64; x.len()];
    for indices in members.values() {
        let values: Vec<f64> = indices.iter().map(|&index| y[index]).collect();
        let centre = median(&values);
        let spread = mad(&values) + f64::EPSILON;
        for (&index, value) in indices.iter().zip(&values) {
            score[index] = (value - centre) / spread;
        }
    }
    score
}

/// sctransform's `is_outlier(y, x, th = 10)`.
///
/// `x` is the log10 geometric mean of every step-one gene -- the full set, not
/// the surviving one, because upstream scores before it filters.
pub fn is_outlier(y: &[f64], x: &[f64], threshold: f64) -> Vec<bool> {
    let n = x.len();
    if n == 0 {
        return Vec::new();
    }
    let (mut min, mut max) = (f64::INFINITY, f64::NEG_INFINITY);
    for &value in x {
        min = min.min(value);
        max = max.max(value);
    }
    // The unadjusted bandwidth: `bw_adjust` belongs to the smoother, not here.
    let bin_width = (max - min) * bw_sj(x, 1.0) / 2.0;
    if !(bin_width > 0.0) {
        return vec![false; n];
    }

    let breaks1 = r_seq(min - EPS_10, max + bin_width, bin_width);
    let breaks2 = r_seq(min - EPS_10 - bin_width / 2.0, max + bin_width, bin_width);
    let score1 = robust_scale_binned(y, x, &breaks1);
    let score2 = robust_scale_binned(y, x, &breaks2);

    score1
        .iter()
        .zip(&score2)
        .map(|(a, b)| a.abs().min(b.abs()) > threshold)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_and_mad_follow_r() {
        assert_eq!(median(&[3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median(&[4.0, 1.0, 3.0, 2.0]), 2.5);
        // R: mad(c(1,2,3,4,100)) == 1.4826
        assert!((mad(&[1.0, 2.0, 3.0, 4.0, 100.0]) - 1.4826).abs() < 1e-12);
    }

    /// R counts the breaks with a 1e-10 fudge and then clamps the last one to
    /// `to`. The obvious alternative -- accumulate `by` while the value is at
    /// most `to` -- drops a whole bin whenever the division lands just below an
    /// integer, and bins are what the outlier scores are computed within.
    ///
    /// `seq(0, 0.6, by = 0.2)` is exactly that case: three additions of 0.2
    /// give 0.6000000000000001, which is greater than 0.6, so accumulation
    /// stops at three breaks. R returns four, the last clamped back to 0.6.
    #[test]
    fn r_seq_counts_breaks_the_way_r_does() {
        let got = r_seq(0.0, 0.6, 0.2);
        assert_eq!(got.len(), 4, "got {got:?}");
        assert_eq!(got[3], 0.6, "last break should be clamped to `to`");
        assert!(
            0.2 + 0.2 + 0.2 > 0.6,
            "the premise of this test no longer holds on this platform"
        );

        // A case with no clamping, to pin the ordinary behaviour too.
        let plain = r_seq(0.0, 1.0, 0.3);
        assert_eq!(plain.len(), 4, "got {plain:?}");
        assert!(
            (plain[3] - 0.9).abs() < 1e-15,
            "last break was {}",
            plain[3]
        );
    }

    /// Intervals are left-open and right-closed, and a value exactly on the
    /// first break is outside every bin.
    #[test]
    fn cut_is_left_open_and_right_closed() {
        let breaks = [0.0, 1.0, 2.0];
        assert_eq!(bin_index(0.0, &breaks), None);
        assert_eq!(bin_index(0.5, &breaks), Some(0));
        assert_eq!(bin_index(1.0, &breaks), Some(0));
        assert_eq!(bin_index(1.5, &breaks), Some(1));
        assert_eq!(bin_index(2.0, &breaks), Some(1));
        assert_eq!(bin_index(2.5, &breaks), None);
    }

    /// A single wild value among well-behaved ones must be caught, and its
    /// neighbours must not be. The two-grid rule is what makes the second part
    /// true: on one gridding the outlier can drag its bin's median far enough
    /// to make a neighbour look extreme, and the other gridding vetoes it.
    #[test]
    fn a_lone_extreme_value_is_the_only_outlier() {
        let x: Vec<f64> = (0..200).map(|i| i as f64 / 100.0).collect();
        let mut y: Vec<f64> = x.iter().map(|v| 0.5 * v).collect();
        y[100] = 50.0;
        let flags = is_outlier(&y, &x, 10.0);
        assert!(flags[100], "the planted outlier was not caught");
        assert_eq!(
            flags.iter().filter(|&&f| f).count(),
            1,
            "caught {} outliers, expected exactly one",
            flags.iter().filter(|&&f| f).count()
        );
    }

    /// A constant column has zero spread, so every score is 0/eps = 0 and
    /// nothing is an outlier. sctransform relies on this: it scores the
    /// `log_umi` column, which is `log(10)` for every gene.
    #[test]
    fn a_constant_column_yields_no_outliers() {
        let x: Vec<f64> = (0..100).map(|i| i as f64 / 50.0).collect();
        let y = vec![std::f64::consts::LN_10; 100];
        assert!(is_outlier(&y, &x, 10.0).iter().all(|&f| !f));
    }
}
