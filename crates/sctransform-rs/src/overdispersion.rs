// Cox-Reid adjusted overdispersion MLE, matching glmGamPoi's estimator.
//
// Copyright (C) 2026 ORIC Labs.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, version 3 only.
//
// Derived from glmGamPoi `src/overdispersion.cpp`
// (https://github.com/const-ae/glmGamPoi), specifically
// `conventional_loglikelihood_fast`, `conventional_score_function_fast` and
// `conventional_overdispersion_mle`. That file carries its own notice: the code
// was originally taken from DESeq2's `src/DESeq2.cpp` by Michael I. Love and
// Constantin Ahlmann-Eltze and is marked "License: LGPL (>= 3)". LGPL-3
// section 2 permits conveying such a work under GPL-3, which is what this
// repository does. See NOTICE.md.
//
//! What sctransform's v2 path actually fits.
//!
//! `fit_glmGamPoi_offset` calls `glmGamPoi::glm_gp(design = ~1, offset =
//! log(total_umi), size_factors = FALSE)` and reports `theta = 1 /
//! overdispersion`. So the model is one intercept per gene over a fixed
//! per-cell offset, and the number that matters is the Cox-Reid adjusted
//! maximum-likelihood overdispersion given that intercept.

use crate::ln_gamma;

/// glmGamPoi's `cr_correction_factor`.
///
/// Its own comment explains why it is not 1.0: for some combinations of `y`,
/// `mu` and `X` the `lgamma(1/theta)` term and the log-determinant term cancel
/// exactly at large theta, and the estimate runs away to infinity. Shaving the
/// Cox-Reid term by 1% breaks the cancellation. It is a deliberate bias, and
/// dropping it would change every estimate.
const CR_CORRECTION_FACTOR: f64 = 0.99;

/// Bounds `nlminb` is called with, on the log scale.
const LOG_THETA_LOWER: f64 = -36.841_361_487_904_734; // ln(1e-16)
const LOG_THETA_UPPER: f64 = 36.841_361_487_904_734; // ln(1e16)

/// The probe point for "is there any maximum at all".
const LOG_THETA_FAR_LEFT: f64 = -18.420_680_743_952_367; // ln(1e-8)

/// Digamma, accurate to about 1e-15.
///
/// The crate's other `digamma` is a smaller approximation that is fine for the
/// bisection it serves. This one is not interchangeable with it: the estimate
/// returned here is compared against R at 1e-8, and the reference evaluates
/// `Rf_digamma` at full double precision, so an estimator built on a 1e-9
/// digamma could not reach the target however well it converged.
///
/// Recurrence up to `x >= 10`, then the asymptotic series with Bernoulli
/// coefficients.
pub(crate) fn digamma_accurate(x: f64) -> f64 {
    if !x.is_finite() || x <= 0.0 {
        return f64::NAN;
    }
    let mut x = x;
    let mut result = 0.0;
    while x < 10.0 {
        result -= 1.0 / x;
        x += 1.0;
    }
    let inv = 1.0 / x;
    let inv2 = inv * inv;
    // ln(x) - 1/(2x) - sum B_2n / (2n x^2n)
    result + x.ln()
        - 0.5 * inv
        - inv2
            * (1.0 / 12.0
                - inv2
                    * (1.0 / 120.0
                        - inv2
                            * (1.0 / 252.0
                                - inv2
                                    * (1.0 / 240.0
                                        - inv2 * (1.0 / 132.0 - inv2 * 691.0 / 32760.0)))))
}

/// Distinct counts and their frequencies, or `None` when there are too many.
///
/// Mirrors `make_table_if_small(y, stop_if_larger = length(y) / 2)`. On UMI
/// data the table is always built -- a gene takes a handful of distinct values
/// across thousands of cells -- and it changes the arithmetic, not just the
/// speed: the lgamma and digamma sums are then accumulated over distinct
/// values weighted by frequency rather than over cells.
///
/// Upstream collects these in a `std::unordered_map`, so its iteration order is
/// the hash order of the table and neither sorted nor reproducible across
/// implementations. Sorting here is a deliberate divergence: it costs a
/// last-bit difference in the sums that no optimum is sensitive to, and buys a
/// result that does not depend on a hash seed.
fn count_table(y: &[f64]) -> Option<Vec<(f64, f64)>> {
    let limit = y.len() / 2;
    let mut keys: Vec<i64> = Vec::new();
    let mut counts: std::collections::HashMap<i64, f64> = std::collections::HashMap::new();
    for &value in y {
        let key = value as i64;
        let entry = counts.entry(key).or_insert_with(|| {
            keys.push(key);
            0.0
        });
        *entry += 1.0;
        if counts.len() > limit {
            return None;
        }
    }
    keys.sort_unstable();
    Some(keys.into_iter().map(|k| (k as f64, counts[&k])).collect())
}

/// The Cox-Reid term and its first derivative, for a single-column design.
///
/// With `X` a column of ones, `b = X' W X` is the scalar `sum(w)`, its
/// determinant is itself, and the LU decomposition upstream performs is the
/// identity. That collapses the whole matrix path to two sums -- but the
/// details still have to be kept: the `1e-50` floor before the log, and the
/// `1e-6` ridge added before inverting.
struct CoxReid {
    /// `-0.5 * log(det(b)) * cr_correction_factor`
    term: f64,
    /// `-0.5 * trace(b^-1 db) * cr_correction_factor`
    derivative: f64,
}

fn cox_reid(mu: &[f64], theta: f64) -> CoxReid {
    let mut b = 0.0;
    let mut db = 0.0;
    for &m in mu {
        let w = 1.0 / (1.0 / m + theta);
        b += w;
        db -= w * w;
    }
    let log_det = if b < 1e-50 { 1e-50_f64.ln() } else { b.ln() };
    CoxReid {
        term: -0.5 * log_det * CR_CORRECTION_FACTOR,
        derivative: -0.5 * (db / (b + 1e-6)) * CR_CORRECTION_FACTOR,
    }
}

/// `conventional_loglikelihood_fast` for a single-column design.
pub fn conventional_loglikelihood(
    y: &[f64],
    mu: &[f64],
    log_theta: f64,
    table: Option<&[(f64, f64)]>,
    do_cr_adjustment: bool,
) -> f64 {
    let theta = log_theta.exp();
    let theta_neg1 = 1.0 / theta;
    let cr_term = if do_cr_adjustment {
        cox_reid(mu, theta).term
    } else {
        0.0
    };

    let mut lgamma_term: f64 = match table {
        Some(entries) => entries
            .iter()
            .map(|&(count, frequency)| frequency * ln_gamma(count + theta_neg1))
            .sum(),
        None => y.iter().map(|&count| ln_gamma(count + theta_neg1)).sum(),
    };
    lgamma_term -= y.len() as f64 * ln_gamma(theta_neg1);

    let mut ll_part = 0.0;
    for (&count, &m) in y.iter().zip(mu) {
        ll_part += (-count - theta_neg1) * (m + theta_neg1).ln();
    }
    ll_part -= y.len() as f64 * theta_neg1 * theta.ln();

    lgamma_term + ll_part + cr_term
}

/// `conventional_score_function_fast` for a single-column design.
///
/// The guards in the `mu * theta` loop are upstream's and are load-bearing near
/// zero, where `log(1 + x) - x/(1 + x)` loses every significant digit to
/// cancellation; upstream brackets the value with the Taylor bounds instead.
/// The `digamma_term` branch is the same idea from the other side: for very
/// large `1/theta` the digamma difference reaches `sum(y)` before the rest of
/// the expression does, so upstream clamps it there minus the first Laurent
/// term.
pub fn conventional_score(
    y: &[f64],
    mu: &[f64],
    log_theta: f64,
    table: Option<&[(f64, f64)]>,
    do_cr_adjustment: bool,
) -> f64 {
    let theta = log_theta.exp();
    let theta_neg1 = 1.0 / theta;
    let cr_term = if do_cr_adjustment {
        cox_reid(mu, theta).derivative
    } else {
        0.0
    };

    let mut digamma_term = 0.0;
    let mut max_y = 0.0f64;
    let mut sum_y = 0.0;
    let mut sum_prod_y = 0.0;
    match table {
        Some(entries) => {
            for &(count, frequency) in entries {
                digamma_term += frequency * digamma_accurate(count + theta_neg1);
                sum_y += frequency * count;
                sum_prod_y += frequency * (count - 1.0) * count;
                max_y = max_y.max(count);
            }
        }
        None => {
            for &count in y {
                digamma_term += digamma_accurate(count + theta_neg1);
                sum_y += count;
                sum_prod_y += (count - 1.0) * count;
                max_y = max_y.max(count);
            }
        }
    }
    let correction = if theta_neg1 > 1e5 {
        sum_prod_y / (2.0 * theta_neg1)
    } else {
        0.0
    };
    if max_y * 1e6 < theta_neg1 {
        digamma_term = sum_y - correction;
    } else {
        digamma_term -= y.len() as f64 * digamma_accurate(theta_neg1);
        digamma_term *= theta_neg1;
        digamma_term = digamma_term.min(sum_y - correction);
    }

    let mut ll_part = 0.0;
    for (&count, &m) in y.iter().zip(mu) {
        let mu_theta = m * theta;
        if mu_theta < 1e-10 {
            ll_part += mu_theta * mu_theta * (1.0 / (1.0 + mu_theta) - 0.5);
        } else if mu_theta < 1e-4 {
            let inv = 1.0 / (1.0 + mu_theta);
            let upper_bound = mu_theta * mu_theta * inv;
            let lower_bound = mu_theta * mu_theta * (inv - 0.5);
            let suggest = (1.0 + mu_theta).ln() - m / (m + theta_neg1);
            ll_part += suggest.min(upper_bound).max(lower_bound);
        } else {
            ll_part += (1.0 + mu_theta).ln() - m / (m + theta_neg1);
        }
        ll_part += count / (m + theta_neg1);
    }
    ll_part *= theta_neg1;

    ll_part - digamma_term + cr_term * theta
}

/// One gene's offset-model fit: what `fit_glmGamPoi_offset` reports.
#[derive(Debug, Clone, Copy)]
pub struct OffsetFit {
    /// `1 / overdispersion`, infinite when none was detected.
    pub theta: f64,
    /// The intercept, on the natural-log scale.
    pub intercept: f64,
}

/// `glm_gp(design = ~1, offset = log(total_umi), size_factors = FALSE)`, for
/// one gene.
///
/// The three stages are upstream's and must run in this order, because each
/// one's output is the next one's input: a moment estimate of the dispersion,
/// a Newton-Raphson fit of the intercept holding that dispersion fixed, and
/// the Cox-Reid MLE of the overdispersion holding the resulting mu fixed.
///
/// The intercept returned here is the *first* one. `glm_gp` fits beta again
/// after shrinking the dispersions and returns that second fit, so this is not
/// yet the intercept sctransform records -- but it is exactly the intercept
/// that defines the mu the overdispersion was estimated from, which is the
/// number this stage exists to get right.
///
/// `counts` is sparse: `(cell, count)` for the non-zeros only. `cell_totals`
/// is dense and gives the offset for every cell, zeros included.
pub fn fit_offset_model(counts: &[(usize, f64)], cell_totals: &[f64]) -> OffsetFit {
    let n_cells = cell_totals.len();
    if n_cells == 0 || counts.is_empty() {
        return OffsetFit {
            theta: f64::INFINITY,
            intercept: f64::NEG_INFINITY,
        };
    }

    // `estimate_dispersions_by_moment`: (rowVars(Y) - bm/mean(colSums)) / bm^2,
    // negatives and NaNs floored to zero by `estimate_dispersions_roughly`.
    let mut sum = 0.0;
    let mut sum_squares = 0.0;
    for &(_, count) in counts {
        sum += count;
        sum_squares += count * count;
    }
    let cells = n_cells as f64;
    let bm = sum / cells;
    let bv = if n_cells > 1 {
        (sum_squares - cells * bm * bm) / (cells - 1.0)
    } else {
        0.0
    };
    let mean_total = cell_totals.iter().sum::<f64>() / cells;
    let moment = (bv - bm / mean_total) / (bm * bm);
    let dispersion_init = if moment.is_finite() && moment > 0.0 {
        moment
    } else {
        0.0
    };

    // `estimate_betas_roughly_group_wise`: log(mean(Y / exp(offset))).
    let offsets: Vec<f64> = cell_totals.iter().map(|total| total.ln()).collect();
    let normalised: f64 = counts
        .iter()
        .map(|&(cell, count)| count / cell_totals[cell])
        .sum::<f64>()
        / cells;
    let mut beta = normalised.ln();
    if !beta.is_finite() {
        return OffsetFit {
            theta: f64::INFINITY,
            intercept: f64::NEG_INFINITY,
        };
    }

    // `fitBeta_one_group`: Newton-Raphson, tolerance 1e-8, at most 100 steps.
    // The zero cells cannot be skipped -- they contribute `-mu / (1 + mu
    // theta)` to the score, which is most of it for a sparse gene.
    let mut dense = vec![0.0f64; n_cells];
    for &(cell, count) in counts {
        dense[cell] = count;
    }
    for _ in 0..100 {
        let mut dl = 0.0;
        let mut ddl = 0.0;
        for (cell, &offset) in offsets.iter().enumerate() {
            let count = dense[cell];
            let mu = (beta + offset).exp();
            let denominator = 1.0 + mu * dispersion_init;
            dl += (count - mu) / denominator;
            ddl += mu * (1.0 + count * dispersion_init) / denominator / denominator;
        }
        let step = dl / ddl;
        beta += step;
        if step.abs() < 1e-8 || beta.is_nan() {
            break;
        }
    }
    beta = beta.max(-1e8);

    let mu: Vec<f64> = offsets.iter().map(|offset| (beta + offset).exp()).collect();
    let estimate = conventional_overdispersion_mle(&dense, &mu, true);

    OffsetFit {
        theta: estimate.theta(),
        intercept: beta,
    }
}

/// The objective at `log_theta`, prepared exactly as the estimator prepares it.
///
/// Callers comparing two estimates must evaluate them on the same function, and
/// that is easier to get wrong than it looks: the zero-mean substitution and
/// the count table both change the arithmetic, so evaluating with a
/// hand-rolled table -- or none -- scores the answers against a function
/// neither optimiser used.
pub fn objective_at(
    y: &[f64],
    mean_vector: &[f64],
    log_theta: f64,
    do_cox_reid_adjustment: bool,
) -> f64 {
    let mu: Vec<f64> = mean_vector
        .iter()
        .map(|&m| if m == 0.0 { 1e-6 } else { m })
        .collect();
    let table = count_table(y);
    conventional_loglikelihood(y, &mu, log_theta, table.as_deref(), do_cox_reid_adjustment)
}

/// What the estimator returned, and why.
#[derive(Debug, Clone, PartialEq)]
pub struct OverdispersionEstimate {
    /// glmGamPoi's `overdispersions` entry. Zero means no overdispersion was
    /// identified, which sctransform turns into an infinite theta.
    pub overdispersion: f64,
    pub message: &'static str,
}

impl OverdispersionEstimate {
    /// `theta = 1 / overdispersion`, infinite when the estimate is zero.
    pub fn theta(&self) -> f64 {
        if self.overdispersion > 0.0 {
            1.0 / self.overdispersion
        } else {
            f64::INFINITY
        }
    }
}

/// `conventional_overdispersion_mle`.
///
/// Upstream optimises with `nlminb`, R's interface to the PORT library's
/// bounded quasi-Newton routine. This does not reproduce PORT's iterate
/// sequence, and does not need to: the maximum is a property of the objective,
/// not of the search. What must match exactly is the objective, which is why
/// the two functions above are transcribed rather than re-derived. Root-finding
/// on the score converges far tighter than `nlminb`'s default `x.tol` of
/// 1.5e-8, so the residual disagreement with R is R's own convergence slack.
///
/// The early exits are upstream's and are not optimisations: all-zero counts
/// and "no maximum even at tiny theta" both return an overdispersion of exactly
/// zero, which is how a gene comes to have infinite theta and be dropped from
/// the regularisation fit.
pub fn conventional_overdispersion_mle(
    y: &[f64],
    mean_vector: &[f64],
    do_cox_reid_adjustment: bool,
) -> OverdispersionEstimate {
    debug_assert_eq!(y.len(), mean_vector.len());
    if y.iter().all(|&count| count == 0.0) {
        return OverdispersionEstimate {
            overdispersion: 0.0,
            message: "All counts y are 0.",
        };
    }

    // Upstream mutates its copy: `mean_vector[mean_vector == 0] <- 1e-06`.
    let mu: Vec<f64> = mean_vector
        .iter()
        .map(|&m| if m == 0.0 { 1e-6 } else { m })
        .collect();
    let table = count_table(y);
    let table = table.as_deref();

    let far_left = conventional_score(y, &mu, LOG_THETA_FAR_LEFT, table, do_cox_reid_adjustment);
    if far_left < 0.0 {
        return OverdispersionEstimate {
            overdispersion: 0.0,
            message: "Even for very small theta, no maximum identified",
        };
    }

    // The score is positive at the far-left probe. If it is still positive at
    // the upper bound the optimum is the bound itself, which is what `nlminb`
    // would report too.
    let far_right = conventional_score(y, &mu, LOG_THETA_UPPER, table, do_cox_reid_adjustment);
    if far_right > 0.0 {
        return OverdispersionEstimate {
            overdispersion: LOG_THETA_UPPER.exp(),
            message: "Optimum at upper bound",
        };
    }

    // Bisect the score's sign change. The bracket is guaranteed by the two
    // probes above, so this cannot fail to converge; 200 halvings take the
    // interval far below the last bit of a double at this scale.
    let (mut low, mut high) = (LOG_THETA_FAR_LEFT.max(LOG_THETA_LOWER), LOG_THETA_UPPER);
    for _ in 0..200 {
        let middle = 0.5 * (low + high);
        if middle <= low || middle >= high {
            break;
        }
        if conventional_score(y, &mu, middle, table, do_cox_reid_adjustment) > 0.0 {
            low = middle;
        } else {
            high = middle;
        }
    }

    // The estimate is the score's root, not a maximiser of the loglikelihood
    // found by function values. That distinction was measured, not assumed.
    //
    // `nlminb` is handed an analytic gradient, so it converges where that
    // gradient vanishes. Polishing the root by golden-section maximisation of
    // the loglikelihood -- the obvious "more correct" thing to do, since the
    // objective is what `nlminb` nominally minimises -- moved the median
    // disagreement with glmGamPoi from 1.3e-10 to 7.5e-7, four thousand times
    // worse. Upstream's score carries clamps and Taylor brackets that make it a
    // deliberately inexact derivative of its own loglikelihood, so the two have
    // slightly different stationary points, and the reference sits at the
    // gradient's.
    OverdispersionEstimate {
        overdispersion: (0.5 * (low + high)).exp(),
        message: "",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// R's `digamma` at points chosen to exercise both the recurrence and the
    /// asymptotic series. Expectations are R's values.
    #[test]
    fn digamma_matches_r_to_full_precision() {
        for (x, want) in [
            (0.1_f64, -10.423_754_940_411_076_f64),
            (1.0, -0.577_215_664_901_532_31),
            (2.5, 0.703_156_640_645_243_41),
            // Just below the recurrence threshold and just on it, so a change
            // to that boundary cannot pass unnoticed.
            (9.999, 2.251_647_417_205_735),
            (10.0, 2.251_752_589_066_721_4),
            (1000.0, 6.907_255_195_648_812),
        ] {
            let got = digamma_accurate(x);
            let tolerance = 1e-14 * want.abs().max(1.0);
            assert!(
                (got - want).abs() < tolerance,
                "digamma({x}) = {got}, R says {want}, off by {}",
                got - want
            );
        }
    }

    /// The table must summarise counts, and must refuse to when there are too
    /// many distinct values -- that refusal is what keeps the estimator correct
    /// on non-count data rather than merely fast on counts.
    #[test]
    fn the_count_table_summarises_only_when_it_is_small() {
        let y = [0.0, 0.0, 1.0, 2.0, 1.0, 0.0, 3.0, 0.0];
        let table = count_table(&y).expect("four distinct values in eight is small enough");
        assert_eq!(table, vec![(0.0, 4.0), (1.0, 2.0), (2.0, 1.0), (3.0, 1.0)]);

        let distinct: Vec<f64> = (0..10).map(|value| value as f64).collect();
        assert!(
            count_table(&distinct).is_none(),
            "ten distinct values in ten is over the length/2 limit"
        );
    }

    /// An all-zero gene has no overdispersion to find, and upstream says so
    /// before touching the optimiser.
    #[test]
    fn an_all_zero_gene_reports_zero_overdispersion() {
        let y = vec![0.0; 32];
        let mu = vec![0.5; 32];
        let estimate = conventional_overdispersion_mle(&y, &mu, true);
        assert_eq!(estimate.overdispersion, 0.0);
        assert!(estimate.theta().is_infinite());
        assert_eq!(estimate.message, "All counts y are 0.");
    }

    /// The score must actually vanish at the returned estimate. This is the
    /// property the bisection claims and the one that makes matching `nlminb`'s
    /// path unnecessary.
    #[test]
    fn the_returned_estimate_is_a_stationary_point() {
        let mu: Vec<f64> = (0..500).map(|cell| 1.0 + (cell % 7) as f64 * 0.3).collect();
        let y: Vec<f64> = (0..500)
            .map(|cell| ((cell * 37 % 11) as f64 - 2.0).max(0.0))
            .collect();
        let estimate = conventional_overdispersion_mle(&y, &mu, true);
        assert!(
            estimate.overdispersion > 0.0,
            "expected a finite estimate, got {estimate:?}"
        );
        let table = count_table(&y);
        let score = conventional_score(
            &y,
            &mu,
            estimate.overdispersion.ln(),
            table.as_deref(),
            true,
        );
        assert!(
            score.abs() < 1e-6,
            "score at the estimate was {score}, so it is not stationary"
        );
    }
}
