// sctransform-rs baseline engine.
//
// Copyright (c) 2024 ORIC Labs (oriclabs), originally contributed to BioLang
// under the MIT license reproduced in THIRD_PARTY_LICENSES.md.
// Modifications Copyright (C) 2026 ORIC Labs.
//
// This file is distributed as part of sctransform-rs under GNU GPL version 3
// only. It began as the BioLang control implementation; identified sampling,
// density, and filtering sections now contain GPL-compatible translations from
// upstream sctransform and R sources, with provenance recorded in comments and
// NOTICE.md.
//
//! Regularized negative-binomial normalization for UMI count matrices.
//!
//! Method references: Hafemeister and Satija (2019), Genome Biology 20:296;
//! Choudhary and Satija (2022), Genome Biology 23:27.

use std::thread;

/// Options, defaulting to the published defaults.
#[derive(Debug, Clone)]
pub struct SctOptions {
    /// Genes to estimate theta on before smoothing. The curve is smooth, so
    /// sampling it does not cost accuracy, and estimating theta is the
    /// expensive step.
    pub genes_for_fit: usize,
    /// Cells used while estimating the model parameters. Seurat's SCTransform
    /// wrapper fits on 5,000 cells by default, then evaluates residuals for the
    /// full dataset.
    pub cells_for_fit: usize,
    /// Multiplies the kernel bandwidth. The reference uses 3; a wider kernel
    /// smooths harder, which is the conservative direction for a regularizer.
    pub bandwidth_adjust: f64,
    /// Residual clip. `None` means `sqrt(n_cells / 30)`, the published default.
    pub clip: Option<f64>,
    /// Genes must be seen in at least this many cells to enter the transform.
    /// SCTransform v2 defaults to five and drops genes below the threshold.
    pub min_cells: usize,
    /// Keep only this many genes, ranked by residual variance. `None` keeps all.
    pub n_variable_features: Option<usize>,
    /// Worker threads. 0 asks the runtime.
    pub threads: usize,
    /// Optional cell-level covariates removed from the Pearson residuals with
    /// an ordinary, non-regularized linear model. Each inner vector is one
    /// covariate and must contain one value per cell. This is deliberately a
    /// second stage: it matches the documented `vars.to.regress` contract
    /// without changing the negative-binomial depth model or HVG ranking.
    pub latent_covariates: Vec<Vec<f64>>,
    /// Center each returned residual column. SCTransform's scale-data contract
    /// centers by default but does not scale to unit variance.
    pub center: bool,
}

impl Default for SctOptions {
    fn default() -> Self {
        Self {
            genes_for_fit: 2000,
            cells_for_fit: 5000,
            bandwidth_adjust: 3.0,
            clip: None,
            min_cells: 5,
            n_variable_features: None,
            threads: 0,
            latent_covariates: Vec::new(),
            center: true,
        }
    }
}

/// A gene-major sparse view: for each gene, the `(cell, count)` pairs that are
/// non-zero. Gene-major because every step here works down a gene at a time.
///
/// Stored flat rather than as a vector per gene. A `Vec<Vec<_>>` costs one
/// allocation and a 24-byte header per gene, and -- the part that actually shows
/// on a UMI matrix -- every column keeps whatever growth slack `push` left it,
/// up to twice the bytes it holds. Sizing the arrays exactly, from a counting
/// pass, removes both. `u32` cell indices halve the index array; a dataset with
/// four billion cells is not the constraint worth planning for here.
pub struct GeneColumns {
    pub n_cells: usize,
    /// `starts[gene]..starts[gene + 1]` indexes `cells` and `counts`.
    /// Length `n_genes + 1`.
    pub starts: Vec<usize>,
    /// Cell index per non-zero, ascending within each gene.
    pub cells: Vec<u32>,
    /// Count per non-zero, parallel to `cells`.
    pub counts: Vec<f64>,
}

impl GeneColumns {
    pub fn n_genes(&self) -> usize {
        self.starts.len().saturating_sub(1)
    }

    /// The non-zeros of one gene: cells ascending, and the counts beside them.
    #[inline]
    pub fn column(&self, gene: usize) -> (&[u32], &[f64]) {
        let (from, to) = (self.starts[gene], self.starts[gene + 1]);
        (&self.cells[from..to], &self.counts[from..to])
    }

    /// How many cells detected this gene.
    #[inline]
    pub fn detected(&self, gene: usize) -> usize {
        self.starts[gene + 1] - self.starts[gene]
    }

    /// Transpose a cell-major scan into this layout.
    ///
    /// `scan` is invoked twice: once to count each gene's non-zeros and once to
    /// place them. Two cheap passes over the source beat one pass that has to
    /// guess every column's size, which is the whole reason this layout exists.
    /// The second pass must emit the same entries in the same order, and cells
    /// must be visited in ascending order so each column comes out sorted.
    pub fn from_cell_major<F>(n_cells: usize, n_genes: usize, mut scan: F) -> Self
    where
        F: FnMut(&mut dyn FnMut(usize, usize, f64)),
    {
        let mut starts = vec![0usize; n_genes + 1];
        scan(&mut |_cell, gene, _count| starts[gene + 1] += 1);
        for gene in 0..n_genes {
            starts[gene + 1] += starts[gene];
        }
        let total = starts[n_genes];
        let mut cells = vec![0u32; total];
        let mut counts = vec![0.0f64; total];
        let mut cursor = starts.clone();
        scan(&mut |cell, gene, count| {
            let position = cursor[gene];
            cells[position] = cell as u32;
            counts[position] = count;
            cursor[gene] = position + 1;
        });
        Self {
            n_cells,
            starts,
            cells,
            counts,
        }
    }
}

#[derive(Debug, Clone)]
pub struct SctResult {
    /// Row-major `n_cells * kept_genes.len()`.
    pub residuals: Vec<f64>,
    /// Indices into the original gene axis, ascending.
    pub kept_genes: Vec<usize>,
    /// Variance of the clipped Pearson residuals before centering and optional
    /// covariate regression, per kept gene, in the same order. This is the
    /// ranking statistic, not a description of the final residual matrix when
    /// second-stage regression is requested.
    pub residual_variance: Vec<f64>,
    /// Regularized theta, per kept gene.
    pub theta: Vec<f64>,
    /// Regularized log-scale intercept, where
    /// `mu_cell,gene = cell_total * exp(intercept_gene)`.
    pub intercept: Vec<f64>,
    /// Original gene indices ordered by decreasing residual variance. Unlike
    /// `kept_genes`, this preserves feature rank and is not matrix-column order.
    pub ranked_genes: Vec<usize>,
    /// Original gene indices used for the raw v2 parameter fit. Exposed for
    /// reproducibility and clean-room conformance diagnostics.
    pub fit_genes: Vec<usize>,
    /// Step-one genes eligible for density-weighted sampling, on the original
    /// feature axis. Exposed for GPL-port conformance diagnostics.
    pub fit_candidates: Vec<usize>,
    /// Unnormalized inverse-density weight parallel to `fit_candidates`.
    pub fit_candidate_weights: Vec<f64>,
    /// Effective `stats::density(..., bw = "nrd")` bandwidth.
    pub sampling_bandwidth: f64,
    /// Original cell indices used for the raw v2 parameter fit.
    pub fit_cells: Vec<usize>,
    /// Unregularized theta per entry of `fit_genes`, as estimated before any
    /// smoothing. `f64::INFINITY` where the gene showed no overdispersion,
    /// matching the convention `theta` uses for Poisson genes.
    ///
    /// Exported because the per-gene estimator and the regularizer are
    /// otherwise confounded: a systematic offset in the returned `theta` could
    /// come from either, and there is no way to tell them apart from the
    /// smoothed values alone.
    pub raw_theta: Vec<f64>,
    /// Unregularized log-scale intercept per entry of `fit_genes`.
    pub raw_intercept: Vec<f64>,
    /// log10 geometric mean per kept gene, the abscissa the regularizer is
    /// fitted against.
    ///
    /// Exported because `theta` is a badly conditioned way to compare two
    /// implementations at low overdispersion. What v2 actually smooths is the
    /// overdispersion factor `od = log10(1 + gmean/theta)`, and inverting it
    /// gives `d ln(theta)/d od = -(10^od * ln10)/(10^od - 1)`, which is about
    /// -25 at the od values a UMI matrix produces. A 0.003 difference in the
    /// smoothed od - negligible in the model - therefore shows up as a 7%
    /// difference in theta. Comparing od directly measures the fit; comparing
    /// theta measures the amplification.
    pub log_geometric_mean: Vec<f64>,
}

/// Center residual columns and optionally regress cell-level covariates out.
///
/// Covariates are centered first, so the intercept is exactly the column mean.
/// A tiny diagonal ridge is used only to make duplicate/constant covariates
/// harmless; it is many orders below the data scale and is not a statistical
/// regularizer.
fn residualize(
    values: &mut [f64],
    n_cells: usize,
    n_genes: usize,
    covariates: &[Vec<f64>],
    center: bool,
) {
    if n_cells == 0 || n_genes == 0 {
        return;
    }
    let mut x: Vec<Vec<f64>> = covariates
        .iter()
        .filter(|column| column.len() == n_cells && column.iter().all(|v| v.is_finite()))
        .cloned()
        .collect();
    for column in &mut x {
        let mean = column.iter().sum::<f64>() / n_cells as f64;
        for value in column {
            *value -= mean;
        }
    }
    let p = x.len();
    let inverse = if p == 0 {
        Vec::new()
    } else {
        let mut augmented = vec![vec![0.0; p * 2]; p];
        let trace = (0..p)
            .map(|j| x[j].iter().map(|value| value * value).sum::<f64>())
            .sum::<f64>();
        let ridge = (trace / p as f64).max(1.0) * 1e-12;
        for row in 0..p {
            for column in 0..p {
                augmented[row][column] = (0..n_cells)
                    .map(|cell| x[row][cell] * x[column][cell])
                    .sum::<f64>();
            }
            augmented[row][row] += ridge;
            augmented[row][p + row] = 1.0;
        }
        for pivot in 0..p {
            let best = (pivot..p)
                .max_by(|&a, &b| {
                    augmented[a][pivot]
                        .abs()
                        .total_cmp(&augmented[b][pivot].abs())
                })
                .unwrap_or(pivot);
            augmented.swap(pivot, best);
            let divisor = augmented[pivot][pivot];
            if divisor.abs() <= 1e-20 {
                continue;
            }
            for column in 0..p * 2 {
                augmented[pivot][column] /= divisor;
            }
            for row in 0..p {
                if row == pivot {
                    continue;
                }
                let factor = augmented[row][pivot];
                for column in 0..p * 2 {
                    augmented[row][column] -= factor * augmented[pivot][column];
                }
            }
        }
        (0..p)
            .map(|row| augmented[row][p..].to_vec())
            .collect::<Vec<_>>()
    };

    if p == 0 && !center {
        return;
    }

    // Two row-major sweeps rather than three strided ones per gene.
    //
    // The gene axis is the *minor* axis of this matrix, so the obvious
    // gene-at-a-time formulation reads every column with a stride of `n_genes`
    // -- a fresh cache line for all but one value in eight, three times over,
    // for each of thousands of genes. Accumulating every gene's statistics
    // together during one sequential pass touches each line once instead.
    //
    // Because the covariates were centered above, `sum_cell x[p][cell]` is zero,
    // so the mean correction drops out of X'y exactly and the two statistics can
    // share a pass.
    let mut means = vec![0.0f64; n_genes];
    let mut xty = vec![0.0f64; p * n_genes];
    for cell in 0..n_cells {
        let row = &values[cell * n_genes..(cell + 1) * n_genes];
        for (mean, value) in means.iter_mut().zip(row) {
            *mean += *value;
        }
        for (predictor, column) in x.iter().enumerate() {
            let weight = column[cell];
            if weight == 0.0 {
                continue;
            }
            let target = &mut xty[predictor * n_genes..(predictor + 1) * n_genes];
            for (accumulator, value) in target.iter_mut().zip(row) {
                *accumulator += weight * *value;
            }
        }
    }
    for mean in &mut means {
        *mean /= n_cells as f64;
    }

    // beta[predictor][gene], laid out to match the sweep below.
    let mut beta = vec![0.0f64; p * n_genes];
    for row in 0..p {
        for (source, coefficient) in inverse[row].iter().enumerate() {
            if *coefficient == 0.0 {
                continue;
            }
            let (target, from) = (row * n_genes, source * n_genes);
            for gene in 0..n_genes {
                beta[target + gene] += coefficient * xty[from + gene];
            }
        }
    }
    drop(xty);

    for cell in 0..n_cells {
        let row = &mut values[cell * n_genes..(cell + 1) * n_genes];
        for (predictor, column) in x.iter().enumerate() {
            let weight = column[cell];
            if weight == 0.0 {
                continue;
            }
            let coefficients = &beta[predictor * n_genes..(predictor + 1) * n_genes];
            for (value, coefficient) in row.iter_mut().zip(coefficients) {
                *value -= weight * coefficient;
            }
        }
        if center {
            for (value, mean) in row.iter_mut().zip(&means) {
                *value -= *mean;
            }
        }
    }
}

/// Digamma, the derivative of log-gamma.
///
/// Recurrence up to a large argument, then the asymptotic series. Needed
/// because the derivative of the negative binomial log-likelihood with respect
/// to theta is written in digammas, and solving for theta means evaluating it
/// once per observation per iteration.
#[cfg(test)]
fn digamma(mut x: f64) -> f64 {
    let mut result = 0.0;
    while x < 6.0 {
        result -= 1.0 / x;
        x += 1.0;
    }
    let inverse = 1.0 / x;
    let square = inverse * inverse;
    result + x.ln()
        - 0.5 * inverse
        - square
            * (1.0 / 12.0
                - square * (1.0 / 120.0 - square * (1.0 / 252.0 - square * (1.0 / 240.0))))
}

/// Derivative of the negative binomial log-likelihood with respect to theta.
///
/// With `mu` fixed, the likelihood in theta alone is
///
/// ```text
///   sum_i [ lgamma(y_i + t) - lgamma(t) + t*log(t) - t*log(t + mu_i)
///           + y_i*log(mu_i) - (y_i + t)*log(t + mu_i) ]      (+ const)
/// ```
///
/// so the score is
///
/// ```text
///   sum_i [ digamma(y_i + t) - digamma(t) + log(t) + 1
///           - log(t + mu_i) - (y_i + t)/(t + mu_i) ]
/// ```
///
/// Zero counts dominate a UMI matrix, and for those `digamma(y + t)` is just
/// `digamma(t)` -- a constant that can be lifted out of the loop rather than
/// recomputed a few hundred million times. That single observation is what
/// makes fitting at this scale practical.
#[cfg(test)]
fn theta_score(counts: &[(usize, f64)], cell_mu: &[f64], theta: f64) -> f64 {
    let digamma_theta = digamma(theta);
    let ln_theta = theta.ln();
    let mut total = 0.0;

    // The part every cell contributes, zero or not.
    for &mu in cell_mu {
        let sum = theta + mu;
        total += ln_theta + 1.0 - sum.ln() - theta / sum;
    }
    // Corrections for the cells that actually observed something.
    for &(cell, count) in counts {
        let sum = theta + cell_mu[cell];
        total += digamma(count + theta) - digamma_theta - count / sum;
    }
    total
}

/// Maximum-likelihood theta for one gene, by bisection on the score.
///
/// Bisection rather than Newton: the score is monotone decreasing in theta over
/// the range that matters, so bisection cannot diverge, and the second
/// derivative costs another digamma evaluation per observation that buys little
/// when the interval is already bracketed.
///
/// Returns `None` when the gene shows no overdispersion at all -- the score
/// stays positive however large theta gets, meaning the counts are Poisson or
/// tighter. That is a real answer, not a failure, and the caller excludes such
/// genes from the smoothing fit rather than letting an arbitrary ceiling drag
/// the curve.
#[cfg(test)]
fn fit_theta(counts: &[(usize, f64)], cell_mu: &[f64]) -> Option<f64> {
    const LOWER: f64 = 1e-3;
    const UPPER: f64 = 1e5;
    const ITERATIONS: usize = 60;

    if theta_score(counts, cell_mu, UPPER) > 0.0 {
        return None;
    }
    if theta_score(counts, cell_mu, LOWER) < 0.0 {
        // Overdispersed beyond the range worth modelling.
        return Some(LOWER);
    }
    let (mut low, mut high) = (LOWER, UPPER);
    for _ in 0..ITERATIONS {
        // Geometric midpoint: theta ranges over orders of magnitude, so
        // bisecting the log is what converges in a predictable number of steps.
        let middle = (low * high).sqrt();
        if theta_score(counts, cell_mu, middle) > 0.0 {
            low = middle;
        } else {
            high = middle;
        }
        if high / low < 1.0 + 1e-6 {
            break;
        }
    }
    Some((low * high).sqrt())
}

/// Gaussian kernel regression, evaluated at `at`.
///
/// The public R `stats::ksmooth` API defines bandwidth by placing the normal
/// kernel quartiles at `+/- bandwidth / 4`; it is not the Gaussian standard
/// deviation.  Since `qnorm(0.75) = 0.67448975...`, the corresponding standard
/// deviation is `bandwidth / (4 * qnorm(0.75))`.
fn kernel_smooth(points: &[(f64, f64)], at: f64, bandwidth: f64) -> f64 {
    const R_KSMOOTH_BANDWIDTH_TO_SIGMA: f64 = 2.697_959_000_784_327;
    let sigma = bandwidth / R_KSMOOTH_BANDWIDTH_TO_SIGMA;
    let mut weight_total = 0.0;
    let mut weighted = 0.0;
    for &(x, y) in points {
        let z = (x - at) / sigma;
        // Match the compact support used by the public R smoother. At four
        // standard deviations the discarded Gaussian weight is below 0.00034.
        if z.abs() > 4.0 {
            continue;
        }
        let weight = (-0.5 * z * z).exp();
        weight_total += weight;
        weighted += weight * y;
    }
    if weight_total > 0.0 {
        weighted / weight_total
    } else {
        // Nothing within reach: fall back to the nearest observation rather
        // than to zero, which would be a silent, wrong answer.
        points
            .iter()
            .min_by(|a, b| (a.0 - at).abs().total_cmp(&(b.0 - at).abs()))
            .map(|&(_, y)| y)
            .unwrap_or(0.0)
    }
}

/// The abscissa `reg_model_pars` evaluates the smoother at.
///
/// Upstream clamps to the range of the genes it fitted:
/// `x_points <- pmin(pmax(genes_log_gmean, min(step1)), max(step1))`. Without
/// it, a gene sparser than every fitted gene is smoothed at its own position,
/// out in a tail where the kernel sees only the few nearest observations --
/// upstream instead reads the curve at its endpoint, which is a flat
/// extrapolation rather than a lopsided average.
///
/// On the HBC control this moves 1,025 of 13,799 genes, 7.43%, all from below,
/// by a median of 0.12 in log10 geometric mean against a bandwidth of 0.55.
#[inline]
fn clamp_to_fitted_range(at: f64, low: f64, high: f64) -> f64 {
    at.max(low).min(high)
}

/// Median of the non-zero UMI entries used by the v2 residual variance floor.
///
/// UMI matrices are integer-valued, so a compact histogram avoids cloning the
/// entire sparse value array. Fractional inputs are accepted by the surrounding
/// matrix API; those rare inputs use an exact comparison-based fallback.
fn nonzero_umi_median(values: &[f64]) -> f64 {
    let positive: usize = values.iter().filter(|value| **value > 0.0).count();
    if positive == 0 {
        return 0.0;
    }
    let integral_max =
        values
            .iter()
            .filter(|value| **value > 0.0)
            .try_fold(0usize, |maximum, value| {
                if value.is_finite()
                    && value.fract() == 0.0
                    && *value <= 1_000_000.0
                    && *value <= usize::MAX as f64
                {
                    Some(maximum.max(*value as usize))
                } else {
                    None
                }
            });
    if let Some(maximum) = integral_max {
        let mut histogram = vec![0usize; maximum + 1];
        for &value in values.iter().filter(|value| **value > 0.0) {
            histogram[value as usize] += 1;
        }
        let lower_rank = (positive - 1) / 2;
        let upper_rank = positive / 2;
        let mut seen = 0usize;
        let mut lower = 0usize;
        let mut upper = 0usize;
        for (value, count) in histogram.into_iter().enumerate().skip(1) {
            let next = seen + count;
            if seen <= lower_rank && lower_rank < next {
                lower = value;
            }
            if seen <= upper_rank && upper_rank < next {
                upper = value;
                break;
            }
            seen = next;
        }
        return (lower as f64 + upper as f64) * 0.5;
    }

    let mut sorted: Vec<f64> = values
        .iter()
        .copied()
        .filter(|value| value.is_finite() && *value > 0.0)
        .collect();
    sorted.sort_by(f64::total_cmp);
    let middle = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[middle - 1] + sorted[middle]) * 0.5
    } else {
        sorted[middle]
    }
}

#[inline]
fn sample_variance_from_moments(sum: f64, sum_squares: f64, count: usize) -> f64 {
    if count < 2 {
        0.0
    } else {
        ((sum_squares - sum * sum / count as f64) / (count - 1) as f64).max(0.0)
    }
}

/// log10 of a gene's geometric mean, from the running sum of `log(count + 1)`
/// over its non-zero cells. Zeros contribute `log(1) = 0` and are skipped.
///
/// The last step is written `exp(m) - 1` rather than `expm1(m)` deliberately,
/// and the choice is not a matter of taste. `expm1` is the *more* accurate of
/// the two for the small `m` a sparse gene produces, which is precisely why it
/// cannot be used: upstream `row_gmean_dgcmatrix` computes `exp(sum/ncol) -
/// eps`, so the more accurate form is the one that disagrees with the
/// reference.
///
/// A one-ULP disagreement here is not cosmetic. This value is the abscissa of
/// the step-one density estimate, so it moves the `bw.nrd` bandwidth, moves
/// every interpolated density, and reweights a sequential sample *without
/// replacement* -- where one flipped draw changes the population every later
/// draw sees. Measured on the HBC control: `expm1` reproduced R's log
/// geometric mean for 8.32% of candidate genes and selected 91.55% of R's fit
/// genes; `exp(m) - 1` reproduces 99.73% and selects 100%.
///
/// GPL-port change derived from upstream commit
/// 49e35b5aeb76a602910207cbfde1561093340be3, `src/utils.cpp`
/// `row_gmean_dgcmatrix`.
#[inline]
fn log10_geometric_mean(log_total: f64, n_cells: usize) -> f64 {
    let geometric_mean = (log_total / n_cells as f64).exp() - 1.0;
    geometric_mean.max(1e-30).log10()
}

/// Current SCTransform regularization works on the overdispersion factor
/// `1 + mean/theta`, on a log10 scale, rather than directly on log10(theta).
#[inline]
fn theta_to_log10_od_factor(log10_mean: f64, theta: f64) -> f64 {
    (1.0 + 10f64.powf(log10_mean) / theta).log10()
}

#[inline]
fn log10_od_factor_to_theta(log10_mean: f64, log10_factor: f64) -> f64 {
    let excess = 10f64.powf(log10_factor) - 1.0;
    if excess > 0.0 && excess.is_finite() {
        10f64.powf(log10_mean) / excess
    } else {
        f64::INFINITY
    }
}

/// Density-dependent gene sample described by Hafemeister and Satija.
///
/// Candidate genes arrive sorted by log10 geometric mean. A Gaussian density
/// estimate gives each gene weight `1 / density(expression)`, so sparse tails
/// are represented alongside the crowded middle. The public R sampling
/// contract applies those weights sequentially without replacement.
struct DensityWeightedSample {
    selected: Vec<usize>,
    candidates: Vec<usize>,
    weights: Vec<f64>,
    bandwidth: f64,
}

fn density_weighted_sample(
    candidates: &[usize],
    log_means: &[f64],
    wanted: usize,
    rng: &mut RMersenneTwister,
) -> DensityWeightedSample {
    if candidates.len() <= wanted {
        return DensityWeightedSample {
            selected: candidates.to_vec(),
            candidates: candidates.to_vec(),
            weights: vec![1.0; candidates.len()],
            bandwidth: f64::NAN,
        };
    }
    let n = candidates.len();
    let candidate_x: Vec<f64> = candidates.iter().map(|gene| log_means[*gene]).collect();
    let mut x = candidate_x.clone();
    x.sort_by(f64::total_cmp);
    let mean = x.iter().sum::<f64>() / n as f64;
    let sd = (x.iter().map(|value| (value - mean).powi(2)).sum::<f64>() / (n - 1) as f64)
        .max(0.0)
        .sqrt();
    let quantile = |probability: f64| -> f64 {
        let position = probability * (n - 1) as f64;
        let low = position.floor() as usize;
        let high = position.ceil() as usize;
        let fraction = position - low as f64;
        x[low] * (1.0 - fraction) + x[high] * fraction
    };
    let iqr = quantile(0.75) - quantile(0.25);
    // Upstream explicitly asks stats::density(..., bw = "nrd"), whose rule is
    // 1.06 * min(sd, IQR / 1.34) * n^(-1/5). The MIT baseline used the nrd0
    // constants (0.9 and 1.349), which changes the weighted sample even though
    // the resulting density curves look nearly identical.
    let spread = if iqr > 0.0 { sd.min(iqr / 1.34) } else { sd };
    let bandwidth = 1.06 * spread * (n as f64).powf(-0.2);
    if bandwidth <= 0.0 || !bandwidth.is_finite() {
        let selected = (0..wanted)
            .map(|index| candidates[index * n / wanted])
            .collect();
        return DensityWeightedSample {
            selected,
            candidates: candidates.to_vec(),
            weights: vec![1.0; candidates.len()],
            bandwidth,
        };
    }

    // stats::density defaults to 512 returned grid points extending three
    // bandwidths beyond the data, after which vst.R linearly interpolates the
    // grid back at each gene. Reproducing that operation matters for sampling:
    // evaluating the KDE directly at every observation is a visually
    // indistinguishable curve but changes enough probability mass to select a
    // different 2,000-gene set under the same random stream.
    //
    // The linear binning below is a Rust translation of R's GPL-compatible
    // stats::density/BinDist path (r-source src/library/stats/R/density.R and
    // src/library/stats/src/massdist.c). At only 512 bins, direct Gaussian
    // convolution is smaller and clearer than carrying an FFT implementation,
    // while evaluating the same mathematical circular convolution.
    const GRID: usize = 512;
    const NORMALIZER: f64 = 0.398_942_280_401_432_7;
    let from = x[0] - 3.0 * bandwidth;
    let to = x[n - 1] + 3.0 * bandwidth;
    let lower = from - 4.0 * bandwidth;
    let upper = to + 4.0 * bandwidth;
    let bin_step = (upper - lower) / (GRID - 1) as f64;
    let mut bins = vec![0.0f64; GRID];
    // Bin in the original gene order, as R's C_BinDist does. Reordering these
    // additions is mathematically harmless but changes the last bits of bins;
    // sequential weighted sampling can amplify those differences.
    for value in &candidate_x {
        let position = (value - lower) / bin_step;
        let bin = position.floor() as usize;
        let fraction = position - bin as f64;
        let weight = 1.0 / n as f64;
        if bin < GRID - 1 {
            bins[bin] += (1.0 - fraction) * weight;
            bins[bin + 1] += fraction * weight;
        } else if bin == GRID - 1 {
            bins[bin] += (1.0 - fraction) * weight;
        }
    }
    let convolved: Vec<f64> = (0..GRID)
        .map(|at| {
            bins.iter()
                .enumerate()
                .map(|(source, weight)| {
                    let z = (source as isize - at as isize) as f64 * bin_step / bandwidth;
                    weight * NORMALIZER / bandwidth * (-0.5 * z * z).exp()
                })
                .sum::<f64>()
                .max(0.0)
        })
        .collect();
    let interpolate = |values: &[f64], origin: f64, step: f64, at: f64| -> f64 {
        let position = ((at - origin) / step).clamp(0.0, (values.len() - 1) as f64);
        let left = position.floor() as usize;
        let right = (left + 1).min(values.len() - 1);
        let fraction = position - left as f64;
        values[left] * (1.0 - fraction) + values[right] * fraction
    };
    let output_step = (to - from) / (GRID - 1) as f64;
    let density_grid: Vec<f64> = (0..GRID)
        .map(|index| {
            interpolate(
                &convolved,
                lower,
                bin_step,
                from + index as f64 * output_step,
            )
        })
        .collect();

    let mut weighted_candidates = Vec::with_capacity(n);
    for (index, value) in candidate_x.iter().copied().enumerate() {
        let density = interpolate(&density_grid, from, output_step, value);
        weighted_candidates.push((1.0 / density.max(f64::MIN_POSITIVE), candidates[index]));
    }

    let diagnostic_weights: Vec<f64> = weighted_candidates
        .iter()
        .map(|(weight, _gene)| *weight)
        .collect();
    let mut selected = weighted_sample_without_replacement(weighted_candidates, wanted, rng);
    selected.sort_by(|left, right| {
        log_means[*left]
            .total_cmp(&log_means[*right])
            .then(left.cmp(right))
    });
    DensityWeightedSample {
        selected,
        candidates: candidates.to_vec(),
        weights: diagnostic_weights,
        bandwidth,
    }
}

/// R's unequal-probability sampling without replacement. This follows
/// `ProbSampleNoReplace` from R's GPL-compatible `src/main/random.c`: sort
/// probabilities descending, draw from remaining mass, and remove the chosen
/// entry by shifting the tail left.
fn weighted_sample_without_replacement(
    mut weighted_candidates: Vec<(f64, usize)>,
    wanted: usize,
    rng: &mut RMersenneTwister,
) -> Vec<usize> {
    // `sample.int(..., prob=)` first normalizes in the original order, then
    // applies R's `revsort` heapsort to probabilities and identities in
    // parallel. Both details are observable here: scaling an unnormalized
    // total is mathematically equivalent but accumulates enough rounding over
    // 2,000 removals to select a different tail of genes.
    let normalization: f64 = weighted_candidates.iter().map(|item| item.0).sum();
    for item in &mut weighted_candidates {
        item.0 /= normalization;
    }
    r_revsort(&mut weighted_candidates);

    let mut total_weight = 1.0;
    let mut selected = Vec::with_capacity(wanted);
    for _ in 0..wanted {
        let target = rng.uniform() * total_weight;
        let mut cumulative = 0.0;
        let mut selected_slot = weighted_candidates.len() - 1;
        for (slot, &(weight, _)) in weighted_candidates
            .iter()
            .take(weighted_candidates.len() - 1)
            .enumerate()
        {
            cumulative += weight;
            if target <= cumulative {
                selected_slot = slot;
                break;
            }
        }
        let (weight, gene) = weighted_candidates.remove(selected_slot);
        total_weight -= weight;
        selected.push(gene);
    }
    selected
}

/// Descending heapsort used by R's public probability-sampling path.
///
/// This is intentionally not Rust's standard sort: equal density weights are
/// common for sparse genes and the exact tie permutation changes the seeded
/// without-replacement sample. Translation of GPL R `revsort` from
/// `src/main/sort.c`.
fn r_revsort(values: &mut [(f64, usize)]) {
    let n = values.len();
    if n <= 1 {
        return;
    }

    let mut l = (n >> 1) + 1;
    let mut ir = n;
    loop {
        let carried;
        if l > 1 {
            l -= 1;
            carried = values[l - 1];
        } else {
            carried = values[ir - 1];
            values[ir - 1] = values[0];
            ir -= 1;
            if ir == 1 {
                values[0] = carried;
                return;
            }
        }

        let mut i = l;
        let mut j = l << 1;
        while j <= ir {
            if j < ir && values[j - 1].0 > values[j].0 {
                j += 1;
            }
            if carried.0 > values[j - 1].0 {
                values[i - 1] = values[j - 1];
                i = j;
                j += i;
            } else {
                j = ir + 1;
            }
        }
        values[i - 1] = carried;
    }
}

/// R-compatible Mersenne-Twister stream used only for the public, fixed-seed
/// SCTransform sampling contract. The generator is the published MT19937
/// recurrence; initialization and discrete draws are pinned by black-box
/// observations from public R APIs rather than package implementation source.
struct RMersenneTwister {
    state: [u32; 624],
    index: usize,
}

impl RMersenneTwister {
    fn new(mut seed: u32) -> Self {
        for _ in 0..50 {
            seed = seed.wrapping_mul(69_069).wrapping_add(1);
        }
        // R's public state has a position word followed by 624 MT words; one
        // LCG result is skipped before those state words are filled.
        seed = seed.wrapping_mul(69_069).wrapping_add(1);
        let mut state = [0u32; 624];
        for word in &mut state {
            seed = seed.wrapping_mul(69_069).wrapping_add(1);
            *word = seed;
        }
        Self { state, index: 624 }
    }

    fn next_u32(&mut self) -> u32 {
        if self.index >= 624 {
            for index in 0..624 {
                let value = (self.state[index] & 0x8000_0000)
                    | (self.state[(index + 1) % 624] & 0x7fff_ffff);
                self.state[index] = self.state[(index + 397) % 624]
                    ^ (value >> 1)
                    ^ if value & 1 == 0 { 0 } else { 0x9908_b0df };
            }
            self.index = 0;
        }
        let mut value = self.state[self.index];
        self.index += 1;
        value ^= value >> 11;
        value ^= (value << 7) & 0x9d2c_5680;
        value ^= (value << 15) & 0xefc6_0000;
        value ^= value >> 18;
        value
    }

    fn uniform(&mut self) -> f64 {
        self.next_u32() as f64 / 4_294_967_296.0
    }

    fn uniform_index(&mut self, upper: usize) -> usize {
        debug_assert!(upper > 0);
        let bits = usize::BITS - (upper - 1).leading_zeros();
        let mask = if bits == usize::BITS {
            usize::MAX
        } else {
            (1usize << bits) - 1
        };
        loop {
            let mut value = 0usize;
            let mut generated = 0u32;
            while generated < bits {
                value = (value << 16) | ((self.next_u32() >> 16) as usize);
                generated += 16;
            }
            value &= mask;
            if value < upper {
                return value;
            }
        }
    }

    fn sample_indices(&mut self, count: usize, wanted: usize) -> Vec<usize> {
        let mut pool: Vec<usize> = (0..count).collect();
        let mut selected = Vec::with_capacity(wanted);
        for remaining in (count - wanted + 1..=count).rev() {
            let slot = self.uniform_index(remaining);
            selected.push(pool[slot]);
            pool[slot] = pool[remaining - 1];
        }
        selected
    }
}

/// Run the transform.
pub fn sctransform(data: &GeneColumns, options: &SctOptions) -> SctResult {
    let n_cells = data.n_cells;
    let n_genes = data.n_genes();
    if n_cells == 0 || n_genes == 0 {
        return SctResult {
            residuals: vec![],
            kept_genes: vec![],
            residual_variance: vec![],
            theta: vec![],
            intercept: vec![],
            ranked_genes: vec![],
            fit_genes: vec![],
            fit_candidates: vec![],
            fit_candidate_weights: vec![],
            sampling_bandwidth: f64::NAN,
            raw_theta: vec![],
            raw_intercept: vec![],
            log_geometric_mean: vec![],
            fit_cells: vec![],
        };
    }

    // -- Margins -----------------------------------------------------------
    //
    // One pass produces every per-gene summary the rest of the routine needs:
    // the total, the sum of squares that decides whether a gene is Poisson, and
    // the log-scale expression level the regularizer is fitted against. The
    // alternative -- recomputing each on demand -- is what a closure invites, and
    // one of these is consumed inside a sort comparator, where "on demand"
    // means O(n log n) rescans of the whole sparse column.
    let mut cell_totals = vec![0.0f64; n_cells];
    let mut gene_totals = vec![0.0f64; n_genes];
    let mut poisson_like = vec![false; n_genes];
    let mut regularization_poisson = vec![false; n_genes];
    // The published regularizer uses the geometric mean of log1p counts. On a
    // sparse column this is exp(sum(log1p(nonzero))/n)-1; zeros contribute
    // exactly zero to the log sum and need not be materialized.
    let mut log_means = vec![0.0f64; n_genes];
    for gene in 0..n_genes {
        let (cells, counts) = data.column(gene);
        let mut total = 0.0;
        let mut squares = 0.0;
        let mut log_total = 0.0;
        for (&cell, &count) in cells.iter().zip(counts) {
            cell_totals[cell as usize] += count;
            total += count;
            squares += count * count;
            log_total += (count + 1.0).ln();
        }
        gene_totals[gene] = total;
        let mean = total / n_cells as f64;
        let variance = sample_variance_from_moments(total, squares, n_cells);
        // Two different Poisson rules, applied at two different stages, and
        // conflating them costs genes in both directions.
        //
        // Step-one sampling uses `variance - mean > 0` alone; adding a low-mean
        // rule there removes rare but genuinely overdispersed genes from the
        // population the density sample draws from.
        //
        // Regularization uses more: `reg_model_pars` builds `all_poisson_genes`
        // as the union of `variance - mean <= 0` and `mean < 0.001`, and fits
        // those with an analytic offset model instead of letting them into the
        // smoother. Measured on the HBC control, the low-mean rule alone
        // accounts for 29 of the 68 genes this port was smoothing that R was
        // not.
        poisson_like[gene] = variance <= mean;
        regularization_poisson[gene] = poisson_like[gene] || mean < 0.001;
        log_means[gene] = log10_geometric_mean(log_total, n_cells);
    }
    let grand_total: f64 = gene_totals.iter().sum::<f64>().max(1e-12);
    let min_variance = (nonzero_umi_median(&data.counts) / 5.0).powi(2);

    // SCT v2 estimates parameters from a fixed-size cell sample. Use a local
    // R-compatible stream for the public fixed-seed sampling contract; it does
    // not touch BioLang's global RNG or depend on thread order. Sorting the
    // retained indices restores sparse-column traversal order for likelihoods.
    let mut sampling_rng = RMersenneTwister::new(1_448_145);
    let fit_cells: Vec<usize> = if options.cells_for_fit == 0 || n_cells <= options.cells_for_fit {
        (0..n_cells).collect()
    } else {
        let mut indices = sampling_rng.sample_indices(n_cells, options.cells_for_fit);
        indices.sort_unstable();
        indices
    };
    let mut fit_cell_lookup = vec![usize::MAX; n_cells];
    for (local, &cell) in fit_cells.iter().enumerate() {
        fit_cell_lookup[cell] = local;
    }
    let fit_cell_totals: Vec<f64> = fit_cells.iter().map(|&cell| cell_totals[cell]).collect();

    // The public vst contract removes genes detected in fewer than min_cells.
    // Keeping them and merely excluding them from the parameter fit changes
    // both the returned feature axis and downstream variable-feature labels.
    let modelled: Vec<usize> = (0..n_genes)
        .filter(|&gene| data.detected(gene) >= options.min_cells)
        .collect();
    if modelled.is_empty() {
        return SctResult {
            residuals: vec![],
            kept_genes: vec![],
            residual_variance: vec![],
            theta: vec![],
            intercept: vec![],
            ranked_genes: vec![],
            fit_genes: vec![],
            fit_candidates: vec![],
            fit_candidate_weights: vec![],
            sampling_bandwidth: f64::NAN,
            raw_theta: vec![],
            raw_intercept: vec![],
            log_geometric_mean: vec![],
            fit_cells: vec![],
        };
    }
    // -- Which genes to estimate theta on ---------------------------------
    //
    // Upstream vst.R 0.4.3 constructs the step-one gene population *after*
    // sampling cells: a gene must be detected in min_cells of the sampled
    // cells, and v2 removes genes whose full-matrix variance does not exceed
    // their mean before density sampling. This is intentionally different from
    // the MIT baseline, which sampled from every full-matrix modelled gene and
    // removed Poisson genes only after fitting. The ordering here is a
    // GPL-port change derived from upstream commit
    // 49e35b5aeb76a602910207cbfde1561093340be3, R/vst.R:217-266.
    let fit_detected: Vec<usize> = (0..n_genes)
        .map(|gene| {
            let (cells, _) = data.column(gene);
            cells
                .iter()
                .filter(|cell| fit_cell_lookup[**cell as usize] != usize::MAX)
                .count()
        })
        .collect();

    // Sample inversely to the density of log-expression, so the regularizer is
    // informed across the full abundance range rather than being dominated by
    // the crowded middle.
    let fit_candidates: Vec<usize> = modelled
        .iter()
        .copied()
        .filter(|&gene| fit_detected[gene] >= options.min_cells && !poisson_like[gene])
        .collect();
    let gene_sample = density_weighted_sample(
        &fit_candidates,
        &log_means,
        options.genes_for_fit,
        &mut sampling_rng,
    );
    let fit_genes = gene_sample.selected;

    // -- Estimate theta on those genes ------------------------------------
    let fitted = parallel_map(options.threads, fit_genes.clone(), {
        let fit_cell_totals = fit_cell_totals.clone();
        let fit_cell_lookup = fit_cell_lookup.clone();
        let columns: Vec<Vec<(usize, f64)>> = fit_genes
            .iter()
            .map(|&gene| {
                let (cells, counts) = data.column(gene);
                cells
                    .iter()
                    .zip(counts)
                    .filter_map(|(&cell, &count)| {
                        let local = fit_cell_lookup[cell as usize];
                        (local != usize::MAX).then_some((local, count))
                    })
                    .collect()
            })
            .collect();
        let totals: Vec<f64> = fit_genes.iter().map(|&g| gene_totals[g]).collect();
        let poisson: Vec<bool> = fit_genes.iter().map(|&gene| poisson_like[gene]).collect();
        move |slot: usize, _gene: usize| {
            let scale = totals[slot] / grand_total;
            let fallback = scale.max(1e-30).ln();
            if poisson[slot] {
                return (None, fallback);
            }
            // `fit_glmGamPoi_offset`, rather than this crate's own Cox-Reid
            // fit. Handed identical inputs the two disagreed by 1.3% at the
            // median; against glmGamPoi the ported chain agrees to 1.3e-10.
            let fit = crate::overdispersion::fit_offset_model(&columns[slot], &fit_cell_totals);
            if fit.theta.is_finite() && fit.intercept.is_finite() {
                (Some(fit.theta), fit.intercept)
            } else {
                (
                    None,
                    if fit.intercept.is_finite() {
                        fit.intercept
                    } else {
                        fallback
                    },
                )
            }
        }
    });

    // -- Regularize --------------------------------------------------------
    //
    // Genes with no detectable overdispersion are dropped from the fit rather
    // than pinned to a ceiling: an arbitrary large value is not an observation,
    // and letting a run of them into the smoother bends the curve toward a
    // number nobody measured.
    //
    // Upstream additionally drops outliers, and it scores them on the *full*
    // step-one set before any filtering -- so the scoring vectors below are
    // built over every fit gene, with an infinite theta contributing an
    // od-factor of exactly zero, which is what `log10(1 + gmean/Inf)` gives.
    // `reg_model_pars` scores every column of `model_pars` and takes `any`; the
    // third column, `log_umi`, is `log(10)` for every gene, so its spread is
    // zero and it can never flag anything.
    let step1_x: Vec<f64> = fit_genes.iter().map(|&gene| log_means[gene]).collect();
    let step1_dispersion: Vec<f64> = fit_genes
        .iter()
        .zip(&fitted)
        .map(|(&gene, &(theta, _))| match theta {
            Some(value) => theta_to_log10_od_factor(log_means[gene], value),
            None => 0.0,
        })
        .collect();
    let step1_intercept: Vec<f64> = fitted.iter().map(|&(_, intercept)| intercept).collect();
    let outlier_flags: Vec<bool> = {
        let by_dispersion = crate::outlier::is_outlier(&step1_dispersion, &step1_x, 10.0);
        let by_intercept = crate::outlier::is_outlier(&step1_intercept, &step1_x, 10.0);
        by_dispersion
            .iter()
            .zip(&by_intercept)
            .map(|(a, b)| *a || *b)
            .collect()
    };
    let is_outlier_gene: std::collections::HashSet<usize> = fit_genes
        .iter()
        .zip(&outlier_flags)
        .filter(|(_, &flagged)| flagged)
        .map(|(&gene, _)| gene)
        .collect();

    let observations: Vec<(f64, f64)> = fit_genes
        .iter()
        .zip(&fitted)
        .filter(|(gene, _)| !regularization_poisson[**gene] && !is_outlier_gene.contains(*gene))
        .filter_map(|(&gene, &(theta, _))| {
            theta.map(|value| {
                (
                    log_means[gene],
                    theta_to_log10_od_factor(log_means[gene], value),
                )
            })
        })
        .collect();

    let intercept_observations: Vec<(f64, f64)> = fit_genes
        .iter()
        .zip(&fitted)
        .filter(|(gene, (theta, _))| {
            !regularization_poisson[**gene] && theta.is_some() && !is_outlier_gene.contains(*gene)
        })
        .map(|(&gene, &(_, intercept))| (log_means[gene], intercept))
        .collect();

    let regularized_theta: Vec<f64> = if observations.len() < 2 {
        // Nothing to smooth against. Poisson is the honest fallback: it is what
        // "no overdispersion could be estimated" actually means.
        vec![f64::INFINITY; modelled.len()]
    } else {
        let x_values: Vec<f64> = observations.iter().map(|&(x, _)| x).collect();
        let bandwidth = crate::bandwidth::bw_sj(&x_values, options.bandwidth_adjust);
        let low = x_values.iter().copied().fold(f64::INFINITY, f64::min);
        let high = x_values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        modelled
            .iter()
            .map(|&gene| {
                if regularization_poisson[gene] {
                    f64::INFINITY
                } else {
                    let at = clamp_to_fitted_range(log_means[gene], low, high);
                    log10_od_factor_to_theta(
                        log_means[gene],
                        kernel_smooth(&observations, at, bandwidth),
                    )
                }
            })
            .collect()
    };

    let regularized_intercept: Vec<f64> = if intercept_observations.len() < 2 {
        modelled
            .iter()
            .map(|&gene| (gene_totals[gene] / grand_total).max(1e-30).ln())
            .collect()
    } else {
        let x_values: Vec<f64> = intercept_observations.iter().map(|&(x, _)| x).collect();
        let bandwidth = crate::bandwidth::bw_sj(&x_values, options.bandwidth_adjust);
        let low = x_values.iter().copied().fold(f64::INFINITY, f64::min);
        let high = x_values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        modelled
            .iter()
            .map(|&gene| {
                // `vst_out_offset`: the analytic model upstream substitutes for
                // every gene in `all_poisson_genes`, which is the wider rule.
                // `gene_total / grand_total` is `amean / mean_cell_sum`, so this
                // is upstream's `log(genes_amean) - log(mean_cell_sum)`.
                if regularization_poisson[gene] {
                    (gene_totals[gene] / grand_total).max(1e-30).ln()
                } else {
                    let at = clamp_to_fitted_range(log_means[gene], low, high);
                    kernel_smooth(&intercept_observations, at, bandwidth)
                }
            })
            .collect()
    };

    // -- Residuals ---------------------------------------------------------
    let clip = options
        .clip
        .unwrap_or_else(|| (n_cells as f64 / 30.0).sqrt());

    let n_modelled = modelled.len();

    // First pass: rank every evaluable gene from scalar accumulators only.
    // Keeping this separate from materialization is what makes a feature cap a
    // bound on the dense allocation itself rather than only on the return value.
    let residual_variance = parallel_map(options.threads, modelled.clone(), |slot, gene| {
        let scale = regularized_intercept[slot].exp();
        let theta = regularized_theta[slot];

        // Sum the whole column as if it were empty, then correct the cells that
        // observed something. The dense part is then a straight-line loop over
        // `cell_totals` with no per-cell test against a sparse cursor and no
        // per-element `theta` branch, which is what lets it vectorize; the
        // correction costs two evaluations per non-zero, and non-zeros are a few
        // percent of a UMI matrix.
        //
        // A zero count always gives a negative residual, so only the lower clip
        // can bind, and `-mu/s` is the same value `(0 - mu)/s` produced before.
        let mut sum = 0.0;
        let mut sum_squares = 0.0;
        if theta.is_finite() {
            for &cell_total in &cell_totals {
                let mu = (cell_total * scale).max(1e-12);
                let residual = (-mu / (mu + mu * mu / theta).max(min_variance).sqrt()).max(-clip);
                sum += residual;
                sum_squares += residual * residual;
            }
        } else {
            for &cell_total in &cell_totals {
                let mu = (cell_total * scale).max(1e-12);
                let residual = (-mu / mu.max(min_variance).sqrt()).max(-clip);
                sum += residual;
                sum_squares += residual * residual;
            }
        }
        let (cells, counts) = data.column(gene);
        for (&cell, &count) in cells.iter().zip(counts) {
            let mu = (cell_totals[cell as usize] * scale).max(1e-12);
            let deviation = if theta.is_finite() {
                (mu + mu * mu / theta).max(min_variance).sqrt()
            } else {
                mu.max(min_variance).sqrt()
            };
            let empty = (-mu / deviation).max(-clip);
            let observed = ((count - mu) / deviation).clamp(-clip, clip);
            sum += observed - empty;
            sum_squares += observed * observed - empty * empty;
        }
        sample_variance_from_moments(sum, sum_squares, n_cells)
    });

    // -- Optionally keep only the most variable genes ---------------------
    let mut ranked: Vec<usize> = (0..n_modelled).collect();
    ranked.sort_by(|&a, &b| {
        residual_variance[b]
            .total_cmp(&residual_variance[a])
            .then(a.cmp(&b))
    });

    if let Some(wanted) = options.n_variable_features.filter(|&n| n < n_modelled) {
        ranked.truncate(wanted);
    }
    let ranked_genes: Vec<usize> = ranked.iter().map(|&slot| modelled[slot]).collect();
    let mut selected = ranked;
    selected.sort_unstable();
    let width = selected.len();

    // Second pass: materialize only the selected columns.
    //
    // Split by *cell*, not by gene. The output is row-major, so a thread that
    // owns a run of cells owns whole rows: it writes them sequentially, the
    // ranges are disjoint at the element and the cache line, and `chunks_mut`
    // hands them out without an unsafe block. Building a column at a time
    // instead means scattering each of its values `width` doubles apart, which
    // on a 3,000-gene result is a fresh cache line for every single store.
    //
    // Cells are visited in ascending order, so one cursor per selected gene
    // walks its sparse column forward and never seeks. The live window is
    // `width` cursors plus the cache lines they point at -- a few hundred
    // kilobytes -- rather than a residual column per worker.
    let scales: Vec<f64> = selected
        .iter()
        .map(|&slot| regularized_intercept[slot].exp())
        .collect();
    let thetas: Vec<f64> = selected
        .iter()
        .map(|&slot| regularized_theta[slot])
        .collect();
    let spans: Vec<(usize, usize)> = selected
        .iter()
        .map(|&slot| {
            let gene = modelled[slot];
            (data.starts[gene], data.starts[gene + 1])
        })
        .collect();

    let mut residuals = vec![0.0f64; n_cells * width];
    if width > 0 {
        let workers = worker_count(options.threads, n_cells);
        let rows_per_worker = n_cells.div_ceil(workers);
        let fill = |first_cell: usize, block: &mut [f64]| {
            // Where each gene's column enters this cell range.
            let mut cursor: Vec<usize> = spans
                .iter()
                .map(|&(from, to)| {
                    from + data.cells[from..to].partition_point(|&c| (c as usize) < first_cell)
                })
                .collect();
            for (offset, row) in block.chunks_mut(width).enumerate() {
                let cell = (first_cell + offset) as u32;
                let cell_total = cell_totals[cell as usize];
                for gene in 0..width {
                    let position = cursor[gene];
                    let count = if position < spans[gene].1 && data.cells[position] == cell {
                        cursor[gene] = position + 1;
                        data.counts[position]
                    } else {
                        0.0
                    };
                    let mu = (cell_total * scales[gene]).max(1e-12);
                    let theta = thetas[gene];
                    let variance = if theta.is_finite() {
                        mu + mu * mu / theta
                    } else {
                        mu
                    }
                    .max(min_variance);
                    row[gene] = ((count - mu) / variance.sqrt()).clamp(-clip, clip);
                }
            }
        };
        if workers == 1 {
            fill(0, &mut residuals);
        } else {
            thread::scope(|scope| {
                for (block_index, block) in
                    residuals.chunks_mut(rows_per_worker * width).enumerate()
                {
                    let fill = &fill;
                    scope.spawn(move || fill(block_index * rows_per_worker, block));
                }
            });
        }
    }
    residualize(
        &mut residuals,
        n_cells,
        width,
        &options.latent_covariates,
        options.center,
    );
    SctResult {
        residuals,
        kept_genes: selected.iter().map(|&slot| modelled[slot]).collect(),
        residual_variance: selected
            .iter()
            .map(|&slot| residual_variance[slot])
            .collect(),
        theta: selected
            .iter()
            .map(|&slot| regularized_theta[slot])
            .collect(),
        intercept: selected
            .iter()
            .map(|&slot| regularized_intercept[slot])
            .collect(),
        ranked_genes,
        raw_theta: fitted
            .iter()
            .map(|&(theta, _)| theta.unwrap_or(f64::INFINITY))
            .collect(),
        raw_intercept: fitted.iter().map(|&(_, intercept)| intercept).collect(),
        log_geometric_mean: selected
            .iter()
            .map(|&slot| log_means[modelled[slot]])
            .collect(),
        fit_genes,
        fit_candidates: gene_sample.candidates,
        fit_candidate_weights: gene_sample.weights,
        sampling_bandwidth: gene_sample.bandwidth,
        fit_cells,
    }
}

/// Map `f` over `items`, in order, across threads.
///
/// A few lines of `std::thread` rather than a dependency: `bio-core` carries
/// only `serde`, and one embarrassingly parallel loop does not justify changing
/// that. Work is split into contiguous chunks, which suits this because the
/// items are sorted by expression and adjacent genes cost about the same.
fn worker_count(threads: usize, count: usize) -> usize {
    let requested = if threads == 0 {
        thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    } else {
        threads
    };
    requested.min(count).max(1)
}

fn parallel_map<T, R, F>(threads: usize, items: Vec<T>, f: F) -> Vec<R>
where
    T: Send + Sync + Copy,
    R: Send,
    F: Fn(usize, T) -> R + Sync,
{
    let count = items.len();
    if count == 0 {
        return vec![];
    }
    let threads = worker_count(threads, count);

    if threads == 1 {
        return items
            .into_iter()
            .enumerate()
            .map(|(index, item)| f(index, item))
            .collect();
    }

    let chunk = count.div_ceil(threads);
    let chunks = thread::scope(|scope| {
        let mut handles = Vec::with_capacity(threads);
        for worker in 0..threads {
            let start = worker * chunk;
            if start >= count {
                break;
            }
            let end = (start + chunk).min(count);
            let items = &items;
            let f = &f;
            handles.push(scope.spawn(move || {
                (start..end)
                    .map(|index| f(index, items[index]))
                    .collect::<Vec<R>>()
            }));
        }
        handles
            .into_iter()
            .map(|handle| match handle.join() {
                Ok(values) => values,
                Err(payload) => std::panic::resume_unwind(payload),
            })
            .collect::<Vec<_>>()
    });
    let mut out = Vec::with_capacity(count);
    for values in chunks {
        out.extend(values);
    }
    debug_assert_eq!(out.len(), count);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The geometric mean must use `exp(m) - 1`, not `expm1(m)`.
    ///
    /// This test exists because the wrong version looks better. `expm1` is the
    /// numerically superior choice on its own terms, a reviewer would be right
    /// to suggest it, and nothing else in this crate would notice the change --
    /// every other test here passes with either form. What it costs is fit-gene
    /// agreement with R: 100% becomes 91.55%, because these last bits steer a
    /// sequential sample without replacement.
    ///
    /// The expectations are R's values for the same inputs, so the assertion
    /// fails on any silent switch back.
    #[test]
    fn geometric_mean_uses_the_upstream_inaccurate_form() {
        // One gene detected in five of ten thousand cells, count 1 each.
        let n_cells = 10_000;
        let log_total = 5.0 * 2.0f64.ln();
        let m = log_total / n_cells as f64;

        // The two forms genuinely disagree at this scale; if they ever stop
        // disagreeing this test proves nothing and should be revisited.
        assert_ne!(
            m.exp() - 1.0,
            m.exp_m1(),
            "the two forms agree here, so this test no longer guards anything"
        );

        assert_eq!(
            log10_geometric_mean(log_total, n_cells),
            -3.460_129_274_946_374,
            "expected R's exp(m) - 1 result; expm1 would give {}",
            m.exp_m1().log10()
        );

        // A second case with mixed counts, to catch a change that happens to
        // land on the right answer for the uniform one.
        let log_total: f64 = [1.0f64, 2.0, 3.0, 4.0, 5.0]
            .iter()
            .map(|count| (count + 1.0).ln())
            .sum();
        assert_eq!(
            log10_geometric_mean(log_total, n_cells),
            -3.181_680_656_395_853_7
        );
    }

    #[test]
    fn digamma_matches_known_values() {
        // digamma(1) = -gamma, digamma(0.5) = -gamma - 2 ln 2.
        const EULER_MASCHERONI: f64 = 0.577_215_664_901_532_9;
        assert!((digamma(1.0) + EULER_MASCHERONI).abs() < 1e-9);
        assert!(
            (digamma(0.5) + EULER_MASCHERONI + 2.0 * 2.0f64.ln()).abs() < 1e-9,
            "digamma(0.5) = {}",
            digamma(0.5)
        );
        // Recurrence: digamma(x+1) - digamma(x) = 1/x.
        for x in [0.3, 1.7, 4.2, 11.0] {
            assert!((digamma(x + 1.0) - digamma(x) - 1.0 / x).abs() < 1e-9);
        }
    }

    /// Poisson counts have no overdispersion, so theta must run away rather
    /// than settle on some finite value. Returning None here is what keeps such
    /// genes out of the smoothing fit.
    #[test]
    fn a_poisson_gene_reports_no_overdispersion() {
        let n_cells = 400;
        let mu = vec![5.0; n_cells];
        // Deterministic counts with variance close to the mean.
        let counts: Vec<(usize, f64)> = (0..n_cells)
            .map(|cell| {
                let jitter = ((cell as f64 * 0.7).sin() * 2.2).round();
                (cell, (5.0 + jitter).max(0.0))
            })
            .filter(|&(_, count)| count > 0.0)
            .collect();
        assert_eq!(fit_theta(&counts, &mu), None);
    }

    /// A gene whose counts swing far wider than Poisson must produce a small,
    /// finite theta -- small theta means heavy overdispersion.
    #[test]
    fn an_overdispersed_gene_gets_a_small_theta() {
        let n_cells = 400;
        let mu = vec![5.0; n_cells];
        // Half the cells at zero, half at 10: variance 25 against mean 5.
        let counts: Vec<(usize, f64)> = (0..n_cells)
            .filter(|cell| cell % 2 == 0)
            .map(|cell| (cell, 10.0))
            .collect();
        let theta = fit_theta(&counts, &mu).expect("should detect overdispersion");
        assert!(
            theta < 20.0,
            "theta {theta} is too large for counts this overdispersed"
        );
    }

    #[test]
    fn kernel_smoothing_follows_the_trend() {
        // y = 2x, sampled; the smoother should track it in the interior.
        let points: Vec<(f64, f64)> = (0..100)
            .map(|i| {
                let x = i as f64 / 10.0;
                (x, 2.0 * x)
            })
            .collect();
        for at in [2.0, 5.0, 7.0] {
            let got = kernel_smooth(&points, at, 0.5);
            assert!(
                (got - 2.0 * at).abs() < 0.2,
                "at {at}: smoothed to {got}, expected about {}",
                2.0 * at
            );
        }
    }

    #[test]
    fn kernel_smoothing_matches_public_r_ksmooth_bandwidth_scale() {
        // Public R 4.5.2 observation:
        // ksmooth(-2:2, c(4,1,0,1,4), "normal", bandwidth=1,
        //         x.points=0)$y == 0.04990962160424732.
        let points = [(-2.0, 4.0), (-1.0, 1.0), (0.0, 0.0), (1.0, 1.0), (2.0, 4.0)];
        let actual = kernel_smooth(&points, 0.0, 1.0);
        assert!((actual - 0.049_909_621_604_247_32).abs() < 1e-5);
    }

    #[test]
    fn sheather_jones_matches_public_r_stats_observations() {
        // Expected values were generated through the public stats::bw.SJ API in
        // R 4.5.2. The direct pair-sum equations avoid R's binned numerical
        // approximation, so agreement is asserted to 0.6%, not bit-for-bit.
        let linear: Vec<f64> = (0..101).map(|i| -2.0 + 5.0 * i as f64 / 100.0).collect();
        let tied: Vec<f64> = (0..20)
            .map(|_| -1.0)
            .chain((0..50).map(|i| -0.5 + i as f64 / 49.0))
            .chain((0..30).map(|_| 2.0))
            .collect();
        for (values, expected) in [
            (linear, 0.545_357_608_768_253_8),
            (tied, 0.086_854_392_248_167_17),
        ] {
            let actual = crate::bandwidth::bw_sj(&values, 1.0);
            assert!(
                (actual / expected - 1.0).abs() < 1e-12,
                "bw.SJ was {actual}, expected {expected}"
            );
        }
    }

    #[test]
    fn density_weighted_sampling_represents_sparse_expression_tails() {
        let candidates: Vec<usize> = (0..1000).collect();
        let mut expression = Vec::with_capacity(1000);
        expression.extend((0..50).map(|_| -5.0));
        expression.extend((0..900).map(|index| -0.5 + index as f64 / 899.0));
        expression.extend((0..50).map(|_| 5.0));
        let mut rng = RMersenneTwister::new(1_448_145);
        let mut rng_again = RMersenneTwister::new(1_448_145);
        let selected = density_weighted_sample(&candidates, &expression, 100, &mut rng);
        let selected_again = density_weighted_sample(&candidates, &expression, 100, &mut rng_again);
        assert_eq!(
            selected.selected, selected_again.selected,
            "the seeded sample changed"
        );
        let tail_count = selected
            .selected
            .iter()
            .filter(|gene| **gene < 50 || **gene >= 950)
            .count();
        assert!(
            tail_count >= 18,
            "density weighting selected only {tail_count} tail genes"
        );
    }

    #[test]
    fn r_sampling_stream_matches_public_r_observations() {
        // Public R 4.5.2 observations from set.seed(1448145).
        let expected_uniforms = [
            0.046_777_348_732_575_774,
            0.741_889_093_304_053,
            0.993_382_635_992_020_4,
            0.127_942_691_324_278_7,
        ];
        let mut rng = RMersenneTwister::new(1_448_145);
        for expected in expected_uniforms {
            assert_eq!(rng.uniform(), expected);
        }

        // sample.int(20, 10, replace=FALSE), converted to zero-based indices.
        let mut rng = RMersenneTwister::new(1_448_145);
        assert_eq!(
            rng.sample_indices(20, 10),
            vec![12, 14, 0, 13, 19, 18, 3, 4, 1, 6]
        );
    }

    #[test]
    fn weighted_sampling_matches_public_r_observation_after_cell_sample() {
        // Public R 4.5.2 observation:
        // set.seed(1448145); sample.int(14847, 5000)
        // sample.int(100, 20, prob=seq_len(100))
        let expected = [
            65, 41, 48, 89, 53, 52, 35, 61, 59, 51, 45, 67, 82, 70, 88, 75, 98, 86, 83, 62,
        ];
        let mut rng = RMersenneTwister::new(1_448_145);
        let _cells = rng.sample_indices(14_847, 5_000);
        let weighted = (1..=100).map(|value| (value as f64, value)).collect();
        let observed = weighted_sample_without_replacement(weighted, 20, &mut rng);
        assert_eq!(observed, expected);
    }

    #[test]
    fn nonzero_umi_median_handles_integer_and_fractional_counts() {
        assert_eq!(nonzero_umi_median(&[0.0, 1.0, 1.0, 2.0, 9.0]), 1.5);
        assert_eq!(nonzero_umi_median(&[0.0, 0.5, 1.5, 4.0]), 1.5);
        assert_eq!(nonzero_umi_median(&[0.0, 0.0]), 0.0);
    }

    #[test]
    fn sample_variance_uses_n_minus_one() {
        // sample variance of [1, 2, 3] is exactly one
        assert!((sample_variance_from_moments(6.0, 14.0, 3) - 1.0).abs() < 1e-12);
        assert_eq!(sample_variance_from_moments(5.0, 25.0, 1), 0.0);
    }

    #[test]
    fn od_factor_transform_round_trips_theta() {
        for (mean, theta) in [(0.01_f64, 0.2_f64), (1.0, 5.0), (25.0, 100.0)] {
            let log_mean = mean.log10();
            let transformed = theta_to_log10_od_factor(log_mean, theta);
            let recovered = log10_od_factor_to_theta(log_mean, transformed);
            assert!((recovered / theta - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn genes_below_min_cells_are_removed_from_the_feature_axis() {
        let data = GeneColumns::from_cell_major(10, 2, |emit| {
            for cell in 0..4 {
                emit(cell, 0, 1.0);
            }
            for cell in 0..5 {
                emit(cell, 1, 1.0);
            }
        });
        let result = sctransform(
            &data,
            &SctOptions {
                min_cells: 5,
                center: false,
                ..Default::default()
            },
        );
        assert_eq!(result.kept_genes, vec![1]);
        assert_eq!(result.residuals.len(), 10);
    }

    fn two_population_data() -> GeneColumns {
        // 200 cells, 60 genes. Genes 0-9 are on in the first half only.
        let n_cells = 200;
        GeneColumns::from_cell_major(n_cells, 60, |emit| {
            for cell in 0..n_cells {
                for gene in 0..60 {
                    let base = if gene < 10 && cell < n_cells / 2 {
                        20.0
                    } else {
                        3.0
                    };
                    let count = base + ((cell * (gene + 3)) % 5) as f64;
                    if count > 0.0 {
                        emit(cell, gene, count);
                    }
                }
            }
        })
    }

    /// The flat layout is only correct if the counting pass and the filling
    /// pass agree, and if each column comes out sorted by cell -- everything
    /// downstream walks columns with a forward-only cursor and would silently
    /// read zeros if they did not.
    #[test]
    fn gene_columns_transpose_round_trips() {
        let data = two_population_data();
        assert_eq!(data.n_genes(), 60);
        assert_eq!(data.starts[60], data.cells.len());
        assert_eq!(data.cells.len(), data.counts.len());
        for gene in 0..60 {
            let (cells, counts) = data.column(gene);
            assert_eq!(cells.len(), data.detected(gene));
            assert!(
                cells.windows(2).all(|pair| pair[0] < pair[1]),
                "gene {gene} came out unsorted"
            );
            for (&cell, &count) in cells.iter().zip(counts) {
                let base = if gene < 10 && (cell as usize) < 100 {
                    20.0
                } else {
                    3.0
                };
                let expected = base + ((cell as usize * (gene + 3)) % 5) as f64;
                assert_eq!(count, expected, "gene {gene} cell {cell}");
            }
        }
    }

    #[test]
    fn residuals_are_dense_and_shaped_correctly() {
        let data = two_population_data();
        let result = sctransform(&data, &SctOptions::default());
        assert_eq!(result.kept_genes.len(), 60);
        assert_eq!(result.residuals.len(), 200 * 60);
        // A zero count still has a non-zero residual; that is the whole point.
        assert!(
            result.residuals.iter().any(|&r| r != 0.0),
            "every residual was zero"
        );
    }

    /// The genes that actually separate the two populations must rank above
    /// the ones that do not. This is the property variable-feature selection
    /// depends on.
    #[test]
    fn structured_genes_outrank_flat_ones() {
        let data = two_population_data();
        let options = SctOptions {
            n_variable_features: Some(10),
            ..Default::default()
        };
        let result = sctransform(&data, &options);
        assert_eq!(result.kept_genes.len(), 10);
        assert!(
            result.kept_genes.iter().all(|&g| g < 10),
            "selection missed the structured genes: {:?}",
            result.kept_genes
        );
    }

    #[test]
    fn residuals_respect_the_clip() {
        let data = two_population_data();
        let options = SctOptions {
            clip: Some(1.5),
            center: false,
            ..Default::default()
        };
        let result = sctransform(&data, &options);
        assert!(
            result.residuals.iter().all(|&r| r.abs() <= 1.5 + 1e-12),
            "a residual escaped the clip"
        );
    }

    #[test]
    fn threading_does_not_change_the_answer() {
        let data = two_population_data();
        let single = sctransform(
            &data,
            &SctOptions {
                threads: 1,
                ..Default::default()
            },
        );
        let many = sctransform(
            &data,
            &SctOptions {
                threads: 8,
                ..Default::default()
            },
        );
        assert_eq!(single.kept_genes, many.kept_genes);
        for (a, b) in single.residuals.iter().zip(&many.residuals) {
            assert!((a - b).abs() < 1e-12, "{a} vs {b}");
        }
    }

    #[test]
    fn residuals_are_centered_by_default() {
        let data = two_population_data();
        let result = sctransform(&data, &SctOptions::default());
        let width = result.kept_genes.len();
        for gene in 0..width {
            let mean = (0..data.n_cells)
                .map(|cell| result.residuals[cell * width + gene])
                .sum::<f64>()
                / data.n_cells as f64;
            assert!(mean.abs() < 1e-10, "gene {gene} mean was {mean}");
        }
    }

    #[test]
    fn second_stage_covariate_regression_removes_linear_signal() {
        let data = two_population_data();
        let covariate: Vec<f64> = (0..data.n_cells)
            .map(|cell| if cell < data.n_cells / 2 { 0.0 } else { 1.0 })
            .collect();
        let result = sctransform(
            &data,
            &SctOptions {
                n_variable_features: Some(10),
                latent_covariates: vec![covariate.clone()],
                ..Default::default()
            },
        );
        let width = result.kept_genes.len();
        for gene in 0..width {
            let dot = (0..data.n_cells)
                .map(|cell| covariate[cell] * result.residuals[cell * width + gene])
                .sum::<f64>();
            assert!(dot.abs() < 1e-8, "gene {gene} retained covariate dot {dot}");
        }
    }

    #[test]
    fn ranked_genes_preserve_variance_order() {
        let data = two_population_data();
        let result = sctransform(
            &data,
            &SctOptions {
                n_variable_features: Some(10),
                ..Default::default()
            },
        );
        assert_eq!(result.ranked_genes.len(), 10);
        assert!(result.ranked_genes.iter().all(|gene| *gene < 10));
        assert_ne!(result.ranked_genes, result.kept_genes);
    }
}
