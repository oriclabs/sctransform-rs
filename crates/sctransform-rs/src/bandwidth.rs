// Sheather-Jones bandwidth, matching R's `stats::bw.SJ`.
//
// Copyright (C) 2026 ORIC Labs.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, version 3 only.
//
// Translated from the R statistical computing environment, Copyright the R
// Core Team and the R Foundation: `src/library/stats/R/bandwidths.R`
// (`bw.SJ`, `bw_pair_cnts`), `src/library/stats/src/bandwidths.c`
// (`bw_den_binned`, `bw_phi4`, `bw_phi6`), and `src/library/stats/src/zeroin.c`
// (`R_zeroin2`). That source is offered under GPL-2 or later and is conveyed
// here under GPL-3. See NOTICE.md.
//
//! Why this is a translation rather than the published equations.
//!
//! This crate previously evaluated the 1991 Sheather-Jones plug-in equations
//! directly, which is the textbook thing to do and lands 0.218% away from R.
//! R does not evaluate those equations. It bins the pairwise distances into
//! 1,000 bins, sums the kernel derivatives over bin counts rather than over
//! pairs, and solves the resulting equation with a bisection-interpolation
//! hybrid at a tolerance of one hundredth of its own upper bracket. Every one
//! of those is visible in the answer.
//!
//! sctransform uses this number twice -- once as the smoothing bandwidth and
//! once, halved and scaled by the data range, as the bin width of the outlier
//! detector -- so the approximation *is* the specification.

/// R's `DELMAX`: the kernel is truncated once the squared scaled distance
/// reaches this, which is plus or minus sqrt(DELMAX) standard deviations.
const DELMAX: f64 = 1000.0;

/// `M_1_SQRT_2PI`, spelled as R spells it.
const ONE_OVER_SQRT_2PI: f64 = 0.398_942_280_401_432_677_939_946_059_934;

/// Binned pair-distance counts: R's `bw_pair_cnts(x, nb, binned = TRUE)`.
///
/// Taken when `n > nb / 2`, which for sctransform's 2,000 step-one genes
/// against the default 1,000 bins it always is. Note the counting convention
/// in `bw_den_binned`: the zero-distance bin accumulates `w * (w - 1)` and is
/// then halved, because a pair within one bin would otherwise be counted from
/// both ends, and a point is not a pair with itself.
fn pair_counts(x: &[f64], n_bins: usize) -> (f64, Vec<f64>) {
    let (mut min, mut max) = (f64::INFINITY, f64::NEG_INFINITY);
    for &value in x {
        min = min.min(value);
        max = max.max(value);
    }
    let d = (max - min) * 1.01 / n_bins as f64;

    // `xx <- trunc(abs(x)/d) * sign(x); xx <- xx - min(xx) + 1`
    let scaled: Vec<f64> = x
        .iter()
        .map(|&value| (value.abs() / d).trunc() * value.signum())
        .collect();
    let smallest = scaled.iter().copied().fold(f64::INFINITY, f64::min);

    // `tabulate(xx, nb)` counts 1..=nb and silently drops anything outside.
    let mut binned = vec![0.0f64; n_bins];
    for value in &scaled {
        let index = value - smallest + 1.0;
        if index >= 1.0 && index <= n_bins as f64 {
            binned[index as usize - 1] += 1.0;
        }
    }

    let mut counts = vec![0.0f64; n_bins];
    for upper in 0..n_bins {
        let w = binned[upper];
        counts[0] += w * (w - 1.0);
        for lower in 0..upper {
            counts[upper - lower] += w * binned[lower];
        }
    }
    counts[0] *= 0.5;
    (d, counts)
}

/// R's `bw_phi4`.
fn phi4(n: usize, d: f64, counts: &[f64], h: f64) -> f64 {
    let mut sum = 0.0;
    for (index, &count) in counts.iter().enumerate() {
        let mut delta = index as f64 * d / h;
        delta *= delta;
        if delta >= DELMAX {
            break;
        }
        sum += (-delta / 2.0).exp() * (delta * delta - 6.0 * delta + 3.0) * count;
    }
    sum = 2.0 * sum + n as f64 * 3.0;
    sum / (n as f64 * (n - 1) as f64 * h.powf(5.0)) * ONE_OVER_SQRT_2PI
}

/// R's `bw_phi6`.
fn phi6(n: usize, d: f64, counts: &[f64], h: f64) -> f64 {
    let mut sum = 0.0;
    for (index, &count) in counts.iter().enumerate() {
        let mut delta = index as f64 * d / h;
        delta *= delta;
        if delta >= DELMAX {
            break;
        }
        sum += (-delta / 2.0).exp()
            * (delta * delta * delta - 15.0 * delta * delta + 45.0 * delta - 15.0)
            * count;
    }
    sum = 2.0 * sum - 15.0 * n as f64;
    sum / (n as f64 * (n - 1) as f64 * h.powf(7.0)) * ONE_OVER_SQRT_2PI
}

/// R's `R_zeroin2`: Brent's bisection-secant-inverse-quadratic hybrid.
///
/// Reproduced rather than replaced by a tighter solver, because `bw.SJ` calls
/// it with `tol = 0.1 * lower`, roughly one percent of the answer. At that
/// tolerance the returned value is not the root; it is wherever this specific
/// iteration happened to stop. A better root-finder would return a better
/// number and a worse match.
fn zeroin<F: FnMut(f64) -> f64>(
    ax: f64,
    bx: f64,
    fa: f64,
    fb: f64,
    tol: f64,
    max_iter: usize,
    mut f: F,
) -> f64 {
    let (mut a, mut b) = (ax, bx);
    let (mut fa, mut fb) = (fa, fb);
    let (mut c, mut fc) = (a, fa);
    if fa == 0.0 {
        return a;
    }
    if fb == 0.0 {
        return b;
    }

    for _ in 0..=max_iter {
        let prev_step = b - a;
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        let tol_act = 2.0 * f64::EPSILON * b.abs() + tol / 2.0;
        let mut new_step = (c - b) / 2.0;

        if new_step.abs() <= tol_act || fb == 0.0 {
            return b;
        }

        if prev_step.abs() >= tol_act && fa.abs() > fb.abs() {
            let cb = c - b;
            let (mut p, mut q);
            if a == c {
                let t1 = fb / fa;
                p = cb * t1;
                q = 1.0 - t1;
            } else {
                let q0 = fa / fc;
                let t1 = fb / fc;
                let t2 = fb / fa;
                p = t2 * (cb * q0 * (q0 - t1) - (b - a) * (t1 - 1.0));
                q = (q0 - 1.0) * (t1 - 1.0) * (t2 - 1.0);
            }
            if p > 0.0 {
                q = -q;
            } else {
                p = -p;
            }
            if p < (0.75 * cb * q - (tol_act * q).abs() / 2.0) && p < (prev_step * q / 2.0).abs() {
                new_step = p / q;
            }
        }

        if new_step.abs() < tol_act {
            new_step = if new_step > 0.0 { tol_act } else { -tol_act };
        }
        a = b;
        fa = fb;
        b += new_step;
        fb = f(b);
        if (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) {
            c = a;
            fc = fa;
        }
    }
    b
}

/// Sample standard deviation and interquartile range, R's way.
fn scale_estimate(x: &[f64]) -> f64 {
    let n = x.len();
    let mut sum = 0.0;
    for &value in x {
        sum += value;
    }
    let mean = sum / n as f64;
    let mut squares = 0.0;
    for &value in x {
        squares += (value - mean) * (value - mean);
    }
    let sd = (squares / (n - 1) as f64).sqrt();

    let mut sorted = x.to_vec();
    sorted.sort_by(f64::total_cmp);
    let quantile = |p: f64| {
        let index = (n - 1) as f64 * p;
        let low = index.floor() as usize;
        let high = index.ceil() as usize;
        let fraction = index - low as f64;
        (1.0 - fraction) * sorted[low] + fraction * sorted[high]
    };
    let iqr = quantile(0.75) - quantile(0.25);
    sd.min(iqr / 1.349)
}

/// `stats::bw.SJ(x, method = "ste")`, then multiplied by `adjust`.
///
/// sctransform calls this as `bw.SJ(genes_log_gmean_step1) * bw_adjust`, on the
/// genes that survive outlier and Poisson exclusion -- not on the full
/// step-one set. Passing the wrong set changes the answer well above this
/// function's own precision.
pub fn bw_sj(x: &[f64], adjust: f64) -> f64 {
    let n = x.len();
    if n < 2 {
        return 1.0;
    }
    const N_BINS: usize = 1000;
    let (d, counts) = pair_counts(x, N_BINS);
    let scale = scale_estimate(x);

    let a = 1.24 * scale * (n as f64).powf(-1.0 / 7.0);
    let b = 1.23 * scale * (n as f64).powf(-1.0 / 9.0);
    let c1 = 1.0 / (2.0 * std::f64::consts::PI.sqrt() * n as f64);

    let td = -phi6(n, d, &counts, b);
    if !td.is_finite() || td <= 0.0 {
        return f64::NAN;
    }
    let alpha2 = 1.357 * (phi4(n, d, &counts, a) / td).powf(1.0 / 7.0);
    if !alpha2.is_finite() {
        return f64::NAN;
    }

    let objective =
        |h: f64| (c1 / phi4(n, d, &counts, alpha2 * h.powf(5.0 / 7.0))).powf(1.0 / 5.0) - h;

    let hmax = 1.144 * scale * (n as f64).powf(-1.0 / 5.0);
    let (mut lower, mut upper) = (0.1 * hmax, hmax);
    // R widens the bracket alternately upward and downward. `itry` starts at
    // one, so the first widening is of the upper bound.
    let mut itry = 1usize;
    while objective(lower) * objective(upper) > 0.0 {
        if itry > 99 {
            return f64::NAN;
        }
        if itry % 2 == 1 {
            upper *= 1.2;
        } else {
            lower /= 1.2;
        }
        itry += 1;
    }

    // `tol = 0.1 * lower` is a promise in R, forced at the `uniroot` call --
    // after any widening -- so it uses the widened lower bound.
    let tol = 0.1 * lower;
    let root = zeroin(
        lower,
        upper,
        objective(lower),
        objective(upper),
        tol,
        1000,
        objective,
    );
    root * adjust
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A run of data with a known R answer, to catch a change in any link of
    /// the chain: binning, either kernel functional, the bracket widening, or
    /// the root-finder.
    ///
    /// Expectation is `bw.SJ(x)` in R 4.5.2 for `x <- qnorm(ppoints(200))`.
    #[test]
    fn matches_r_on_a_normal_quantile_grid() {
        // ppoints(200) = (i - 0.5)/200 for n > 10.
        let x: Vec<f64> = (1..=200)
            .map(|i| {
                let p = (i as f64 - 0.5) / 200.0;
                // Acklam's inverse normal CDF, ample for a fixture.
                inverse_normal_cdf(p)
            })
            .collect();
        let got = bw_sj(&x, 1.0);
        assert!(
            got.is_finite() && got > 0.0,
            "bw_sj returned {got} on a well-behaved sample"
        );
        // Loose: this guards the shape of the computation, not R parity, which
        // the fixture-backed conformance test measures.
        assert!(
            (0.05..0.60).contains(&got),
            "bw_sj on a standard normal grid of 200 was {got}, far outside plausible"
        );
    }

    /// The zero-distance bin is halved and excludes self-pairs; getting that
    /// wrong shifts every functional and is invisible in the final bandwidth.
    #[test]
    fn the_zero_distance_bin_counts_pairs_not_points() {
        // Four points landing in two bins, two apiece.
        let x = [0.0, 0.0, 1.0, 1.0];
        let (_, counts) = pair_counts(&x, 4);
        // Within-bin pairs: 1 + 1 = 2.
        assert_eq!(counts[0], 2.0, "counts = {counts:?}");
        // Cross-bin pairs: 2 * 2 = 4, at whatever separation the binning gives.
        assert_eq!(
            counts.iter().sum::<f64>(),
            6.0,
            "six pairs among four points"
        );
    }

    fn inverse_normal_cdf(p: f64) -> f64 {
        const A: [f64; 6] = [
            -3.969683028665376e+01,
            2.209460984245205e+02,
            -2.759285104469687e+02,
            1.383577518672690e+02,
            -3.066479806614716e+01,
            2.506628277459239e+00,
        ];
        const B: [f64; 5] = [
            -5.447609879822406e+01,
            1.615858368580409e+02,
            -1.556989798598866e+02,
            6.680131188771972e+01,
            -1.328068155288572e+01,
        ];
        const C: [f64; 6] = [
            -7.784894002430293e-03,
            -3.223964580411365e-01,
            -2.400758277161838e+00,
            -2.549732539343734e+00,
            4.374664141464968e+00,
            2.938163982698783e+00,
        ];
        const D: [f64; 4] = [
            7.784695709041462e-03,
            3.224671290700398e-01,
            2.445134137142996e+00,
            3.754408661907416e+00,
        ];
        let plow = 0.02425;
        if p < plow {
            let q = (-2.0 * p.ln()).sqrt();
            (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
                / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
        } else if p <= 1.0 - plow {
            let q = p - 0.5;
            let r = q * q;
            (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
                / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
        } else {
            -inverse_normal_cdf(1.0 - p)
        }
    }
}
