# Parity status

Where the port stands against R `sctransform` 0.4.3 with `glmGamPoi` 1.22.0,
and what is left. Measurements are on the HBC control matrix: 14,847 cells,
14,065 input genes, 13,799 modelled genes.

## Reachable target

Bit-identical output is not reachable and is not the goal. R sums through its
own BLAS in its own order, so the last bits will differ however faithful the
port is. The target is agreement to roughly 1e-10 relative, at which every
acceptance gate passes and the ranked feature lists are identical gene for
gene. No downstream analysis can distinguish two results that close.

## Step 1 — step-one gene sampling (done)

**Fit-gene overlap with R: 91.55% → 100.00%.** Bandwidth is now bit-identical
to `stats::density(..., bw = "nrd")`: `0.1113703131988174`.

The cause was not where the code comments predicted. The port's density
estimate uses direct Gaussian convolution where R uses linear binning and an
FFT, and that was assumed to be the source of the divergence. Measured at full
precision, it is not: the sampling weights already agreed with R to a median
relative difference of 1.6e-15, and the candidate population was already
identical at 11,075 genes.

The first link that actually differed was one step further upstream. The port
computed each gene's geometric mean with `expm1(m)`; upstream
`row_gmean_dgcmatrix` computes `exp(m) - eps`. `expm1` is the more accurate of
the two — which is exactly why it was wrong here. Only 8.32% of genes matched
R's log geometric mean. Switching to the upstream form took that to 99.73%,
and because that value is the abscissa of the whole density estimate, the
bandwidth, the interpolated weights and the sample all fell into place behind
it: 100% of R's 2,000 fit genes, from a single arithmetic form.

The lesson generalises to the remaining steps. A sequential sample *without
replacement* has no error budget: one flipped draw changes the population every
later draw sees. So divergence must be chased to the first link in the chain,
not the loudest one, and "more accurate than the reference" is a defect.

`validation/export_density_exact.R` is the probe that located it. It prints
every intermediate at `%.17g` so a comparison measures the disagreement rather
than the printer's 15 digits.

Guarded by `geometric_mean_uses_the_upstream_inaccurate_form`, which was
mutation-tested: reverting to `expm1` fails that test and no other.

## What step 1 did not fix

Downstream accuracy is essentially unchanged, and that is the informative part.

| Measurement | before | after |
|---|---:|---:|
| Fit-gene overlap | 91.550% | **100.000%** |
| Unregularized theta rel. error, median | 1.258% | 1.323% |
| Unregularized theta rel. error, p90 | 6.925% | 7.120% |
| Unregularized theta regression slope | 1.022193 | 1.022154 |
| Median regularized theta rel. error | 6.304% | 6.489% |
| Intercept RMSE | 0.138172 | 0.136828 |
| Top-3,000 feature overlap | 98.500% | 98.433% |
| Residual RMSE / R residual SD | 1.476% | 1.453% |

Gene selection was not the bottleneck. What step 1 bought, besides the overlap
itself, is a clean measurement: raw theta was previously compared on the 1,650
genes both implementations happened to choose, a self-selected subset. It is
now compared on all 2,000 of the same genes, and the estimator disagreement
stands at 1.32% median with a systematic slope of 1.022.

That residual bias is the entire remaining story. The port fits its own
Cox-Reid negative binomial; R runs glmGamPoi's C++ `glm_gp`. A 1.3% difference
in raw theta becomes 6.5% after regularization because `d ln(theta)/d od` is
about -25 at these UMI values, so the od-factor inversion amplifies it roughly
fivefold.

## Step 2 — the estimator (done, wired in)

**Median relative theta error against glmGamPoi: 1.3e-2 → 1.3e-10.** Measured
on 300 genes sampled evenly across the expression range, handed the exact
`(y, mu)` glmGamPoi optimised over.

| | value |
|---|---:|
| Genes with measurable overdispersion | 260 |
| Genes both call unoverdispersed | 40 |
| Median relative theta error | 1.255e-10 |
| P90 | 2.367e-7 |
| Max | 1.133e-5 |
| Worst relative objective gap at the two estimates | 1.059e-11 |

That last row is the one that carries the claim. Evaluated at either
implementation's answer, the objective returns the same value to eleven
figures, which a mistranscribed term could not do. What remains is not an
estimator difference but a stopping rule: `nlminb` gets an analytic gradient
and halts where that gradient is numerically zero.

Two findings worth keeping.

**The reference stops on the gradient, not the objective.** Polishing the
score's root by golden-section maximisation of the loglikelihood is the obvious
correction -- `nlminb` nominally minimises the objective -- and it made
agreement four thousand times *worse*, median 1.3e-10 to 7.5e-7. Upstream's
score carries clamps and Taylor brackets that make it a deliberately inexact
derivative of its own loglikelihood, so the two have slightly different
stationary points, and the reference sits at the gradient's. The polish was
removed.

**glm_gp's returned `Beta` is not the one the estimator used.** It fits beta,
estimates overdispersion from that `Mu`, shrinks the dispersions, then fits
beta again. Reconstructing mu from the returned `Beta` scores a port against
inputs the reference never optimised over.
`validation/export_overdispersion_fixture.R` walks the stages explicitly and
asserts the walk reproduces `glm_gp`'s own output.

### Wired into the pipeline

The beta stage (`estimate_dispersions_by_moment`,
`estimate_betas_roughly_group_wise`, `fitBeta_one_group`) is ported too, and
`sctransform()` now fits through this chain rather than the crate's own
Cox-Reid fit, which has been deleted.

| Measurement | baseline | after step 1 | after step 2 |
|---|---:|---:|---:|
| Fit-gene overlap | 91.550% | 100.000% | 100.000% |
| **Unregularized theta rel. error, median** | **1.258%** | 1.323% | **5.6e-9%** |
| Unregularized theta rel. error, p90 | 6.925% | 7.120% | 9.1e-6% |
| Unregularized theta, port / oracle | — | — | 1.000000000046 |
| Median *regularized* theta rel. error | 6.304% | 6.489% | 6.709% |
| od-factor median absolute difference | 0.002784 | 0.002962 | 0.002980 |
| Intercept slope | 0.978152 | 0.978716 | 0.981348 |
| Top-3,000 feature overlap | 98.500% | 98.433% | 98.733% |
| Feature-rank Spearman | 0.999867 | 0.999870 | 0.999886 |
| Transform time | 2.95 s | 2.78 s | 2.04 s |

One gate flipped to passing (`intercept_slope_0_98_to_1_02`), leaving four.

## What this proves about the rest

The raw estimator now agrees with glmGamPoi to eleven figures **and the
regularized theta did not move**. Those two facts together are the point of
this step. The remaining 6.7% is not an estimator difference, not a gene-sample
difference, and not a beta difference -- all three of those are now measured at
or near machine precision. It is entirely in the regularization stage, and the
od-factor difference of 0.00298 is where it lives.

That difference is small in its own units and large after inversion:
`d ln(theta) / d od` is about -25 at these UMI values, so 0.003 in od becomes
roughly 7% in theta. Nothing downstream of the regularizer can be improved
without fixing the regularizer first.

## Step 3 — regularization (started; carries all remaining error)

`reg_model_pars` is a longer chain than "smooth the od-factor", and each link
is a candidate:

1. `dispersion_par = log10(1 + 10^log_gmean_step1 / theta)`;
2. `is_outlier` on **every** column of `model_pars`, scored by
   `robust_scale_binned` against two bin grids offset by half a bin width, with
   a threshold of 10 on the smaller absolute score -- outliers are dropped from
   the fit;
3. `exclude_poisson` additionally drops genes with `theta = Inf`, with
   `variance <= mean`, and with `mean < 0.001`, and replaces them with an
   analytic offset model rather than a smoothed one;
4. `bw <- bw.SJ(genes_log_gmean_step1) * bw_adjust`;
5. `x_points` clamped to the step-1 range, then
   `ksmooth(kernel = "normal", bandwidth = bw)`.

Point 5 hides two constants in R's C source: `ksmooth` rescales the bandwidth
by `0.3706506` and **truncates the kernel at four scaled bandwidths**. A
Gaussian smoother without that truncation is a different estimator.

### First measurement

R's `bw.SJ` on the 2,000 step-one genes is `0.184475430376767`. This crate's
`sheather_jones_bandwidth` returns `0.18407329286868537` -- **0.218% low**.

That is a real divergence but probably not the whole 0.00298 in od-factor, so
it is a starting point rather than an answer. R's `bw.SJ` is not the 1991
equations evaluated directly, which is what this crate implements: it bins the
pair distances into 1,000 bins and solves with `uniroot`, and the binning is
part of the result rather than an optimisation of it.

The links above are listed in the order they must be checked, because the
lesson from steps 1 and 2 is that only the *first* divergence is worth acting
on. In particular the outlier and Poisson-exclusion rules change *which genes
are smoothed*, and no bandwidth work means anything if the fit set is wrong.

### The first divergence was the fit set, not the bandwidth

Exporting every link at `%.17g` settled it. `dispersion_par` itself agrees to
5.1e-11, and the step-one gene sets are identical. What differed is who is
smoothed: **R smooths 1,712 genes, this port smoothed 1,780**, and all 68
extras are genes R excludes — 18 of R's outliers, 29 of R's Poisson genes, and
21 where R's estimator returns an infinite theta and this one does not.

### Poisson exclusion (fixed)

Upstream applies *two different Poisson rules at two different stages*, and
this port applied the narrower one to both. Step-one sampling uses
`variance - mean > 0` alone. Regularization builds `all_poisson_genes` as the
union of that and `mean < 0.001`, and fits those with an analytic offset model
rather than letting them into the smoother.

| Measurement | step 2 | + Poisson rule |
|---|---:|---:|
| **Intercept RMSE** | 0.135231 | **0.031209** |
| **Intercept slope** | 0.981348 | **0.998750** |
| Median regularized theta rel. error | 6.709% | 6.465% |
| od-factor median absolute difference | 0.002980 | 0.002833 |
| Residual-variance slope | 1.005547 | 0.992416 |
| Residual RMSE / R residual SD | 1.599% | 1.628% |

`intercept_rmse_at_most_0_10` now passes. The intercept is essentially done.

The residual numbers moved the wrong way by a little, and the reason is worth
stating rather than hiding: an intercept error was partly cancelling the theta
error in the residuals, and removing the intercept error removed the
cancellation. The same arithmetic moved `theta_to_residual_attenuation` from
0.229 to 0.258 and so failed a gate that had been passing -- theta error on the
probe genes fell from 4.43% to 3.96% while residual error stayed at 1.02%, so
the *ratio* rose. Nothing got worse except a ratio with a smaller denominator.

All four remaining gate failures are now theta. There is no longer an intercept
problem to hide behind.

### Still to do in step 3

- `is_outlier`: two bin grids offset by half a bin width, `robust_scale` per
  bin as `(y - median) / (mad + eps)` with `mad`'s 1.4826 constant, threshold
  10 on the smaller absolute score. 18 genes.
- The 21 genes where R's estimator returns infinite theta and this one returns
  a finite value -- a boundary-condition difference in the MLE's early exits,
  not in the objective.
- `bw.SJ` proper: R bins pair distances into 1,000 bins and solves with
  `uniroot`; this crate evaluates the 1991 equations directly and lands 0.218%
  low. Needed twice, since `is_outlier` derives its bin width from it.
- `ksmooth`'s `0.3706506` rescaling and four-bandwidth truncation.

### Licensing, which changed under us

glmGamPoi relicensed GPL-3 to MIT on 26 May 2026 (`a9eeed642`), and the four
C++ files are byte-identical across that change. That would put the estimator
within reach of MIT BioLang -- except that `src/overdispersion.cpp` still
carries an in-file notice marking it LGPL (>= 3) and attributing it to DESeq2,
and the relicensing commit did not touch it. This port treats that file as
LGPL-3 and conveys it under GPL-3, which LGPL-3 section 2 permits.
`beta_estimation.cpp` carries no such notice and is MIT at HEAD.

Residual generation, clipping, residual variance and feature ranking need no
work: they already pass every gate they are measured against.
