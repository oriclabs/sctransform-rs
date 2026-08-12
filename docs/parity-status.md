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

## Step 2 — the estimator (not started)

Port glmGamPoi's offset-model overdispersion estimator. GPL-3, so permitted in
this repository. This is the substantial one: a C++ IRLS fitter with its own
moment initialisation, Cox-Reid adjustment and convergence policy, all of which
must match, not merely converge to the same neighbourhood.

Target: `unregularized_theta_relative_error_median` below 1e-9, which should
carry `theta_raw_*` and both intercept gates with it.

## Step 3 — regularization (not started)

R's `ksmooth` with the Sheather-Jones bandwidth and `bw_adjust = 3`, over the
od-factor scale. Only measurable once step 2 lands, since the amplification
makes it impossible to attribute error between the two before then.

Residual generation, clipping, residual variance and feature ranking need no
work: they already pass every gate they are measured against.
