# SCTransform revalidation — 2026-08-12

The standalone GPL Rust provider now passes every scale-sensitive acceptance
gate against fresh runs of the original GPL R packages on a synthetic boundary
fixture and both real HBC matrices.

## Versions and method

- GPL Rust: working tree based on commit `83030c5`, `sctransform-rs` 0.1.0.
- BioLang built-in comparison: commit `52b72de`.
- R 4.5.2, `sctransform` 0.4.3, `glmGamPoi` 1.22.0.
- Seed 1448145; 5,000 fit cells; 2,000 fit genes; `min_cells = 5`.
- Real data: HBC control (14,847 cells) and stimulated (14,782 cells).
- Synthetic boundary fixture: 480 cells and 120 input genes.

The R oracle and GPL Rust executable ran as separate processes. Runs were
fresh, sequential, on the same Windows host. The comparator uses original-scale
slopes and relative-error percentiles as well as correlations. The Rust
workspace passed all 39 tests.

## Final GPL Rust result

| Measurement | HBC control | HBC stimulated | Synthetic |
|---|---:|---:|---:|
| Unregularized theta median relative error | 5.55e-11 | 5.96e-11 | 4.80e-10 |
| Fit-gene overlap | 100.0% | 96.8% | 100.0% |
| Regularized theta median relative error | 7.90e-8 | 0.0858% | 6.45e-8 |
| Regularized theta p90 relative error | 9.89e-8 | 1.837% | 1.02e-7 |
| Regularized theta slope | 1.000000007 | 1.000248 | 1.000000034 |
| Intercept RMSE | 1.44e-5 | 0.001362 | 0.000103 |
| Top-feature overlap | 100.0% | 99.933% | 100.0% |
| Feature-rank Spearman | 1.000000 | 1.000000 | 1.000000 |
| Residual RMSE / R residual SD | 1.56e-8 | 0.1297% | 0.00260% |
| Residual slope | 0.999999995 | 1.000487 | 1.0000078 |
| Transform time, R / Rust | 33.05 / 1.845 s | 33.10 / 1.635 s | 3.67 / 0.0089 s |
| Process wall time, R / Rust | 38.94 / 3.46 s | 38.65 / 3.29 s | 4.64 / 0.15 s |
| Peak working set, R / Rust | 5.527 / 1.672 GiB | 5.453 / 1.654 GiB | 0.697 / 0.0057 GiB |

All 23 applicable gates pass on each fixture. On the real matrices the
standalone provider used about 30% of R's peak memory and completed the
measured process 11.2–11.7 times faster.

## What prevented portable bit identity

The remaining stimulated-data difference is isolated to density-weighted gene
selection, not the random generator or sampling algorithm:

1. The provider now uses R Core's GPL FFT, R's two interpolation passes, its
   probability normalisation, descending heapsort, seeded Mersenne Twister, and
   without-replacement sampler.
2. The 5,000 sampled cells match R in exact order.
3. Giving R the provider-exported weights reproduces the provider's 2,000
   sampled genes in exact order.
4. Last-bit density differences from compiler and math-library arithmetic
   reorder near-tied weights before heapsort. On stimulated HBC this changes 64
   of 2,000 fit genes.

This is why portable, cross-build bit identity cannot be promised. It is not a
remaining algorithmic omission: common-gene unregularized theta is effectively
exact, top-feature overlap is 99.933%, and residual error is 0.1297% of the R
residual standard deviation. Control and synthetic sampling are exact.

The former `theta_bias_attenuated_below_0_25` acceptance check was removed.
It divided residual error by theta error and therefore became unstable when
theta error approached floating-point noise. The attenuation ratio remains a
diagnostic metric; absolute theta slope/error and per-gene residual error are
the acceptance criteria.

## BioLang built-in comparison

The MIT clean-room built-in remains useful when a GPL component is not wanted,
but the optional GPL provider is closer to the R implementation on both real
datasets:

| Measurement | Engine | HBC control | HBC stimulated | Synthetic |
|---|---|---:|---:|---:|
| Regularized theta median error | built-in | 7.258% | 3.914% | 3.356% |
|  | GPL Rust | 7.90e-8 | 0.0858% | 6.45e-8 |
| Top-feature overlap | built-in | 98.200% | 97.767% | 100.0% |
|  | GPL Rust | 100.0% | 99.933% | 100.0% |
| Residual RMSE / R SD | built-in | 1.925% | 2.926% | 0.204% |
|  | GPL Rust | 1.56e-8 | 0.1297% | 0.00260% |
| Transform time | built-in | 2.757 s | 2.835 s | 0.017 s |
|  | GPL Rust | 1.845 s | 1.635 s | 0.0089 s |
| Peak working set | built-in | 2.292 GiB | 2.276 GiB | 0.048 GiB |
|  | GPL Rust | 1.672 GiB | 1.654 GiB | 0.0057 GiB |

Raw artifacts are under ignored local `validation-results/final-pass-*`
directories. Each contains the provider outputs, resource measurement, and
comparison JSON; the independently generated R artifacts remain under
`revalidation-2026-08-12*`.
