# Independent reproduction — 2026-08-12

This record is a from-scratch re-measurement of the claims in
[validation-hbc-control.md](validation-hbc-control.md), performed without
reusing any previously written artifact. The fixture was rebuilt from the raw
10x directories, the R oracle was re-run, both Rust engines were re-run, and the
headline metrics were additionally recomputed by a second, independent script
reading the raw CSVs directly.

## What was executed

| step | command |
|---|---|
| fixture | `prepare_hbc_sctransform_fixture.R ctrl_raw stim_raw <fixture>` |
| oracle | `sctransform_oracle.R tenx <fixture>/ctrl <oracle>` |
| GPL engine | `bl-sctransform-gpl run --input <fixture>/ctrl --output <gpl> --probe-genes 3000 --probe-cells 64` |
| built-in | `bl run sctransform_biolang.bl` with `BIOLANG_GPU=0` |
| comparison | `compare_sctransform_results.py <oracle> <engine> <json> 3000` |

The harness lives in the MIT `biolang-workflows` repository under
`validation/single-cell/`. It is not vendored here; see
[Reproducibility gap](#reproducibility-gap).

## Environment

R 4.5.2, `sctransform` 0.4.3, `glmGamPoi` 1.22.0, Matrix 1.7.4. The oracle
manifest confirms the effective arguments were `method = glmGamPoi_offset`,
`exclude_poisson = TRUE`, `theta_regularization = od_factor`,
`min_variance = umi_median`, `bw_adjust = 3`, `n_cells = 5000`,
`n_genes = 2000`, `min_cells = 5`, seed 1448145. GPU execution disabled.

The fixture script's built-in checkpoint passed exactly: 14,847 control cells,
14,782 stimulated cells, 14,065 retained genes. All three engines independently
reported 14,847 cells, 14,065 input genes, and 13,799 modelled genes, and the
modelled-gene sets are identical (Jaccard 1.0).

## Accuracy against R `sctransform` 0.4.3

Both engines measured against the same oracle run.

| Measurement | built-in | GPL executable | claimed built-in | claimed GPL |
|---|---:|---:|---:|---:|
| Median raw theta relative error | 7.258% | 6.304% | 7.258% | 6.304% |
| P90 raw theta relative error | 12.523% | 15.819% | 12.523% | 15.819% |
| Intercept RMSE | 0.064979 | 0.138172 | 0.064979 | 0.138172 |
| Residual-variance slope | 1.026886 | 1.016104 | 1.026886 | 1.016104 |
| Top-3,000 feature overlap | 98.200% | 98.500% | 98.20% | 98.50% |
| Feature-rank Spearman | 0.999816 | 0.999867 | 0.999816 | 0.999867 |
| Residual RMSE / R residual SD | 1.925% | 1.476% | 1.925% | 1.476% |
| Residual slope | 0.995300 | 0.995809 | 0.995300 | 0.995809 |
| Fit-gene overlap | 25.500% | 91.550% | 25.50% | 91.55% |
| Fit-cell overlap | 100.000% | 100.000% | — | — |
| od-factor median absolute difference | 0.002912 | 0.002784 | — | — |

Every previously published accuracy figure reproduced to the digits it was
printed at.

### Second-source check

The four load-bearing numbers were recomputed by a separate script that reads
`ranking.csv`, `residuals.csv`, `genes.csv` and `fit-genes.csv` and does its own
set and moment arithmetic, sharing no code with the comparison harness:

| Measurement | harness | independent |
|---|---:|---:|
| Top-3,000 feature overlap | 98.500% | 2,955 / 3,000 = 98.5000% |
| Residual RMSE / oracle SD | 1.47639% | 1.4764% |
| Median theta relative error | 6.30361% | 6.3036% |
| Fit-gene overlap | 91.550% | 1,831 / 2,000 = 91.5500% |

The residual comparison covers 189,120 of 192,000 probe observations. The
2,880-cell shortfall is exactly the 45 genes by which the two top-3,000 lists
differ, times 64 probe cells — an internally consistent accounting, not
silently dropped data.

## Acceptance gates

Of 21 gates, the GPL executable fails 5 and the built-in fails 4.

| gate | built-in | GPL executable |
|---|---|---|
| `theta_raw_slope_0_98_to_1_02` | fail | fail |
| `theta_raw_median_relative_error_at_most_0_05` | fail | fail |
| `theta_raw_p90_relative_error_at_most_0_10` | fail | fail |
| `variance_slope_0_98_to_1_02` | fail | pass |
| `intercept_slope_0_98_to_1_02` | pass | fail |
| `intercept_rmse_at_most_0_10` | pass | fail |
| all 16 remaining gates | pass | pass |

Both engines pass every residual gate — correlation, slope, RMSE-over-SD, and
large-residual relative error — and every feature-selection gate. The failures
are concentrated in the regularized theta and intercept parameters, which is
the expected signature of the od-factor conditioning: the median absolute
od-factor difference is 0.0028, and `d ln(theta) / d od` is approximately -25 at
these UMI values, which turns that into the ~6% theta spread observed.

The port's own gates therefore behave as designed: they fail on the parameters
that genuinely differ and pass on the quantities that feed PCA.

## Resources

Same host, sequential runs, peak working set sampled from the OS counter over
the whole process tree.

| | R oracle | built-in | GPL executable |
|---|---:|---:|---:|
| Transform time | 76.86 / 79.61 / 80.19 s | 3.223 s | 2.778 / 2.867 s |
| Process wall time | 84.6 / 87.0 s | 25.03 s | 4.10 / 4.23 s |
| Peak working set | 5.518 GiB | 2.285 GiB | 1.636 GiB |

Memory reproduced closely: 5.518 vs 5.526 GiB claimed for R, 1.636 vs 1.644 GiB
for the GPL executable, 2.285 vs 2.292 GiB for the built-in.

### One figure did not reproduce

`validation-hbc-control.md` records the R oracle transform at **40.420 s**.
Three runs on this host measured **76.86 s, 79.61 s, and 80.19 s** — roughly
twice as slow, with no run anywhere near 40 s. The Rust timings from the same
sessions reproduced normally, so this is specific to the R arm rather than a
generally slower host.

The published figure is conservative in the port's disfavour: at the measured
oracle time the GPL executable is about 29x faster on the transform, not the
13x the current document implies. The R number should be re-measured and the
document corrected before it is cited.

## Reproducibility gap

`validation-results/` is `.gitignore`d, so no artifact in this repository backs
the published table, and the comparison harness lives in a different repository.
A checkout of `sctransform-rs` alone cannot reproduce or refute its own
validation claims. Two options, in preference order:

1. commit the small comparison JSONs (about 5 KB each) as evidence, and add a
   `validation/run.ps1` that drives the external harness by path; or
2. vendor the oracle and comparison scripts here. They are MIT, so relicensing
   into this GPL-3 tree is permitted, at the cost of a second copy to maintain.

CI cannot run any of this — it needs R, glmGamPoi and a multi-gigabyte fixture —
so committed evidence is the only mechanism that would catch a regression.
