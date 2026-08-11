# HBC control validation — 2026-08-12

This is a real-data, same-host comparison on the HBC training control matrix:
14,847 cells, 14,065 input genes, and 13,799 modelled genes. GPU execution was
disabled. The oracle was R 4.5.2 with `sctransform` 0.4.3 and the
`glmGamPoi_offset` method. Each implementation ran as a separate process.

| Measurement | R oracle | BioLang built-in | GPL executable |
|---|---:|---:|---:|
| Transform time | 40.420 s | 2.987 s | 2.954 s |
| Process wall time | 46.763 s | 19.185 s | 5.397 s |
| Peak process-tree working set | 5.526 GiB | 2.292 GiB | 1.644 GiB |
| Step-one candidate genes | 11,075 | not exported | 11,075 |
| Fit-gene overlap with R | 100% | 25.50% | 91.55% |
| Top-3,000 feature overlap with R | 100% | 98.20% | 98.50% |
| Feature-rank Spearman with R | 1.0 | 0.999816 | 0.999867 |
| Residual RMSE / R residual SD | 0% | 1.925% | 1.476% |
| Residual slope | 1.0 | 0.995300 | 0.995809 |
| Median raw theta relative error | 0% | 7.258% | 6.304% |
| P90 raw theta relative error | 0% | 12.523% | 15.819% |

Outcome: the GPL executable is closer on the residuals and feature ranking
that feed downstream PCA, and it uses the least memory. The built-in engine is
closer on theta p90 and intercept calibration. Neither Rust implementation yet
passes every scale-sensitive parameter gate, so this is not a claim of full
SCTransform parity.

The GPL executable's candidate set and density bandwidth match R exactly to
the displayed precision. Its sampling probabilities differ from R by about
floating-point roundoff because this milestone uses direct Gaussian
convolution instead of R's FFT path; sequential sampling amplifies those tiny
differences into a 91.55% fit-gene overlap.

The raw run artifacts are intentionally not committed because they include
large fixture-derived outputs. Use the validation scripts and the comparison
contract in [PORTING.md](../PORTING.md) to reproduce them.
