# Porting and validation policy

## License boundary

The complete workspace is GPL-3.0-only. Source derived from the upstream R/C++
implementation may be translated here, but must never be copied back into the
MIT BioLang repository. BioLang integration is limited to the documented
external-process protocol.

## Provenance

- Upstream: `satijalab/sctransform` 0.4.3, commit
  `49e35b5aeb76a602910207cbfde1561093340be3`, GPL-3.
- Baseline Rust engine: BioLang commit
  `52b72de518aa71a74120b5e8c03eb7cf9daff6bf`, originally MIT licensed.
- Rust workspace and CLI integration: Copyright 2026 ORIC Labs.

The baseline is a control, not an upstream-derived parity claim. Every future
upstream-derived change must identify the relevant upstream function and
commit in its commit message or source comment.

## Acceptance measurements

All engines receive the same count matrix, feature names, cell names,
parameters, and seed. Reports must include:

- regularized and raw parameter comparisons;
- Pearson residual slope, offset, RMSE divided by oracle residual SD, and
  relative-error percentiles;
- residual-variance slope and rank correlation;
- top-feature overlap and rank agreement;
- fit-gene and fit-cell sampling agreement;
- wall-clock transform and process time;
- peak process-tree working set;
- exact engine name, version, input checksum, seed, and parameters.

Correlation alone is never a parity gate. Scale-sensitive slopes and error
percentiles are mandatory.

## Porting order

1. sampling and v2 argument normalization;
2. per-gene parameter fitting;
3. outlier detection and kernel regularization;
4. Pearson residual generation and clipping;
5. residual variance and feature ranking;
6. corrected counts and optional secondary regression;
7. performance and bounded-memory materialization.

