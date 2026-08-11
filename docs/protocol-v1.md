# External process protocol v1

Protocol v1 deliberately uses ordinary process arguments and files. It does
not require dynamic linking, shared memory, callbacks, or BioLang-specific
types.

## Request

```text
bl-sctransform-gpl run --input INPUT_DIR --output OUTPUT_DIR [options]
```

`INPUT_DIR` is an uncompressed 10x Matrix Market directory:

- `matrix.mtx`: coordinate matrix, genes by cells;
- `features.tsv`: feature ID in column one, display name in column two when
  available;
- `barcodes.tsv`: one cell barcode per row.

`OUTPUT_DIR` must be absent or empty. This prevents results from different
engines or versions being accidentally combined.

## Response artifacts

- `manifest.csv`: engine identity, version, dimensions, seed, effective
  parameters, timing, and protocol version;
- `genes.csv`: original gene index, regularized theta, intercept, geometric
  mean coordinate, and residual variance;
- `ranking.csv`: genes ordered by decreasing residual variance;
- `fit-genes.csv`: unregularized fit observations;
- `fit-cells.csv`: cells selected for parameter fitting;
- `sampling.csv`: the complete step-one candidate population and its
  unnormalized inverse-density weights;
- `residuals.csv`: validation probe over top-ranked genes and leading cells;
- `matrix.mtx`: optional dense kept residual matrix in genes-by-cells coordinate
  form, written only when `--write-matrix` is requested.

CSV text is UTF-8 with a header. Indices are zero-based. Matrix Market indices
are one-based as required by the format.

## Process contract

- Exit 0: complete, internally consistent output.
- Exit 2: invalid command or arguments.
- Exit 3: invalid input data.
- Exit 4: computation failure.
- Exit 5: output write failure.

Every invocation prints the engine and license identity before computation.
Consumers should capture stdout and stderr in their reproducibility artifact.
