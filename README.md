# sctransform-rs

A GPL-3.0-only Rust library and standalone executable for SCTransform-compatible
normalization of single-cell UMI count matrices.

This repository is intentionally separate from BioLang. `bl.exe` and BioLang's
built-in, paper-derived implementation remain MIT licensed. BioLang may invoke
the executable in this repository as an external process; it must not link this
library into the MIT executable.

## Status

The first milestone now includes GPL-derived compatibility translations for
upstream step-one filtering, R's `density(..., bw = "nrd")` path, R-compatible
weighted sampling, and seeded cell/gene selection. It is **not yet claimed to
be a complete line-by-line Rust port of upstream `sctransform`**, and it does
not yet pass every scale-sensitive theta/intercept parity gate. Changes are
accepted only after three-way comparison with:

1. the original R `sctransform` package,
2. BioLang's MIT built-in engine, and
3. this standalone executable.

See [PORTING.md](PORTING.md) for the provenance and parity policy.
The latest three-dataset comparison is recorded in
[docs/revalidation-2026-08-12.md](docs/revalidation-2026-08-12.md). The earlier
baseline is retained in
[docs/validation-hbc-control.md](docs/validation-hbc-control.md).

## Build

```powershell
cargo build --release
```

The executable is `target/release/bl-sctransform-gpl.exe` on Windows and
`target/release/bl-sctransform-gpl` on Unix.

## Run

Input is an uncompressed 10x Matrix Market directory containing `matrix.mtx`,
`features.tsv`, and `barcodes.tsv`:

```powershell
target/release/bl-sctransform-gpl.exe run `
  --input path/to/filtered_feature_bc_matrix `
  --output validation-results/run-1
```

The output directory must be absent or empty. The command writes model
parameters, ranked features, residual probes, sampling diagnostics, and a
manifest in a documented, versioned format. Add `--write-matrix` when the full
dense kept residual matrix is required.

```powershell
target/release/bl-sctransform-gpl.exe license
target/release/bl-sctransform-gpl.exe --version
```

## BioLang boundary

```text
bl.exe (MIT)
    |
    | files / command line / process exit status
    v
bl-sctransform-gpl (GPL-3.0-only)
```

The process protocol is documented in [docs/protocol-v1.md](docs/protocol-v1.md).
The GPL component's backend identity and version are written into every run
manifest so methods sections cannot silently confuse it with BioLang's native
engine.

## License

Copyright and provenance are listed in [NOTICE.md](NOTICE.md). This repository
is distributed under GNU GPL version 3 only. See [LICENSE](LICENSE).
