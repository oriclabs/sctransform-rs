# Copyright and provenance notices

This work is distributed under GNU GPL version 3 only.

The SCTransform method and upstream implementation are authored by Christoph
Hafemeister, Saket Choudhary, Rahul Satija, and contributors to
`satijalab/sctransform`. The upstream R package is licensed GPL-3.

The Cox-Reid adjusted overdispersion estimator is derived from glmGamPoi
(<https://github.com/const-ae/glmGamPoi>), Copyright Constantin Ahlmann-Eltze.
That package's `src/overdispersion.cpp` states that its likelihood, score and
optimisation routines were originally taken from DESeq2's `src/DESeq2.cpp` by
Michael I. Love, and marks them "License: LGPL (>= 3)".

The glmGamPoi package relicensed from GPL-3 to MIT on 26 May 2026, in commit
`a9eeed642`. That commit changed only `DESCRIPTION`, `LICENSE`, `LICENSE.md`
and `.Rbuildignore`; the in-file LGPL notice on `src/overdispersion.cpp` and
its attribution to DESeq2 were left in place, and the file's contents are
byte-identical before and after. This port therefore treats that file as
LGPL-3, not MIT. GNU LGPL version 3 section 2 permits conveying such a work
under the GNU GPL version 3, which is what this repository does.

The density binning, FFT, interpolation, and unequal-probability sampling
compatibility work is derived from the R statistical computing environment,
Copyright the R Core Team, the R Foundation, Robert Gentleman, Ross Ihaka,
Adrian Baddeley, and other contributors. The vendored FFT routines are R
Core's C translation of Richard Singleton's mixed-radix FFT. The relevant R
source is offered under GPL version 2 or later and is incorporated here under
GPL version 3. The vendored file was taken from R source commit
`a2066dd40b0c7ee16c388c6153f1f884faa50b24`.

The initial Rust normalization engine was taken from the MIT-licensed BioLang
project at commit `52b72de518aa71a74120b5e8c03eb7cf9daff6bf`:

> Copyright (c) 2024 ORIC Labs (oriclabs)

The MIT permission and warranty notice applicable to that original contribution
is preserved in [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md). The combined
workspace is conveyed under GPL-3.0-only.

Rust workspace, CLI, protocol, and subsequent port modifications:

> Copyright (C) 2026 ORIC Labs

This is an independent community port and is not endorsed by or affiliated
with the Satija Lab. “SCTransform” and project names are used only to describe
compatibility and provenance.
