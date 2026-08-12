// Copyright (C) 2026 ORIC Labs. GPL-3.0-only.
//
// R's density implementation uses the R Core FFT. Building that exact GPL
// implementation avoids platform- and algorithm-dependent sampling drift.
fn main() {
    println!("cargo:rerun-if-changed=native/r_fft.c");
    cc::Build::new()
        .file("native/r_fft.c")
        .warnings(false)
        .compile("r_fft");
}
