// sctransform-rs: GPL-3 Rust port and standalone provider.
// Copyright (C) 2026 ORIC Labs.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, version 3 only.

mod bandwidth;
mod outlier;
mod overdispersion;
mod sctransform;

pub use bandwidth::bw_sj;
pub use outlier::is_outlier;
pub use overdispersion::{
    conventional_loglikelihood, conventional_overdispersion_mle, conventional_score,
    fit_offset_model, objective_at, OffsetFit, OverdispersionEstimate,
};
pub use sctransform::{sctransform, GeneColumns, SctOptions, SctResult};

/// Lanczos approximation used by the baseline negative-binomial likelihood.
///
/// Kept crate-local so the port can replace this independently when matching
/// the upstream R numerical backend.
pub(crate) fn ln_gamma(x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    let coefficients = [
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        1.208650973866179e-3,
        -5.395239384953e-6,
    ];
    let mut y = x;
    let temporary = x + 5.5;
    let temporary = temporary - (x + 0.5) * temporary.ln();
    let mut series = 1.000000000190015;
    for coefficient in coefficients {
        y += 1.0;
        series += coefficient / y;
    }
    -temporary + (2.5066282746310005 * series / x).ln()
}
