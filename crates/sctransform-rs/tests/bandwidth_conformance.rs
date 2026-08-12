//! `bw_sj` against R's `stats::bw.SJ` on the real step-one gene set.
//!
//! sctransform uses this number twice and both uses are sensitive: it is the
//! smoothing bandwidth, and half of it times the data range is the bin width
//! of the outlier detector. A bandwidth that is merely close produces a
//! different outlier set, which produces a different gene set, which produces
//! a different bandwidth.
//!
//! Skipped unless `SCTRANSFORM_REG_FIXTURE` points at a directory written by
//! `validation/export_regularization_exact.R`.

use sctransform_rs::bw_sj;
use std::path::PathBuf;

fn read(path: PathBuf) -> Vec<Vec<String>> {
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("cannot read {}: {error}", path.display()));
    text.lines()
        .skip(1)
        .filter(|line| !line.trim().is_empty())
        .map(|line| line.split(',').map(str::to_string).collect())
        .collect()
}

#[test]
fn bandwidth_matches_r_bw_sj() {
    let Ok(directory) = std::env::var("SCTRANSFORM_REG_FIXTURE") else {
        eprintln!(
            "skipped: set SCTRANSFORM_REG_FIXTURE to a directory written by \
             validation/export_regularization_exact.R"
        );
        return;
    };
    let directory = PathBuf::from(directory);

    let genes = read(directory.join("step1-genes.csv"));
    let all: Vec<f64> = genes
        .iter()
        .map(|row| row[1].parse().expect("log_gmean_step1"))
        .collect();
    let surviving: Vec<f64> = genes
        .iter()
        .filter(|row| row[8] == "1")
        .map(|row| row[1].parse().expect("log_gmean_step1"))
        .collect();

    let mut expected_surviving = f64::NAN;
    for row in read(directory.join("manifest.csv")) {
        if row[0] == "bw_sj" {
            expected_surviving = row[1].parse().expect("bw_sj");
        }
    }
    assert!(
        expected_surviving.is_finite(),
        "fixture manifest has no bw_sj"
    );

    let got_all = bw_sj(&all, 1.0);
    let got_surviving = bw_sj(&surviving, 1.0);
    eprintln!(
        "bw.SJ over {} step-one genes: port {got_all:.17e}",
        all.len()
    );
    eprintln!(
        "bw.SJ over {} surviving genes: port {got_surviving:.17e}, R {expected_surviving:.17e}",
        surviving.len()
    );

    let error = (got_surviving - expected_surviving).abs() / expected_surviving;
    eprintln!("relative error: {error:.3e}");

    // R solves this with `uniroot` at `tol = 0.1 * lower`, so its own answer is
    // only located to about a percent of itself. Reproducing the iteration
    // rather than the root is what makes tight agreement possible at all; the
    // residual is floating-point noise through the two kernel functionals.
    assert!(
        error < 1e-12,
        "bw.SJ relative error {error:e}: port {got_surviving}, R {expected_surviving}"
    );
}
