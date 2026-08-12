//! The overdispersion estimator, scored against glmGamPoi on real inputs.
//!
//! This is the test the whole step exists for. Every other measurement of
//! theta in this repository is taken at the end of the pipeline, where a
//! disagreement could equally be the gene sample, the beta stage, or the
//! regularisation. Here the estimator is handed the exact `(y, mu)` glmGamPoi
//! optimised over and asked for the same answer, so a failure has one possible
//! cause.
//!
//! The fixture is large and derived from data that cannot be committed, so the
//! test is skipped unless `SCTRANSFORM_OD_FIXTURE` points at a directory
//! produced by `validation/export_overdispersion_fixture.R`. Skipping silently
//! is deliberate: this must not fail a CI checkout that legitimately has no R.
//! It prints what it did either way.

use sctransform_rs::{conventional_overdispersion_mle, objective_at};
use std::path::PathBuf;

/// Above this, the likelihood is flat enough that the estimate is not a
/// measurement of anything. A gene with theta of a million has no detectable
/// overdispersion; the regulariser turns it into an od-factor of order 1e-7 and
/// the difference between two such estimates cannot reach the residuals. Both
/// implementations must agree that a gene is in this regime, but not on where
/// in it they landed.
const FLAT_THETA: f64 = 1e4;

fn rows(path: &PathBuf) -> Vec<Vec<String>> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|error| panic!("cannot read {}: {error}", path.display()));
    text.lines()
        .skip(1)
        .filter(|line| !line.trim().is_empty())
        .map(|line| line.split(',').map(str::to_string).collect())
        .collect()
}

#[test]
fn overdispersion_matches_glmgampoi() {
    let Ok(directory) = std::env::var("SCTRANSFORM_OD_FIXTURE") else {
        eprintln!(
            "skipped: set SCTRANSFORM_OD_FIXTURE to a directory written by \
             validation/export_overdispersion_fixture.R"
        );
        return;
    };
    let directory = PathBuf::from(directory);

    let log_umi: Vec<f64> = rows(&directory.join("cells.csv"))
        .iter()
        .map(|row| row[3].parse().expect("log_umi"))
        .collect();
    let n_cells = log_umi.len();

    let genes = rows(&directory.join("genes.csv"));
    let mut counts: Vec<Vec<f64>> = vec![vec![0.0; n_cells]; genes.len()];
    for row in rows(&directory.join("counts.csv")) {
        let gene: usize = row[0].parse().expect("gene_row");
        let cell: usize = row[1].parse().expect("cell_col");
        counts[gene][cell] = row[2].parse().expect("count");
    }

    let mut relative: Vec<f64> = Vec::new();
    let mut per_gene: Vec<(f64, String)> = Vec::new();
    let mut agreed_no_overdispersion = 0usize;
    let mut regime_disagreements: Vec<String> = Vec::new();
    let mut objective_gaps: Vec<f64> = Vec::new();

    for (index, row) in genes.iter().enumerate() {
        let name = &row[0];
        let expected_theta: f64 = if row[2] == "Inf" {
            f64::INFINITY
        } else {
            row[2].parse().expect("theta")
        };
        // mu is dense and fully determined: exp(intercept + offset). The
        // intercept is the first-stage one, which is what the estimator saw.
        let intercept: f64 = row[4].parse().expect("mle_intercept");
        let mu: Vec<f64> = log_umi
            .iter()
            .map(|offset| (intercept + offset).exp())
            .collect();

        let estimate = conventional_overdispersion_mle(&counts[index], &mu, true);
        let got = estimate.theta();

        // Both must place the gene in the same regime. That is the part a
        // wrong objective would break; where inside the flat regime each
        // stopped is not.
        let (ours_flat, theirs_flat) = (got >= FLAT_THETA, expected_theta >= FLAT_THETA);
        if ours_flat != theirs_flat {
            regime_disagreements.push(format!(
                "{name}: theta {got:.6e} vs glmGamPoi {expected_theta:.6e} ({})",
                estimate.message
            ));
            continue;
        }
        if ours_flat {
            agreed_no_overdispersion += 1;
            continue;
        }

        let error = (got - expected_theta).abs() / expected_theta.abs();
        relative.push(error);
        per_gene.push((
            error,
            format!("{name}: theta {got:.10e} vs glmGamPoi {expected_theta:.10e}"),
        ));

        // Both answers, on the same function. If the port had transcribed the
        // objective wrongly, its optimum would sit somewhere this function does
        // not agree is an optimum, and the gap would be structural rather than
        // the width of a rounding error.
        let ours = objective_at(&counts[index], &mu, estimate.overdispersion.ln(), true);
        let theirs = objective_at(&counts[index], &mu, (1.0 / expected_theta).ln(), true);
        objective_gaps.push((ours - theirs).abs() / theirs.abs().max(1.0));
    }

    assert!(
        regime_disagreements.is_empty(),
        "{} genes were placed in different overdispersion regimes:\n  {}",
        regime_disagreements.len(),
        regime_disagreements.join("\n  ")
    );

    relative.sort_by(f64::total_cmp);
    objective_gaps.sort_by(f64::total_cmp);
    let count = relative.len();
    assert!(count > 0, "fixture contained no measurable-theta genes");
    let at =
        |values: &[f64], quantile: f64| values[(quantile * (count - 1) as f64).round() as usize];
    let (median, p90, worst) = (at(&relative, 0.5), at(&relative, 0.9), relative[count - 1]);
    let objective_worst = objective_gaps[objective_gaps.len() - 1];

    eprintln!(
        "theta vs glmGamPoi over {count} measurable genes \
         ({agreed_no_overdispersion} agreed to have none): \
         median {median:.3e}, p90 {p90:.3e}, max {worst:.3e}"
    );
    eprintln!("worst relative objective gap at the two estimates: {objective_worst:.3e}");

    // The substantive claim: the two implementations are optimising the same
    // function. Evaluated at either answer it returns the same value to within
    // rounding, which a mistranscribed term could not do.
    assert!(
        objective_worst < 1e-9,
        "objective differs by {objective_worst:e} between the two estimates, \
         so it is not the same function"
    );

    per_gene.sort_by(|a, b| b.0.total_cmp(&a.0));
    let offenders: Vec<String> = per_gene
        .iter()
        .take(8)
        .map(|(error, description)| format!("{error:.3e}  {description}"))
        .collect();

    // glmGamPoi optimises with nlminb, which is handed an analytic gradient and
    // therefore stops where that gradient is numerically zero. Agreement cannot
    // be asked to beat the reference's own stopping rule, so these gates assert
    // the disagreement is that slack and not an estimator difference -- the
    // starting point for this work was a median of 1.3e-2.
    assert!(
        median < 1e-8,
        "median relative theta error was {median:e}\nworst genes:\n  {}",
        offenders.join("\n  ")
    );
    assert!(
        worst < 1e-4,
        "worst relative theta error was {worst:e}\nworst genes:\n  {}",
        offenders.join("\n  ")
    );
}

/// The beta stage, against the intercept glmGamPoi actually optimised from.
///
/// This is separated from the estimator test above because they fail for
/// different reasons. The estimator is scored on inputs taken from R; this is
/// scored on inputs the port derives itself, so it is the one that catches a
/// wrong moment estimate or a Newton loop that stops in a different place.
///
/// It matters more than its size suggests. A gene with no real overdispersion
/// sits where the score at the far-left probe is a cancellation of two
/// quantities each around `sum(y)`, so a small change in mu flips the sign and
/// with it the difference between a finite theta and an infinite one.
#[test]
fn the_beta_stage_matches_glmgampoi() {
    let Ok(directory) = std::env::var("SCTRANSFORM_OD_FIXTURE") else {
        eprintln!("skipped: SCTRANSFORM_OD_FIXTURE is not set");
        return;
    };
    let directory = PathBuf::from(directory);

    let cells = rows(&directory.join("cells.csv"));
    let totals: Vec<f64> = cells
        .iter()
        .map(|row| row[2].parse().expect("total_umi"))
        .collect();

    let genes = rows(&directory.join("genes.csv"));
    let mut sparse: Vec<Vec<(usize, f64)>> = vec![Vec::new(); genes.len()];
    for row in rows(&directory.join("counts.csv")) {
        let gene: usize = row[0].parse().expect("gene_row");
        let cell: usize = row[1].parse().expect("cell_col");
        sparse[gene].push((cell, row[2].parse().expect("count")));
    }
    for column in &mut sparse {
        column.sort_by_key(|&(cell, _)| cell);
    }

    let mut intercept_error: Vec<f64> = Vec::new();
    let mut regime_disagreements: Vec<String> = Vec::new();

    for (index, row) in genes.iter().enumerate() {
        let name = &row[0];
        let expected_intercept: f64 = row[4].parse().expect("mle_intercept");
        let expected_theta: f64 = if row[2].to_lowercase().starts_with("inf") {
            f64::INFINITY
        } else {
            row[2].parse().expect("theta")
        };

        let fit = sctransform_rs::fit_offset_model(&sparse[index], &totals);
        intercept_error.push((fit.intercept - expected_intercept).abs() / expected_intercept.abs());
        if fit.theta.is_infinite() != expected_theta.is_infinite() {
            regime_disagreements.push(format!(
                "{name}: theta {:.6e} vs glmGamPoi {:.6e}",
                fit.theta, expected_theta
            ));
        }
    }

    intercept_error.sort_by(f64::total_cmp);
    let n = intercept_error.len();
    eprintln!(
        "beta stage over {n} genes: intercept relative error median {:.3e}, max {:.3e}",
        intercept_error[n / 2],
        intercept_error[n - 1]
    );
    eprintln!(
        "genes whose overdispersion regime disagrees end to end: {}",
        regime_disagreements.len()
    );
    for line in regime_disagreements.iter().take(8) {
        eprintln!("  {line}");
    }

    assert!(
        intercept_error[n / 2] < 1e-10,
        "median intercept relative error {:e}",
        intercept_error[n / 2]
    );
}
