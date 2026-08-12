// Standalone GPL-3 SCTransform process provider.
// Copyright (C) 2026 ORIC Labs.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, version 3 only.

use sctransform_rs::{sctransform, GeneColumns, SctOptions};
use std::collections::HashMap;
use std::env;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

const VERSION: &str = env!("CARGO_PKG_VERSION");
const ENGINE: &str = "sctransform-rs-gpl-port";
const PROTOCOL: &str = "1";

#[derive(Debug)]
struct RunArguments {
    input: PathBuf,
    output: PathBuf,
    n_features: Option<usize>,
    threads: usize,
    probe_genes: usize,
    probe_cells: usize,
    write_matrix: bool,
    write_binary_matrix: bool,
    regress_mito: bool,
}

struct MexInput {
    matrix: GeneColumns,
    genes: Vec<String>,
    cells: Vec<String>,
}

fn main() -> ExitCode {
    match real_main() {
        Ok(()) => ExitCode::SUCCESS,
        Err((code, message)) => {
            eprintln!("error: {message}");
            ExitCode::from(code)
        }
    }
}

fn real_main() -> Result<(), (u8, String)> {
    let arguments: Vec<String> = env::args().skip(1).collect();
    match arguments.first().map(String::as_str) {
        Some("--version" | "version") => {
            println!("bl-sctransform-gpl {VERSION} ({ENGINE}, GPL-3.0-only)");
            Ok(())
        }
        Some("license") => {
            println!(
                "bl-sctransform-gpl {VERSION}\n\
                 Copyright (C) 2026 ORIC Labs\n\
                 License: GNU GPL version 3 only\n\
                 This program comes with ABSOLUTELY NO WARRANTY.\n\
                 Source: https://github.com/oriclabs/sctransform-rs"
            );
            Ok(())
        }
        Some("run") => run(parse_run_arguments(&arguments[1..])?),
        Some("--help" | "-h" | "help") | None => {
            print_help();
            Ok(())
        }
        Some(command) => Err((2, format!("unknown command '{command}'; use --help"))),
    }
}

fn print_help() {
    println!(
        "bl-sctransform-gpl {VERSION} - GPL-3 external SCTransform provider\n\n\
         USAGE:\n  \
         bl-sctransform-gpl run --input DIR --output DIR [OPTIONS]\n\n\
         OPTIONS:\n  \
         --n-features N    Keep the top N residual-variance features\n  \
         --threads N       Worker threads; 0 chooses automatically [default: 0]\n  \
         --probe-genes N   Genes written to residuals.csv [default: 3000]\n  \
         --probe-cells N   Cells written to residuals.csv [default: 64]\n  \
         --write-matrix    Also write the complete kept residual matrix\n\n\
         --write-binary-matrix  Write a compact row-major f64 residual matrix\n  \
         --regress-mito    Regress the per-cell fraction from MT- genes\n\n\
         COMMANDS:\n  \
         license           Show copyright, license, warranty, and source notice\n  \
         version           Show executable and engine versions"
    );
}

fn parse_run_arguments(values: &[String]) -> Result<RunArguments, (u8, String)> {
    let mut input = None;
    let mut output = None;
    let mut n_features = None;
    let mut threads = 0;
    let mut probe_genes = 3000;
    let mut probe_cells = 64;
    let mut write_matrix = false;
    let mut write_binary_matrix = false;
    let mut regress_mito = false;
    let mut index = 0;
    while index < values.len() {
        let option = values[index].as_str();
        if option == "--write-matrix"
            || option == "--write-binary-matrix"
            || option == "--regress-mito"
        {
            if option == "--write-matrix" {
                write_matrix = true;
            } else if option == "--write-binary-matrix" {
                write_binary_matrix = true;
            } else {
                regress_mito = true;
            }
            index += 1;
            continue;
        }
        let value = values
            .get(index + 1)
            .ok_or_else(|| (2, format!("{option} requires a value")))?;
        match option {
            "--input" => input = Some(PathBuf::from(value)),
            "--output" => output = Some(PathBuf::from(value)),
            "--n-features" => n_features = Some(parse_usize(option, value)?),
            "--threads" => threads = parse_usize(option, value)?,
            "--probe-genes" => probe_genes = parse_usize(option, value)?,
            "--probe-cells" => probe_cells = parse_usize(option, value)?,
            _ => return Err((2, format!("unknown option '{option}'"))),
        }
        index += 2;
    }
    Ok(RunArguments {
        input: input.ok_or_else(|| (2, "--input is required".to_string()))?,
        output: output.ok_or_else(|| (2, "--output is required".to_string()))?,
        n_features,
        threads,
        probe_genes,
        probe_cells,
        write_matrix,
        write_binary_matrix,
        regress_mito,
    })
}

fn parse_usize(option: &str, value: &str) -> Result<usize, (u8, String)> {
    value.parse::<usize>().map_err(|_| {
        (
            2,
            format!("{option} expects a non-negative integer, got '{value}'"),
        )
    })
}

fn run(arguments: RunArguments) -> Result<(), (u8, String)> {
    eprintln!("engine={ENGINE} version={VERSION} license=GPL-3.0-only protocol={PROTOCOL}");
    prepare_output_directory(&arguments.output)?;
    let input = read_mex(&arguments.input).map_err(|message| (3, message))?;
    let latent_covariates = if arguments.regress_mito {
        vec![mitochondrial_fraction(&input)]
    } else {
        Vec::new()
    };
    let options = SctOptions {
        n_variable_features: arguments.n_features,
        threads: arguments.threads,
        latent_covariates,
        ..SctOptions::default()
    };
    let started = Instant::now();
    let result = sctransform(&input.matrix, &options);
    let elapsed = started.elapsed().as_secs_f64();

    write_outputs(&arguments, &input, &options, &result, elapsed)
        .map_err(|message| (5, message))?;
    println!(
        "SCTRANSFORM_GPL_OK engine={ENGINE} cells={} genes={} modelled={} elapsed={elapsed:.6}s",
        input.matrix.n_cells,
        input.matrix.n_genes(),
        result.kept_genes.len()
    );
    Ok(())
}

fn mitochondrial_fraction(input: &MexInput) -> Vec<f64> {
    let mut totals = vec![0.0; input.matrix.n_cells];
    let mut mitochondrial = vec![0.0; input.matrix.n_cells];
    for gene in 0..input.matrix.n_genes() {
        let start = input.matrix.starts[gene];
        let end = input.matrix.starts[gene + 1];
        let is_mitochondrial = input.genes[gene].starts_with("MT-");
        for index in start..end {
            let cell = input.matrix.cells[index] as usize;
            let count = input.matrix.counts[index];
            totals[cell] += count;
            if is_mitochondrial {
                mitochondrial[cell] += count;
            }
        }
    }
    mitochondrial
        .into_iter()
        .zip(totals)
        .map(|(mito, total)| if total > 0.0 { mito / total } else { 0.0 })
        .collect()
}

fn prepare_output_directory(path: &Path) -> Result<(), (u8, String)> {
    if path.exists() {
        let mut entries = fs::read_dir(path)
            .map_err(|error| (5, format!("cannot inspect {}: {error}", path.display())))?;
        if entries.next().is_some() {
            return Err((
                5,
                format!(
                    "output directory must be absent or empty: {}",
                    path.display()
                ),
            ));
        }
    } else {
        fs::create_dir_all(path)
            .map_err(|error| (5, format!("cannot create {}: {error}", path.display())))?;
    }
    Ok(())
}

fn read_mex(path: &Path) -> Result<MexInput, String> {
    let matrix_path = path.join("matrix.mtx");
    let features_path = path.join("features.tsv");
    let barcodes_path = path.join("barcodes.tsv");
    for required in [&matrix_path, &features_path, &barcodes_path] {
        if !required.is_file() {
            return Err(format!(
                "required input file is missing: {}",
                required.display()
            ));
        }
    }

    let genes = read_labels(&features_path, true)?;
    let cells = read_labels(&barcodes_path, false)?;
    let reader = BufReader::new(
        File::open(&matrix_path)
            .map_err(|error| format!("cannot open {}: {error}", matrix_path.display()))?,
    );
    let mut dimensions = None;
    let mut entries: Vec<(usize, usize, f64)> = Vec::new();
    for line in reader.lines() {
        let line =
            line.map_err(|error| format!("cannot read {}: {error}", matrix_path.display()))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('%') {
            continue;
        }
        if dimensions.is_none() {
            let values: Vec<&str> = trimmed.split_whitespace().collect();
            if values.len() != 3 {
                return Err("Matrix Market dimensions must contain rows, columns, and nnz".into());
            }
            let rows = parse_dimension(values[0], "row count")?;
            let columns = parse_dimension(values[1], "column count")?;
            let nnz = parse_dimension(values[2], "non-zero count")?;
            dimensions = Some((rows, columns, nnz));
            entries.reserve(nnz);
            continue;
        }
        let values: Vec<&str> = trimmed.split_whitespace().collect();
        if values.len() != 3 {
            return Err(format!("invalid Matrix Market entry: {trimmed}"));
        }
        let gene = parse_dimension(values[0], "gene index")?;
        let cell = parse_dimension(values[1], "cell index")?;
        if gene == 0 || cell == 0 {
            return Err("Matrix Market indices must be one-based".into());
        }
        let count = values[2]
            .parse::<f64>()
            .map_err(|_| format!("invalid count in Matrix Market entry: {trimmed}"))?;
        if !count.is_finite() || count < 0.0 {
            return Err(format!("counts must be finite and non-negative: {trimmed}"));
        }
        if count != 0.0 {
            entries.push((cell - 1, gene - 1, count));
        }
    }
    let (n_genes, n_cells, declared_nnz) =
        dimensions.ok_or_else(|| "Matrix Market file has no dimensions".to_string())?;
    if genes.len() != n_genes {
        return Err(format!(
            "features.tsv has {} rows but matrix has {n_genes} genes",
            genes.len()
        ));
    }
    if cells.len() != n_cells {
        return Err(format!(
            "barcodes.tsv has {} rows but matrix has {n_cells} cells",
            cells.len()
        ));
    }
    if entries.len() > declared_nnz {
        return Err("Matrix Market contains more entries than declared".into());
    }
    if entries
        .iter()
        .any(|&(cell, gene, _)| cell >= n_cells || gene >= n_genes)
    {
        return Err("Matrix Market entry exceeds declared dimensions".into());
    }
    entries.sort_unstable_by(|left, right| (left.0, left.1).cmp(&(right.0, right.1)));
    let matrix = GeneColumns::from_cell_major(n_cells, n_genes, |emit| {
        for &(cell, gene, count) in &entries {
            emit(cell, gene, count);
        }
    });
    Ok(MexInput {
        matrix,
        genes,
        cells,
    })
}

fn parse_dimension(value: &str, label: &str) -> Result<usize, String> {
    value
        .parse::<usize>()
        .map_err(|_| format!("invalid {label}: {value}"))
}

fn read_labels(path: &Path, prefer_second: bool) -> Result<Vec<String>, String> {
    let reader = BufReader::new(
        File::open(path).map_err(|error| format!("cannot open {}: {error}", path.display()))?,
    );
    reader
        .lines()
        .enumerate()
        .map(|(index, line)| {
            let line = line.map_err(|error| format!("cannot read {}: {error}", path.display()))?;
            let columns: Vec<&str> = line.split('\t').collect();
            let label = if prefer_second && columns.len() >= 2 {
                columns[1]
            } else {
                columns.first().copied().unwrap_or("")
            };
            if label.is_empty() {
                Err(format!("empty label at {}:{}", path.display(), index + 1))
            } else {
                Ok(label.to_string())
            }
        })
        .collect()
}

fn write_outputs(
    arguments: &RunArguments,
    input: &MexInput,
    options: &SctOptions,
    result: &sctransform_rs::SctResult,
    elapsed: f64,
) -> Result<(), String> {
    let output = &arguments.output;
    let mut slot_by_gene = HashMap::with_capacity(result.kept_genes.len());
    for (slot, gene) in result.kept_genes.iter().copied().enumerate() {
        slot_by_gene.insert(gene, slot);
    }

    let mut genes = writer(output.join("genes.csv"))?;
    writeln!(
        genes,
        "gene,gene_index,theta,intercept,residual_variance,log_geometric_mean"
    )
    .map_err(write_error)?;
    for (slot, gene) in result.kept_genes.iter().copied().enumerate() {
        writeln!(
            genes,
            "{},{gene},{},{},{},{}",
            csv(&input.genes[gene]),
            result.theta[slot],
            result.intercept[slot],
            result.residual_variance[slot],
            result.log_geometric_mean[slot]
        )
        .map_err(write_error)?;
    }

    let mut ranking = writer(output.join("ranking.csv"))?;
    writeln!(ranking, "rank,gene,gene_index").map_err(write_error)?;
    for (rank, gene) in result.ranked_genes.iter().copied().enumerate() {
        writeln!(ranking, "{},{},{gene}", rank + 1, csv(&input.genes[gene]))
            .map_err(write_error)?;
    }

    let mut fit_genes = writer(output.join("fit-genes.csv"))?;
    writeln!(fit_genes, "gene,gene_index,raw_theta,raw_intercept").map_err(write_error)?;
    for (slot, gene) in result.fit_genes.iter().copied().enumerate() {
        writeln!(
            fit_genes,
            "{},{gene},{},{}",
            csv(&input.genes[gene]),
            result.raw_theta[slot],
            result.raw_intercept[slot]
        )
        .map_err(write_error)?;
    }

    let mut sampling = writer(output.join("sampling.csv"))?;
    writeln!(sampling, "gene,gene_index,sampling_weight").map_err(write_error)?;
    for (gene, weight) in result
        .fit_candidates
        .iter()
        .copied()
        .zip(result.fit_candidate_weights.iter().copied())
    {
        writeln!(sampling, "{},{gene},{weight}", csv(&input.genes[gene])).map_err(write_error)?;
    }

    let mut fit_cells = writer(output.join("fit-cells.csv"))?;
    writeln!(fit_cells, "cell,cell_index").map_err(write_error)?;
    for cell in result.fit_cells.iter().copied() {
        writeln!(fit_cells, "{},{cell}", csv(&input.cells[cell])).map_err(write_error)?;
    }

    let probe_genes = result
        .ranked_genes
        .iter()
        .filter_map(|gene| slot_by_gene.get(gene).map(|slot| (*gene, *slot)))
        .take(arguments.probe_genes)
        .collect::<Vec<_>>();
    let probe_cells = arguments.probe_cells.min(input.matrix.n_cells);
    let width = result.kept_genes.len();
    let mut residuals = writer(output.join("residuals.csv"))?;
    writeln!(residuals, "gene,cell,residual").map_err(write_error)?;
    for (gene, slot) in &probe_genes {
        for cell in 0..probe_cells {
            writeln!(
                residuals,
                "{},{},{}",
                csv(&input.genes[*gene]),
                csv(&input.cells[cell]),
                result.residuals[cell * width + *slot]
            )
            .map_err(write_error)?;
        }
    }

    if arguments.write_matrix {
        write_residual_matrix(output.join("matrix.mtx"), input, result)?;
    }
    if arguments.write_binary_matrix {
        write_binary_residual_matrix(output.join("matrix.f64"), input, result)?;
    }

    let mut manifest = writer(output.join("manifest.csv"))?;
    writeln!(manifest, "implementation,engine,version,license,protocol_version,cells,genes,modelled_genes,fit_candidate_genes,sampling_bandwidth,clip,seed,cells_for_fit,genes_for_fit,min_cells,residual_probe_strategy,residual_probe_genes,residual_probe_cells,elapsed_seconds,matrix_written,binary_matrix_written,mitochondrial_fraction_regressed")
        .map_err(write_error)?;
    let clip = options
        .clip
        .unwrap_or_else(|| (input.matrix.n_cells as f64 / 30.0).sqrt());
    writeln!(
        manifest,
        "sctransform-rs,{ENGINE},{VERSION},GPL-3.0-only,{PROTOCOL},{},{},{},{},{},{clip},1448145,{},{},{},top_residual_variance,{},{},{elapsed},{},{},{}",
        input.matrix.n_cells,
        input.matrix.n_genes(),
        result.kept_genes.len(),
        result.fit_candidates.len(),
        result.sampling_bandwidth,
        options.cells_for_fit.min(input.matrix.n_cells),
        options.genes_for_fit.min(input.matrix.n_genes()),
        options.min_cells,
        probe_genes.len(),
        probe_cells,
        arguments.write_matrix,
        arguments.write_binary_matrix,
        arguments.regress_mito
    )
    .map_err(write_error)?;
    Ok(())
}

/// Neutral interchange format used at the process boundary:
/// magic `BLMATF64`, little-endian u64 rows, little-endian u64 columns, then
/// row-major IEEE-754 f64 values. BioLang's MIT reader is generic and neither
/// links to nor depends on this executable.
fn write_binary_residual_matrix(
    path: PathBuf,
    input: &MexInput,
    result: &sctransform_rs::SctResult,
) -> Result<(), String> {
    let mut output = writer(path)?;
    output.write_all(b"BLMATF64").map_err(write_error)?;
    output
        .write_all(&(input.matrix.n_cells as u64).to_le_bytes())
        .map_err(write_error)?;
    output
        .write_all(&(result.kept_genes.len() as u64).to_le_bytes())
        .map_err(write_error)?;
    for value in &result.residuals {
        output
            .write_all(&value.to_le_bytes())
            .map_err(write_error)?;
    }
    Ok(())
}

fn write_residual_matrix(
    path: PathBuf,
    input: &MexInput,
    result: &sctransform_rs::SctResult,
) -> Result<(), String> {
    let mut output = writer(path)?;
    let genes = result.kept_genes.len();
    let cells = input.matrix.n_cells;
    writeln!(output, "%%MatrixMarket matrix coordinate real general").map_err(write_error)?;
    writeln!(output, "% genes by cells; generated by {ENGINE} {VERSION}").map_err(write_error)?;
    writeln!(output, "{genes} {cells} {}", genes.saturating_mul(cells)).map_err(write_error)?;
    for gene in 0..genes {
        for cell in 0..cells {
            writeln!(
                output,
                "{} {} {}",
                gene + 1,
                cell + 1,
                result.residuals[cell * genes + gene]
            )
            .map_err(write_error)?;
        }
    }
    Ok(())
}

fn writer(path: PathBuf) -> Result<BufWriter<File>, String> {
    File::create(&path)
        .map(BufWriter::new)
        .map_err(|error| format!("cannot create {}: {error}", path.display()))
}

fn write_error(error: std::io::Error) -> String {
    format!("cannot write output: {error}")
}

fn csv(value: &str) -> String {
    if value.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csv_quotes_only_when_needed() {
        assert_eq!(csv("GAPDH"), "GAPDH");
        assert_eq!(csv("gene,one"), "\"gene,one\"");
        assert_eq!(csv("gene\"one"), "\"gene\"\"one\"");
    }

    #[test]
    fn parse_requires_input_and_output() {
        let error = parse_run_arguments(&[]).unwrap_err();
        assert_eq!(error.0, 2);
        assert!(error.1.contains("--input"));
    }
}
