#!/usr/bin/env Rscript

# GPL-3 validation probe for the upstream sctransform step-one population and
# weighted gene sample. This script is not part of the Rust runtime.

validation_library <- Sys.getenv("BIOLANG_VALIDATION_R_LIB", unset = "")
if (nzchar(validation_library)) {
  .libPaths(c(normalizePath(validation_library, mustWork = TRUE), .libPaths()))
}
suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: Rscript export_upstream_sampling.R INPUT_MEX OUTPUT_DIR")
}
input <- normalizePath(args[[1L]], mustWork = TRUE)
output <- args[[2L]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)

read_labels <- function(path) {
  read.delim(path, header = FALSE, stringsAsFactors = FALSE)
}
umi <- Matrix::readMM(file.path(input, "matrix.mtx"))
features <- read_labels(file.path(input, "features.tsv"))
barcodes <- read_labels(file.path(input, "barcodes.tsv"))
rownames(umi) <- make.unique(as.character(features[[min(2L, ncol(features))]]))
colnames(umi) <- as.character(barcodes[[1L]])
umi <- as(umi, "dgCMatrix")

min_cells <- 5L
genes <- rownames(umi)[Matrix::rowSums(umi > 0) >= min_cells]
umi <- umi[genes, , drop = FALSE]
row_gmean <- getFromNamespace("row_gmean", "sctransform")
row_var <- getFromNamespace("row_var", "sctransform")
genes_log_gmean <- log10(row_gmean(umi, eps = 1))

set.seed(1448145L)
cells_step1 <- sample(colnames(umi), size = min(5000L, ncol(umi)))
detected_step1 <- Matrix::rowSums(umi[, cells_step1, drop = FALSE] > 0)
genes_step1 <- genes[detected_step1 >= min_cells]
genes_amean <- Matrix::rowMeans(umi)
genes_var <- row_var(umi)
genes_step1 <- genes_step1[(genes_var - genes_amean)[genes_step1] > 0]
log_step1 <- genes_log_gmean[genes_step1]

density_fit <- density(log_step1, bw = "nrd", adjust = 1)
sampling_probability <- 1 / (
  approx(density_fit$x, density_fit$y, xout = log_step1)$y + .Machine$double.eps
)
selected <- sample(genes_step1, size = min(2000L, length(genes_step1)),
                   prob = sampling_probability)

write.csv(
  data.frame(
    gene = genes_step1,
    gene_index = match(genes_step1, rownames(umi)) - 1L,
    log_geometric_mean = as.numeric(log_step1),
    sampling_weight = sampling_probability,
    stringsAsFactors = FALSE
  ),
  file.path(output, "candidates.csv"), row.names = FALSE, quote = TRUE
)
write.csv(
  data.frame(
    rank = seq_along(selected),
    gene = selected,
    gene_index = match(selected, rownames(umi)) - 1L,
    stringsAsFactors = FALSE
  ),
  file.path(output, "selected.csv"), row.names = FALSE, quote = TRUE
)
write.csv(
  data.frame(
    cells = ncol(umi),
    genes = nrow(umi),
    candidate_genes = length(genes_step1),
    selected_genes = length(selected),
    seed = 1448145L,
    density_bandwidth = density_fit$bw,
    sctransform_version = as.character(packageVersion("sctransform")),
    stringsAsFactors = FALSE
  ),
  file.path(output, "manifest.csv"), row.names = FALSE, quote = TRUE
)

