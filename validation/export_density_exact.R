#!/usr/bin/env Rscript

# Full-precision export of the upstream step-one density intermediates.
#
# `export_upstream_sampling.R` writes its numbers with write.csv, which stops
# at 15 significant digits. That is not enough to answer the question this
# stage actually poses: whether a difference is a structural mistake or the
# last bit of a double. Every number here is printed with %.17g, so a
# comparison against the Rust engine measures the real disagreement rather
# than the printer's.
#
# The chain is exported at each link -- input coordinate, grid, interpolated
# density, sampling weight -- because a divergence at the end says nothing
# about where it began. Finding that the geometric mean, not the density
# convolution, was the first link to differ is what this script is for.
#
# GPL-3 validation probe. Not part of the Rust runtime.

validation_library <- Sys.getenv("BIOLANG_VALIDATION_R_LIB", unset = "")
if (nzchar(validation_library)) {
  .libPaths(c(normalizePath(validation_library, mustWork = TRUE), .libPaths()))
}
suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: Rscript export_density_exact.R INPUT_MEX OUTPUT_DIR")
}
input <- normalizePath(args[[1L]], mustWork = TRUE)
output <- args[[2L]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)

read_labels <- function(path) read.delim(path, header = FALSE, stringsAsFactors = FALSE)
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

fit <- density(log_step1, bw = "nrd", adjust = 1)
interpolated <- approx(fit$x, fit$y, xout = log_step1)$y
probability <- 1 / (interpolated + .Machine$double.eps)

g <- function(x) sprintf("%.17g", x)

writeLines(
  c("gene,gene_index,log_geometric_mean,density,sampling_weight",
    paste(genes_step1,
          match(genes_step1, rownames(umi)) - 1L,
          g(as.numeric(log_step1)),
          g(interpolated),
          g(probability),
          sep = ",")),
  file.path(output, "density-exact.csv")
)

# The 512-point grid before interpolation, so a disagreement can be attributed
# to the convolution rather than to approx().
writeLines(
  c("index,x,y", paste(seq_along(fit$x) - 1L, g(fit$x), g(fit$y), sep = ",")),
  file.path(output, "grid-exact.csv")
)

writeLines(
  c("key,value",
    paste("bandwidth", g(fit$bw), sep = ","),
    paste("n_grid", length(fit$x), sep = ","),
    paste("candidates", length(genes_step1), sep = ","),
    paste("data_min", g(min(log_step1)), sep = ","),
    paste("data_max", g(max(log_step1)), sep = ","),
    paste("from", g(min(fit$x)), sep = ","),
    paste("to", g(max(fit$x)), sep = ","),
    paste("sctransform_version", as.character(packageVersion("sctransform")), sep = ",")),
  file.path(output, "density-manifest.csv")
)

cat(sprintf("DENSITY_EXACT_OK candidates=%d bw=%s\n", length(genes_step1), g(fit$bw)))
