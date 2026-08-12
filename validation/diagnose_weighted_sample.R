#!/usr/bin/env Rscript

# Copyright (C) 2026 ORIC Labs. GPL-3.0-only.
#
# Isolate R's seeded weighted-sampling stage from density estimation. This
# accepts provider-exported weights, advances the RNG through the same cell
# sample, and records the gene draw order for a direct sampler conformance test.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: Rscript diagnose_weighted_sample.R SAMPLING_CSV N_CELLS N_FIT OUTPUT_CSV")
}
sampling <- read.csv(args[[1L]], stringsAsFactors = FALSE, check.names = FALSE)
n_cells <- as.integer(args[[2L]])
n_fit <- as.integer(args[[3L]])
set.seed(1448145L)
invisible(sample.int(n_cells, min(5000L, n_cells), replace = FALSE))
chosen <- sample.int(
  nrow(sampling), n_fit, replace = FALSE,
  prob = as.numeric(sampling$sampling_weight)
)
write.csv(
  data.frame(gene = sampling$gene[chosen], gene_index = sampling$gene_index[chosen]),
  args[[4L]], row.names = FALSE, quote = TRUE
)
