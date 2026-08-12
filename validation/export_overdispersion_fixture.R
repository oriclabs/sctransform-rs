#!/usr/bin/env Rscript

# Export the overdispersion estimator's inputs and outputs, per gene.
#
# The estimator is the last large disagreement with the reference, and it sits
# behind two stages that can hide it: which genes get fit, and what mu the beta
# stage produced. This script cuts it out of the pipeline. It exports the exact
# (y, mu) glmGamPoi optimised over and the theta it returned, so a port can be
# scored on the estimator alone -- if the Rust answer differs given identical
# inputs, the estimator is wrong, with nothing upstream to blame.
#
# Genes are sampled across the expression range rather than taken as a prefix,
# because theta is easy at high expression and hard in the sparse tail, and a
# prefix would report the easy half.
#
# GPL-3 validation probe. Not part of the Rust runtime.

validation_library <- Sys.getenv("BIOLANG_VALIDATION_R_LIB", unset = "")
if (nzchar(validation_library)) {
  .libPaths(c(normalizePath(validation_library, mustWork = TRUE), .libPaths()))
}
suppressPackageStartupMessages({
  library(Matrix)
  library(glmGamPoi)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop("usage: Rscript export_overdispersion_fixture.R INPUT_MEX OUTPUT_DIR [N_GENES]")
}
input <- normalizePath(args[[1L]], mustWork = TRUE)
output <- args[[2L]]
n_genes_wanted <- if (length(args) >= 3L) as.integer(args[[3L]]) else 200L
dir.create(output, recursive = TRUE, showWarnings = FALSE)

read_labels <- function(path) read.delim(path, header = FALSE, stringsAsFactors = FALSE)
umi <- Matrix::readMM(file.path(input, "matrix.mtx"))
features <- read_labels(file.path(input, "features.tsv"))
barcodes <- read_labels(file.path(input, "barcodes.tsv"))
rownames(umi) <- make.unique(as.character(features[[min(2L, ncol(features))]]))
colnames(umi) <- as.character(barcodes[[1L]])
umi <- as(umi, "dgCMatrix")

# The offset must be computed before the min_cells gene filter, not after.
# `vst` builds `cell_attr` at line 62 of its body and only drops genes at line
# 70, so its `log_umi` is a column sum over every input gene. Computing it on
# the filtered matrix instead gives a fixture that is self-consistent -- a port
# fed the same wrong totals will agree with it perfectly -- while not being the
# pipeline's inputs, which is the failure mode a fixture exists to prevent.
min_cells <- 5L
cell_totals <- Matrix::colSums(umi)
log_umi <- log(cell_totals)
umi <- umi[Matrix::rowSums(umi > 0) >= min_cells, , drop = FALSE]

row_gmean <- getFromNamespace("row_gmean", "sctransform")
log_gmean <- log10(row_gmean(umi, eps = 1))

# Even coverage of the expression range: rank by geometric mean and take
# equally spaced ranks.
ranked <- order(log_gmean)
picks <- unique(round(seq(1, length(ranked), length.out = n_genes_wanted)))
selected <- rownames(umi)[ranked[picks]]
subset <- umi[selected, , drop = FALSE]

# glm_gp's own Beta and Mu cannot be used to reconstruct what the estimator
# saw. It fits beta, estimates the overdispersion from *that* Mu, shrinks the
# dispersions, and then fits beta a second time -- so the Beta and Mu it
# returns are from the second fit, and reconstructing mu from them would score
# the port against inputs the reference never optimised over. Walking the
# stages explicitly exports the first-stage beta, which is the one that
# defines the estimator's mu.
model_matrix <- matrix(1, nrow = ncol(subset), ncol = 1L,
                       dimnames = list(NULL, "Intercept"))
offset_matrix <- matrix(log_umi, nrow = nrow(subset), ncol = ncol(subset), byrow = TRUE)

disp_init <- glmGamPoi:::estimate_dispersions_roughly(subset, model_matrix, offset_matrix)
groups <- glmGamPoi:::get_groups_for_model_matrix(model_matrix)
stopifnot(!is.null(groups))
beta_group_init <- glmGamPoi:::estimate_betas_roughly_group_wise(subset, offset_matrix, groups)
beta_res <- glmGamPoi:::estimate_betas_group_wise(
  subset, offset_matrix = offset_matrix, dispersions = disp_init,
  beta_group_init = beta_group_init, groups = groups, model_matrix = model_matrix
)
Mu <- glmGamPoi:::calculate_mu(beta_res$Beta, model_matrix, offset_matrix)
disp_est <- glmGamPoi::overdispersion_mle(
  subset, Mu, model_matrix = model_matrix,
  do_cox_reid_adjustment = TRUE, subsample = FALSE
)$estimate

fit <- glmGamPoi::glm_gp(
  data = subset,
  design = ~1,
  offset = log_umi,
  size_factors = FALSE
)
# The staged walk must reproduce what glm_gp reports, or the stages above are
# not the ones glm_gp ran.
stopifnot(isTRUE(all.equal(disp_est, unname(fit$overdispersions), tolerance = 1e-12)))

g <- function(x) sprintf("%.17g", x)

# Per gene: the estimate, and the first-stage intercept that defines the mu the
# optimiser actually saw.
writeLines(
  c("gene,gene_index,theta,overdispersion,mle_intercept,disp_init,final_intercept,n_cells",
    paste(selected,
          match(selected, rownames(umi)) - 1L,
          g(1 / disp_est),
          g(disp_est),
          g(beta_res$Beta[, 1]),
          g(disp_init),
          g(fit$Beta[, 1]),
          ncol(subset),
          sep = ",")),
  file.path(output, "genes.csv")
)

# The full (y, mu) the estimator saw, long form. Only non-zero counts are
# listed; mu is dense but fully determined by intercept + offset, so it is
# reconstructible from cells.csv rather than stored per gene.
nz <- Matrix::summary(as(subset, "dgCMatrix"))
writeLines(
  c("gene_row,cell_col,count",
    paste(nz$i - 1L, nz$j - 1L, g(nz$x), sep = ",")),
  file.path(output, "counts.csv")
)
writeLines(
  c("cell,cell_index,total_umi,log_umi",
    paste(colnames(subset), seq_len(ncol(subset)) - 1L,
          g(cell_totals), g(log_umi), sep = ",")),
  file.path(output, "cells.csv")
)
writeLines(
  c("key,value",
    paste("genes", nrow(subset), sep = ","),
    paste("cells", ncol(subset), sep = ","),
    paste("glmGamPoi_version", as.character(packageVersion("glmGamPoi")), sep = ","),
    paste("sctransform_version", as.character(packageVersion("sctransform")), sep = ","),
    paste("cr_correction_factor", g(0.99), sep = ",")),
  file.path(output, "manifest.csv")
)

cat(sprintf("OVERDISPERSION_FIXTURE_OK genes=%d cells=%d\n", nrow(subset), ncol(subset)))
