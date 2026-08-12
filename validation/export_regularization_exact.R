#!/usr/bin/env Rscript

# Full-precision export of `reg_model_pars`'s intermediates, link by link.
#
# The regularization stage now carries all of the port's remaining
# disagreement, and it is a chain: dispersion_par, an outlier rule scored
# against two offset bin grids, three separate Poisson-exclusion criteria, a
# Sheather-Jones bandwidth, and finally `ksmooth`. Comparing only the smoothed
# output says which chain is wrong and nothing about where, and two earlier
# stages of this port were fixed only after the first differing link was found
# rather than the loudest one.
#
# So every link is exported separately, at %.17g: the gene set entering the
# stage, the dispersion_par computed from it, the outlier mask, the set left
# after both filters, the bandwidth derived from *that* set, the clamped
# evaluation points, and the smoothed result.
#
# GPL-3 validation probe. Not part of the Rust runtime.

validation_library <- Sys.getenv("BIOLANG_VALIDATION_R_LIB", unset = "")
if (nzchar(validation_library)) {
  .libPaths(c(normalizePath(validation_library, mustWork = TRUE), .libPaths()))
}
suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: Rscript export_regularization_exact.R INPUT_MEX OUTPUT_DIR")
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

set.seed(1448145L)
clip <- sqrt(ncol(umi) / 30)
oracle <- sctransform::vst(
  umi, vst.flavor = "v2",
  n_cells = min(5000L, ncol(umi)), n_genes = min(2000L, nrow(umi)),
  min_cells = 5L, res_clip_range = c(-clip, clip),
  return_cell_attr = TRUE, return_gene_attr = TRUE,
  return_corrected_umi = FALSE, verbosity = 0
)

g <- function(x) sprintf("%.17g", x)

# vst drops genes below min_cells before modelling; work on the same axis.
kept <- rownames(oracle$model_pars_fit)
umi_kept <- umi[kept, , drop = FALSE]

row_gmean <- getFromNamespace("row_gmean", "sctransform")
row_var <- getFromNamespace("row_var", "sctransform")
is_outlier <- getFromNamespace("is_outlier", "sctransform")

genes_log_gmean <- log10(row_gmean(umi_kept, eps = 1))
step1 <- rownames(oracle$model_pars)
genes_log_gmean_step1 <- oracle$genes_log_gmean_step1
stopifnot(identical(names(genes_log_gmean_step1), step1))

model_pars <- oracle$model_pars

# Link 1: the quantity actually smoothed.
dispersion_par <- log10(1 + 10^genes_log_gmean_step1 / model_pars[, "theta"])

# Link 2: the outlier rule, per column and combined, exactly as reg_model_pars
# assembles it -- theta replaced by dispersion_par before scoring.
scored <- cbind(dispersion_par, model_pars[, colnames(model_pars) != "theta", drop = FALSE])
outlier_columns <- apply(scored, 2, function(y) is_outlier(y, genes_log_gmean_step1))
outliers <- apply(outlier_columns, 1, any)
is_theta_inf <- !is.finite(model_pars[, "theta"])
outliers_v2 <- outliers | is_theta_inf

# Link 3: Poisson exclusion, all three criteria.
genes_amean <- Matrix::rowMeans(umi_kept)
genes_var <- row_var(umi_kept)
all_poisson_genes <- union(
  names(genes_amean)[(genes_var - genes_amean) <= 0],
  names(genes_amean)[genes_amean < 0.001]
)

surviving <- setdiff(step1[!outliers_v2], all_poisson_genes)

# Link 4: the bandwidth, from the surviving set rather than the full one.
bw <- bw.SJ(genes_log_gmean_step1[surviving]) * 3

# Link 5: the evaluation points.
x_points <- pmax(genes_log_gmean, min(genes_log_gmean_step1[surviving]))
x_points <- pmin(x_points, max(genes_log_gmean_step1[surviving]))

writeLines(
  c("gene,log_gmean_step1,theta,intercept,dispersion_par,outlier_any,theta_inf,poisson,survives",
    paste(step1,
          g(as.numeric(genes_log_gmean_step1)),
          g(model_pars[, "theta"]),
          g(model_pars[, "(Intercept)"]),
          g(as.numeric(dispersion_par)),
          as.integer(outliers),
          as.integer(is_theta_inf),
          as.integer(step1 %in% all_poisson_genes),
          as.integer(step1 %in% surviving),
          sep = ",")),
  file.path(output, "step1-genes.csv")
)

writeLines(
  c("gene,log_gmean,x_point,theta_fit,intercept_fit",
    paste(kept,
          g(as.numeric(genes_log_gmean)),
          g(as.numeric(x_points)),
          g(oracle$model_pars_fit[, "theta"]),
          g(oracle$model_pars_fit[, "(Intercept)"]),
          sep = ",")),
  file.path(output, "fitted.csv")
)

writeLines(
  c("key,value",
    paste("kept_genes", length(kept), sep = ","),
    paste("step1_genes", length(step1), sep = ","),
    paste("outliers", sum(outliers), sep = ","),
    paste("theta_inf", sum(is_theta_inf), sep = ","),
    paste("poisson_in_step1", sum(step1 %in% all_poisson_genes), sep = ","),
    paste("poisson_all", length(all_poisson_genes), sep = ","),
    paste("surviving", length(surviving), sep = ","),
    paste("bw_sj_times_adjust", g(bw), sep = ","),
    paste("bw_sj", g(bw / 3), sep = ","),
    paste("ksmooth_scale", g(0.3706506), sep = ","),
    paste("sctransform_version", as.character(packageVersion("sctransform")), sep = ",")),
  file.path(output, "manifest.csv")
)

cat(sprintf(
  "REGULARIZATION_EXACT_OK step1=%d outliers=%d poisson=%d surviving=%d bw=%s\n",
  length(step1), sum(outliers), sum(step1 %in% all_poisson_genes), length(surviving), g(bw)
))
