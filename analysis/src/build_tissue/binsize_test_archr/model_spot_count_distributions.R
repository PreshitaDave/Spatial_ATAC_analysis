#!/usr/bin/env Rscript
# ============================================================================
# model_spot_count_distributions.R
#
# Statistical characterization of per-tile count distributions in spatial ATAC-seq data.
# Analogous to RNA-seq's negative-binomial framing for per-gene counts, but applied to
# spatial ATAC tile-by-spot counts.
#
# For each genomic tile, we examine its count vector across all spots and ask:
# - Is it Poisson (baseline: no overdispersion)?
# - Negative Binomial (overdispersion from biological/technical variance)?
# - Zero-Inflated (excess zeros from dropout)?
#
# CITATIONS (embedded rationale):
# Robinson & Smyth (2007, Bioinformatics) — foundational NB for sequencing counts (edgeR)
# Anders & Huber (2010, Genome Biol) — DESeq NB framework, mean-variance modeling
# Love et al. (2014, Genome Biol) — DESeq2, plotDispEsts diagnostic
# Lun & Smyth (2016, NAR) — csaw, applying NB window-binning to ChIP/ATAC
# Buenrostro et al. (2015, Nature) — original scATAC; diploid sparsity basis
# Risso et al. (2018, Nat Commun) — ZINB-WaVE, zero-inflation in sparse single-cell data
# Svensson (2020, Nat Biotech) — important caveat: NB alone may explain zeros
# Granja et al. (2021, Nat Genet) — ArchR context: TF-IDF/LSI, not explicit count model
#
# Usage:
#   Rscript model_spot_count_distributions.R <tissue> <tilesize>
#   Example: Rscript model_spot_count_distributions.R lowseq_489 500
#
# ============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(Matrix)
  library(fitdistrplus)
  library(MASS)
  library(pscl)
  library(ggplot2)
})

log_msg <- function(tag, msg) {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), tag, msg))
}

# ============================================================================
# Parse arguments and setup paths
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  log_msg("error", "Usage: Rscript model_spot_count_distributions.R <tissue> <tilesize>")
  stop("Missing arguments")
}

tissue <- args[1]
tilesize <- as.integer(args[2])

log_msg("start", sprintf("===== Modeling per-tile count distributions: %s, %dbp =====",
                         tissue, tilesize))

# MUST SET GENOME BEFORE USING ArchR
set.seed(42)
addArchRGenome("hg38")

# Setup paths
project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
arrow_file <- file.path(project_root, "Data/01_inputs/arrow/arrow_not_binarize",
                        sprintf("%s_%dbp.arrow", tissue, tilesize))

output_dir <- file.path(project_root, "analysis/binsize_comparison/distribution_modeling")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 1. Load raw (non-binarized) tile counts from arrow
# ============================================================================

log_msg("step", "Loading non-binarized TileMatrix from arrow file...")

# WHY binarize=FALSE:
# - Binarized matrices (TRUE/FALSE per tile) discard count information
# - All our distributional models (Poisson, NB, ZIP, etc.) require actual counts
# - This is the foundational requirement: we need raw integer counts to fit any model
if (!file.exists(arrow_file)) {
  log_msg("error", sprintf("Arrow file not found: %s", arrow_file))
  stop("Arrow file missing")
}

tryCatch({
  tile_matrix <- ArchR::getMatrixFromArrow(
    ArrowFile = arrow_file,
    useMatrix = "TileMatrix",
    binarize = FALSE  # CRITICAL: must be FALSE to preserve counts
  )
  tile_mat <- assay(tile_matrix)  # Extract as dgCMatrix (sparse)
  log_msg("step", sprintf("Loaded TileMatrix: %d tiles × %d spots",
                          nrow(tile_mat), ncol(tile_mat)))
}, error = function(e) {
  log_msg("error", sprintf("Failed to load TileMatrix: %s", e$message))
  stop(e)
})

# ============================================================================
# 2. Filter to informative tiles (remove all-zero or ultra-sparse tiles)
# ============================================================================

log_msg("step", "Filtering to informative tiles...")

# WHY filtering:
# - All-zero tiles have zero variance, causing degenerate likelihood surfaces
# - Even low-count tiles (sum < 10 across all spots) have unreliable parameter estimates
# - Analogous to edgeR::filterByExpr or DESeq2's automatic low-count filtering
# - This prevents fitting models to essentially noise
row_sums <- Matrix::rowSums(tile_mat)
min_count_threshold <- 10
informative_tiles <- row_sums >= min_count_threshold

n_before <- nrow(tile_mat)
n_after <- sum(informative_tiles)
log_msg("step", sprintf("Filtered from %d tiles to %d tiles (removed %d with sum < %d)",
                        n_before, n_after, n_before - n_after, min_count_threshold))

tile_mat_filt <- tile_mat[informative_tiles, ]

# ============================================================================
# 3. Fast, vectorized characterization across all informative tiles
# ============================================================================

log_msg("step", "Computing per-tile mean, variance, and dispersion...")

# WHY vectorized (no per-tile loop):
# - Matrix operations are orders of magnitude faster than R loops
# - We're computing summary statistics (mean, variance, dispersion) for thousands of tiles
# - This is the "fast path" for getting a quick overview of the whole distribution

# Mean per tile (row)
tile_means <- Matrix::rowMeans(tile_mat_filt)

# Variance per tile using: Var = E[X^2] - E[X]^2
# This avoids densifying the sparse matrix, which would kill memory
tile_mat_sq <- tile_mat_filt^2
tile_means_sq <- Matrix::rowMeans(tile_mat_sq)
tile_vars <- tile_means_sq - tile_means^2

# Index of dispersion: D = Var / Mean
# Poisson predicts D = 1. Values >> 1 suggest overdispersion (NB-like)
dispersion_index <- tile_vars / tile_means

# Chi-squared dispersion test: asymptotically (n-1) * Var / Mean ~ chi-sq(n-1)
# Standard test in edgeR/DESeq2 literature for overdispersion vs Poisson null
n_spots <- ncol(tile_mat_filt)
test_stat <- (n_spots - 1) * dispersion_index
p_values <- pchisq(test_stat, df = n_spots - 1, lower.tail = FALSE)

# Classify tiles
poisson_consistent <- p_values > 0.05 & dispersion_index >= 0.8 & dispersion_index <= 1.2
significantly_overdispersed <- p_values < 0.05 & dispersion_index > 1.2
significantly_underdispersed <- p_values < 0.05 & dispersion_index < 0.8

pct_over <- 100 * sum(significantly_overdispersed) / n_after
pct_under <- 100 * sum(significantly_underdispersed) / n_after
pct_poisson <- 100 * sum(poisson_consistent) / n_after

log_msg("step", sprintf(
  "Dispersion: %.1f%% overdispersed, %.1f%% Poisson-consistent, %.1f%% underdispersed",
  pct_over, pct_poisson, pct_under
))

# ============================================================================
# 4. Mean-variance diagnostic plot (classic DESeq2/edgeR style)
# ============================================================================

log_msg("step", "Generating mean-variance diagnostic plot...")

# WHY this plot (seen in plotDispEsts, edgeR vignettes):
# - Visual confirmation of overdispersion across the range of tile abundances
# - Allows you to see if overdispersion is uniform (constant NB alpha) or mean-dependent
# - Poisson reference line (variance = mean) should sit below most tiles if NB is true

# Fit NB trend via method-of-moments
# For NB: E[Var] = mu + alpha * mu^2, where alpha is the dispersion parameter
# Solve for alpha using least-squares on the mean/variance trend
valid_idx <- tile_means > 0 & tile_vars > 0
mu_valid <- tile_means[valid_idx]
var_valid <- tile_vars[valid_idx]

# Fit: var ~ mu + alpha * mu^2
# Rearrange: (var - mu) ~ alpha * mu^2, so slope alpha = mean((var - mu) / mu^2)
alpha_method_of_moments <- mean((var_valid - mu_valid) / (mu_valid^2))

# Create plotting data frame
plot_df <- data.frame(
  mean = tile_means,
  variance = tile_vars,
  dispersion_class = ifelse(
    poisson_consistent, "Poisson",
    ifelse(significantly_overdispersed, "Overdispersed", "Underdispersed")
  )
)

# Generate the plot
pdf(file.path(output_dir, sprintf("%s_%dbp_meanvar_plot.pdf", tissue, tilesize)),
    width = 12, height = 10)

# Page 1: Mean-variance scatter with reference lines
p1 <- ggplot(plot_df, aes(x = mean, y = variance, color = dispersion_class)) +
  geom_point(alpha = 0.6, size = 1) +
  # Poisson reference line: variance = mean
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed",
              linewidth = 1) +
  # Fitted NB curve: variance = mean + alpha * mean^2
  geom_function(fun = function(x) x + alpha_method_of_moments * x^2,
                color = "red", linetype = "solid", linewidth = 1) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(
    values = c("Poisson" = "gray50", "Overdispersed" = "darkred", "Underdispersed" = "darkblue"),
    name = "Dispersion Class"
  ) +
  labs(
    title = sprintf("%s, %dbp: Mean-Variance Relationship", tissue, tilesize),
    x = "Mean count (log10)",
    y = "Variance (log10)",
    caption = sprintf(
      "Alpha (NB dispersion) = %.4f\nOverdispersed: %.1f%%, Poisson: %.1f%%, Underdispersed: %.1f%%",
      alpha_method_of_moments, pct_over, pct_poisson, pct_under
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.caption = element_text(hjust = 0.5, size = 9)
  )

print(p1)

# Page 2: Histogram of dispersion index
p2 <- ggplot(plot_df, aes(x = dispersion_index)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "black") +
  geom_vline(xintercept = 1, color = "blue", linetype = "dashed", linewidth = 1.5) +
  labs(
    title = sprintf("%s, %dbp: Distribution of Dispersion Index (Var/Mean)", tissue, tilesize),
    x = "Dispersion Index (D = Variance / Mean)",
    y = "Number of Tiles"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p2)

dev.off()

log_msg("step", sprintf("Saved mean-variance plot: %s_%dbp_meanvar_plot.pdf",
                        tissue, tilesize))

# ============================================================================
# 5. Representative tile sampling for detailed per-tile model fits
# ============================================================================

log_msg("step", "Sampling representative tiles for model comparison...")

# WHY sampling (not fitting all tiles):
# - Full MLE on 5 models × 5000+ tiles would be very slow
# - We want a representative sample spanning low/medium/high abundance
# - ~12-20 tiles is enough to see the pattern without excessive computation

# Select tiles by quantiles of mean abundance
n_sampled <- 16
quantile_breaks <- seq(0, 1, length.out = n_sampled + 1)
sampled_indices <- integer(0)

for (i in 1:n_sampled) {
  q_low <- quantile_breaks[i]
  q_high <- quantile_breaks[i + 1]
  mean_low <- quantile(tile_means, probs = q_low)
  mean_high <- quantile(tile_means, probs = q_high)

  # Find tiles in this quantile range
  in_range <- which(tile_means >= mean_low & tile_means <= mean_high)
  if (length(in_range) > 0) {
    # Pick one tile from this quantile (the one closest to the median of the range)
    median_mean <- (mean_low + mean_high) / 2
    closest_idx <- in_range[which.min(abs(tile_means[in_range] - median_mean))]
    sampled_indices <- c(sampled_indices, closest_idx)
  }
}

log_msg("step", sprintf("Sampled %d representative tiles (quantile-stratified)",
                        length(sampled_indices)))

# ============================================================================
# 6. Fit 5 candidate models to each sampled tile
# ============================================================================

log_msg("step", "Fitting 5 models to sampled tiles...")

# WHY these 5 models:
# 1. Poisson: baseline (null hypothesis of no overdispersion)
# 2. Negative Binomial (MASS::glm.nb or fitdistrplus): classic overdispersion model
#    (justified by RNA-seq literature: Robinson & Smyth 2007, Anders & Huber 2010)
# 3. Zero-Inflated Poisson (pscl::zeroinfl): Poisson + separate zero process
#    (captures dropout at single-spot level; Risso et al. 2018)
# 4. Zero-Inflated Negative Binomial (pscl::zeroinfl): NB + separate zero process
#    (combines overdispersion + zero-inflation; Risso et al. 2018)
# 5. Hurdle-NB (pscl::hurdle): alternative to ZINB; separate zero/count processes
#    (may fit better if zero process is qualitatively different; pscl documentation)

model_summary <- list()

for (i in seq_along(sampled_indices)) {
  idx <- sampled_indices[i]

  # Extract count vector for this tile (across all spots)
  counts <- as.numeric(tile_mat_filt[idx, ])

  # Skip if all zeros (shouldn't happen after filtering, but safety check)
  if (sum(counts) == 0) next

  # Tile name: use index position (robust, doesn't depend on rownames preservation)
  tile_name <- sprintf("tile_%d", idx)
  tile_mean <- tile_means[idx]
  tile_var <- tile_vars[idx]

  # Initialize row for this tile (all scalar values, guaranteed same length)
  fit_row <- list(
    tile_index = idx,
    tile_name = tile_name,
    n_spots = length(counts),
    mean_count = tile_mean,
    var_count = tile_var,
    n_zero = sum(counts == 0),
    pct_zero = 100 * sum(counts == 0) / length(counts)
  )

  # Fit each of 5 models and compute AIC
  aic_values <- c()
  bic_values <- c()
  model_names <- c()

  # 1. POISSON
  tryCatch({
    fit_pois <- fitdistrplus::fitdist(counts, "pois")
    aic_values <- c(aic_values, fit_pois$aic)
    bic_values <- c(bic_values, -2 * fit_pois$loglik + 1 * log(length(counts)))
    model_names <- c(model_names, "Poisson")
  }, error = function(e) {
    log_msg("warn", sprintf("Poisson fit failed for tile %s: %s", tile_name, e$message))
  })

  # 2. NEGATIVE BINOMIAL
  tryCatch({
    fit_nb <- fitdistrplus::fitdist(counts, "nbinom")
    aic_values <- c(aic_values, fit_nb$aic)
    bic_values <- c(bic_values, -2 * fit_nb$loglik + 2 * log(length(counts)))
    model_names <- c(model_names, "NegBin")
  }, error = function(e) {
    log_msg("warn", sprintf("NB fit failed for tile %s: %s", tile_name, e$message))
  })

  # 3. ZERO-INFLATED POISSON
  tryCatch({
    fit_zip <- pscl::zeroinfl(counts ~ 1 | 1, dist = "poisson")
    aic_values <- c(aic_values, AIC(fit_zip))
    bic_values <- c(bic_values, BIC(fit_zip))
    model_names <- c(model_names, "ZIP")
  }, error = function(e) {
    log_msg("warn", sprintf("ZIP fit failed for tile %s: %s", tile_name, e$message))
  })

  # 4. ZERO-INFLATED NEGATIVE BINOMIAL
  tryCatch({
    fit_zinb <- pscl::zeroinfl(counts ~ 1 | 1, dist = "negbin")
    aic_values <- c(aic_values, AIC(fit_zinb))
    bic_values <- c(bic_values, BIC(fit_zinb))
    model_names <- c(model_names, "ZINB")
  }, error = function(e) {
    log_msg("warn", sprintf("ZINB fit failed for tile %s: %s", tile_name, e$message))
  })

  # 5. HURDLE NEGATIVE BINOMIAL
  tryCatch({
    fit_hurdle <- pscl::hurdle(counts ~ 1 | 1, dist = "negbin")
    aic_values <- c(aic_values, AIC(fit_hurdle))
    bic_values <- c(bic_values, BIC(fit_hurdle))
    model_names <- c(model_names, "Hurdle-NB")
  }, error = function(e) {
    log_msg("warn", sprintf("Hurdle fit failed for tile %s: %s", tile_name, e$message))
  })

  # Add AIC/BIC results
  best_model_aic <- if (length(model_names) > 0) model_names[which.min(aic_values)] else NA

  # Add model-specific AIC columns (need to do this via named assignment in the list)
  if (length(model_names) > 0) {
    for (j in seq_along(model_names)) {
      col_name <- sprintf("AIC_%s", model_names[j])
      fit_row[[col_name]] <- aic_values[j]
    }
  }

  fit_row$best_model_AIC <- best_model_aic

  model_summary[[length(model_summary) + 1]] <- fit_row
}

# Combine results into a data frame
# Convert list of lists to data frame row by row
if (length(model_summary) > 0) {
  model_fits_df <- do.call(rbind, lapply(model_summary, as.data.frame, stringsAsFactors = FALSE))
  rownames(model_fits_df) <- NULL
} else {
  model_fits_df <- data.frame()
}

log_msg("step", sprintf("Completed model fits for %d sampled tiles", nrow(model_fits_df)))

# Write summary
csv_file <- file.path(output_dir, sprintf("%s_%dbp_tile_model_summary.csv", tissue, tilesize))
write.csv(model_fits_df, csv_file, row.names = FALSE)
log_msg("step", sprintf("Saved model fits summary: %s", csv_file))

# ============================================================================
# 7. Summary report with citations
# ============================================================================

log_msg("step", "Generating final report with citations...")

report_file <- file.path(output_dir, sprintf("%s_distribution_modeling_report.txt",
                                              tissue))

sink(report_file)

cat("================================================================================\n")
cat("SPATIAL ATAC-SEQ PER-TILE COUNT DISTRIBUTION ANALYSIS\n")
cat("================================================================================\n\n")

cat(sprintf("Tissue: %s\n", tissue))
cat(sprintf("Tile Size: %d bp\n", tilesize))
cat(sprintf("Analysis Date: %s\n\n", format(Sys.time(), "%F %T")))

cat("OVERVIEW\n")
cat("-" , rep("", 79), "\n\n")
cat("We characterize the statistical distribution of counts for genomic tiles in the\n")
cat("spatial ATAC-seq data. For each tile, we examine its count vector across all\n")
cat("spots and test whether it follows a Poisson, Negative Binomial, Zero-Inflated,\n")
cat("or Hurdle distribution.\n\n")

cat("MOTIVATION\n")
cat("-" , rep("", 79), "\n\n")
cat("In RNA-seq analysis, per-gene read counts across samples follow a Negative Binomial\n")
cat("(NB) distribution rather than Poisson, due to overdispersion from biological and\n")
cat("technical variance (Robinson & Smyth 2007, Bioinformatics 23:2881-2887; Anders &\n")
cat("Huber 2010, Genome Biol 11:R106; Love et al. 2014, Genome Biol 15:550). This\n")
cat("insight drives the differential-expression methods DESeq2 and edgeR. Here, we ask\n")
cat("the analogous question for spatial ATAC-seq: are tile counts also overdispersed?\n")
cat("Previous work shows NB is appropriate for ChIP/ATAC genomic windows (Lun & Smyth\n")
cat("2016, NAR 44:e45), but spatial ATAC with 1-10 nuclei per spot may introduce\n")
cat("additional zero-inflation from incomplete capture (Risso et al. 2018, Nat Commun\n")
cat("9:284; ZINB-WaVE model). Understanding the count distribution is critical for\n")
cat("choosing tile size (500bp vs 5000bp) and downstream analysis strategy.\n\n")

cat("KEY FINDINGS\n")
cat("-" , rep("", 79), "\n\n")
cat(sprintf("Total informative tiles (sum >= 10): %d\n", n_after))
cat(sprintf("Overdispersed tiles (D > 1.2, p < 0.05): %.1f%%\n", pct_over))
cat(sprintf("Poisson-consistent tiles (D in [0.8,1.2], p >= 0.05): %.1f%%\n", pct_poisson))
cat(sprintf("Underdispersed tiles (D < 0.8, p < 0.05): %.1f%%\n\n", pct_under))

cat(sprintf("Negative Binomial dispersion parameter (alpha, method-of-moments): %.6f\n",
            alpha_method_of_moments))
cat("  (Interpretation: For NB with E[X]=mu, Var[X]=mu + alpha*mu^2)\n\n")

cat(sprintf("Per-tile sample for model comparison: %d tiles (quantile-stratified)\n\n",
            nrow(model_fits_df)))

cat("MODEL COMPARISON ON SAMPLED TILES\n")
cat("-" , rep("", 79), "\n\n")

# Tally which models win
if (nrow(model_fits_df) > 0 && "best_model_AIC" %in% colnames(model_fits_df)) {
  model_wins <- table(model_fits_df$best_model_AIC)
  cat("Best-fit models (by AIC) across sampled tiles:\n")
  for (model in names(sort(model_wins, decreasing = TRUE))) {
    pct <- 100 * model_wins[model] / nrow(model_fits_df)
    cat(sprintf("  %s: %d tiles (%.1f%%)\n", model, model_wins[model], pct))
  }
  cat("\n")
}

cat("REFERENCES\n")
cat("-" , rep("", 79), "\n\n")

cat("1. Robinson MD, Smyth GK. (2007) Moderated statistical tests for assessing\n")
cat("   differences in tag abundance. Bioinformatics 23:2881-2887.\n")
cat("   → Foundational NB modeling of sequencing counts (edgeR).\n\n")

cat("2. Anders S, Huber W. (2010) Differential expression analysis for sequence\n")
cat("   count data. Genome Biol 11:R106.\n")
cat("   → DESeq's NB framework and mean-variance diagnostic approach.\n\n")

cat("3. Love MI, Huber W, Anders S. (2014) Moderated estimation of fold change\n")
cat("   and dispersion for RNA-seq data with DESeq2. Genome Biol 15:550.\n")
cat("   → plotDispEsts-style mean-variance plot, reused here for ATAC tiles.\n\n")

cat("4. Lun ATL, Smyth GK. (2016) csaw: a Bioconductor package for differential\n")
cat("   binding analysis of ChIP-seq data using sliding windows. Nucleic Acids\n")
cat("   Res 44:e45.\n")
cat("   → Direct precedent: applying RNA-seq-style NB modeling to ChIP/ATAC windows.\n\n")

cat("5. Buenrostro JD, et al. (2015) Single-cell chromatin accessibility reveals\n")
cat("   principles of regulatory variation. Nature 523:486-490.\n")
cat("   → Original scATAC-seq paper; basis for diploid sparsity argument.\n\n")

cat("6. Risso D, Perraudeau F, Gribkova S, Dudoit S, Vert JP. (2018) A general\n")
cat("   and flexible method for signal extraction from single-cell RNA-seq data.\n")
cat("   Nat Commun 9:284.\n")
cat("   → ZINB-WaVE: motivates zero-inflated NB for sparse single-cell data.\n\n")

cat("7. Svensson V. (2020) Droplet scRNA-seq is not zero-inflated. Nat Biotechnol 38:147-150.\n")
cat("   → Caveat: NB alone may explain apparent excess zeros in some droplet data.\n\n")

cat("8. Granja JM, et al. (2021) ArchR is a scalable software package for integrative\n")
cat("   single-cell chromatin accessibility analysis. Nat Genet 53:403-411.\n")
cat("   → Context: ArchR uses TF-IDF/LSI, not explicit generative count models.\n\n")

sink()

log_msg("step", sprintf("Saved report: %s", report_file))
log_msg("done", sprintf("Completed successfully: %s %dbp", tissue, tilesize))
