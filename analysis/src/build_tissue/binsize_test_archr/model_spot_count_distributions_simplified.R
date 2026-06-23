#!/usr/bin/env Rscript
# Simplified version: just aggregate statistics (fast), no per-tile fitting

suppressPackageStartupMessages({
  library(ArchR)
  library(Matrix)
  library(ggplot2)
})

log_msg <- function(tag, msg) {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), tag, msg))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  log_msg("error", "Usage: Rscript model_spot_count_distributions.R <tissue> <tilesize>")
  stop("Missing arguments")
}

tissue <- args[1]
tilesize <- as.integer(args[2])
log_msg("start", sprintf("===== Modeling per-tile count distributions: %s, %dbp =====", tissue, tilesize))

set.seed(42)
addArchRGenome("hg38")

project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
arrow_file <- file.path(project_root, "Data/01_inputs/arrow/arrow_not_binarize",
                        sprintf("%s_%dbp.arrow", tissue, tilesize))

output_dir <- file.path(project_root, "analysis/binsize_comparison/distribution_modeling")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("step", "Loading non-binarized TileMatrix from arrow file...")

tile_matrix <- ArchR::getMatrixFromArrow(ArrowFile = arrow_file, useMatrix = "TileMatrix",
                                         binarize = FALSE)
tile_mat <- assay(tile_matrix)
log_msg("step", sprintf("Loaded TileMatrix: %d tiles × %d spots", nrow(tile_mat), ncol(tile_mat)))

# Filter to informative tiles
row_sums <- Matrix::rowSums(tile_mat)
informative_tiles <- row_sums >= 10
tile_mat_filt <- tile_mat[informative_tiles, ]

n_before <- nrow(tile_mat)
n_after <- nrow(tile_mat_filt)
log_msg("step", sprintf("Filtered from %d tiles to %d tiles (removed %d with sum < 10)",
                        n_before, n_after, n_before - n_after))

# Compute mean, variance, dispersion
log_msg("step", "Computing per-tile dispersion statistics...")
tile_means <- Matrix::rowMeans(tile_mat_filt)
tile_mat_sq <- tile_mat_filt^2
tile_means_sq <- Matrix::rowMeans(tile_mat_sq)
tile_vars <- tile_means_sq - tile_means^2
dispersion_index <- tile_vars / tile_means

# Chi-squared test for overdispersion
n_spots <- ncol(tile_mat_filt)
test_stat <- (n_spots - 1) * dispersion_index
p_values <- pchisq(test_stat, df = n_spots - 1, lower.tail = FALSE)

significantly_overdispersed <- p_values < 0.05 & dispersion_index > 1.2
poisson_consistent <- p_values > 0.05 & dispersion_index >= 0.8 & dispersion_index <= 1.2
significantly_underdispersed <- p_values < 0.05 & dispersion_index < 0.8

pct_over <- 100 * sum(significantly_overdispersed) / n_after
pct_under <- 100 * sum(significantly_underdispersed) / n_after
pct_poisson <- 100 * sum(poisson_consistent) / n_after

log_msg("step", sprintf("Dispersion: %.1f%% overdispersed, %.1f%% Poisson-consistent, %.1f%% underdispersed",
                        pct_over, pct_poisson, pct_under))

# Fit NB trend
valid_idx <- tile_means > 0 & tile_vars > 0
mu_valid <- tile_means[valid_idx]
var_valid <- tile_vars[valid_idx]
alpha_method_of_moments <- mean((var_valid - mu_valid) / (mu_valid^2))

# Mean-variance plot
log_msg("step", "Generating mean-variance diagnostic plot...")

plot_df <- data.frame(
  mean = tile_means,
  variance = tile_vars,
  dispersion_class = ifelse(
    poisson_consistent, "Poisson",
    ifelse(significantly_overdispersed, "Overdispersed", "Underdispersed")
  )
)

pdf(file.path(output_dir, sprintf("%s_%dbp_meanvar_plot.pdf", tissue, tilesize)),
    width = 12, height = 10)

p1 <- ggplot(plot_df, aes(x = mean, y = variance, color = dispersion_class)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed", linewidth = 1) +
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
    caption = sprintf("Alpha (NB dispersion) = %.4f\nOverdispersed: %.1f%%, Poisson: %.1f%%, Underdispersed: %.1f%%",
                      alpha_method_of_moments, pct_over, pct_poisson, pct_under)
  ) +
  theme_minimal() + theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold"),
                           plot.caption = element_text(hjust = 0.5, size = 9))
print(p1)

p2 <- ggplot(plot_df, aes(x = dispersion_index)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "black") +
  geom_vline(xintercept = 1, color = "blue", linetype = "dashed", linewidth = 1.5) +
  labs(title = sprintf("%s, %dbp: Dispersion Index Distribution", tissue, tilesize),
       x = "Dispersion Index (D = Variance / Mean)", y = "Number of Tiles") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(p2)

dev.off()
log_msg("step", sprintf("Saved mean-variance plot: %s_%dbp_meanvar_plot.pdf", tissue, tilesize))

# Generate summary report
report_file <- file.path(output_dir, sprintf("%s_distribution_modeling_report.txt", tissue))
sink(report_file)

cat("================================================================================\n")
cat("SPATIAL ATAC-SEQ PER-TILE COUNT DISTRIBUTION ANALYSIS\n")
cat("================================================================================\n\n")

cat(sprintf("Tissue: %s | Tile Size: %d bp | Date: %s\n\n", tissue, tilesize, format(Sys.time(), "%F %T")))

cat("KEY FINDINGS\n")
cat("-" , rep("", 79), "\n\n")
cat(sprintf("Total informative tiles (sum >= 10): %d / %d (%.1f%%)\n", 
            n_after, n_before, 100*n_after/n_before))
cat(sprintf("\nDispersion Classification (Chi-squared test vs Poisson null, p=0.05):\n"))
cat(sprintf("  Overdispersed (D > 1.2): %.1f%% of tiles\n", pct_over))
cat(sprintf("  Poisson-consistent (0.8 <= D <= 1.2): %.1f%% of tiles\n", pct_poisson))
cat(sprintf("  Underdispersed (D < 0.8): %.1f%% of tiles\n\n", pct_under))

cat(sprintf("Negative Binomial Dispersion Parameter (alpha, method-of-moments): %.6f\n", alpha_method_of_moments))
cat("  (Interpretation: For NB, Var[X] = mu + alpha*mu^2)\n\n")

cat("INTERPRETATION\n")
cat("-" , rep("", 79), "\n\n")
if (pct_over > 90) {
  cat("STRONG EVIDENCE FOR NEGATIVE BINOMIAL DISTRIBUTION\n\n")
  cat("Nearly all tiles show statistically significant overdispersion (variance > mean),\n")
  cat("inconsistent with Poisson and highly consistent with Negative Binomial. This\n")
  cat("validates the RNA-seq literature finding (Robinson & Smyth 2007, edgeR;\n")
  cat("Anders & Huber 2010, DESeq; Love et al. 2014, DESeq2) that sequencing count\n")
  cat("data exhibit overdispersion from biological and technical variance.\n\n")
} else if (pct_over > 50) {
  cat("MODERATE EVIDENCE FOR NEGATIVE BINOMIAL DISTRIBUTION\n\n")
  cat("Majority of tiles show overdispersion; mix suggests both NB-like and\n")
  cat("Poisson-like tiles depending on abundance or coverage level.\n\n")
} else {
  cat("WEAK / MIXED EVIDENCE\n\n")
  cat("Count distribution appears mixed or Poisson-like; proceed with caution\n")
  cat("assuming overdispersion in downstream analysis.\n\n")
}

cat("REFERENCES\n")
cat("-" , rep("", 79), "\n\n")
cat("1. Robinson MD, Smyth GK. (2007) Bioinformatics 23:2881-2887.\n")
cat("   Foundational NB modeling of sequencing counts (edgeR).\n\n")
cat("2. Anders S, Huber W. (2010) Genome Biol 11:R106.\n")
cat("   DESeq NB framework and mean-variance diagnostic approach.\n\n")
cat("3. Love MI, Huber W, Anders S. (2014) Genome Biol 15:550.\n")
cat("   DESeq2 and plotDispEsts-style mean-variance diagnostic.\n\n")
cat("4. Lun ATL, Smyth GK. (2016) NAR 44:e45.\n")
cat("   csaw: applying NB window-binning to ChIP/ATAC genomic windows.\n\n")
cat("5. Buenrostro JD, et al. (2015) Nature 523:486-490.\n")
cat("   Original scATAC-seq paper; diploid sparsity basis for NB interpretation.\n\n")

sink()
log_msg("step", sprintf("Saved report: %s", report_file))
log_msg("done", sprintf("Completed: %s %dbp", tissue, tilesize))

