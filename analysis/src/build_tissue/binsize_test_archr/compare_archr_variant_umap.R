#!/usr/bin/env Rscript
# ============================================================================
# compare_archr_variant_umap.R
#
# Safe aggregator (no ArchR, no arrow file access) that reads per-variant
# UMAP embedding CSVs and generates comparison plots showing cluster
# separation across tile-size/binarize variants.
#
# Usage:
#   Rscript compare_archr_variant_umap.R [tissue_filter]
#
# If tissue_filter is provided (e.g., "lowseq_489"), only analyze that tissue
# If not provided, analyze all tissues that have completed runs
#
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(patchwork)
})

log_msg <- function(tag, msg) {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), tag, msg))
}

# Setup paths
project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
comparison_dir <- file.path(project_root, "analysis/binsize_comparison")
umap_scores_dir <- file.path(comparison_dir, "umap_scores")

log_msg("start", "===== Comparing ArchR UMAP variants =====")

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
filter_tissue <- if (length(args) > 0) args[1] else NULL

log_msg("info", sprintf("Tissue filter: %s",
                        if (is.null(filter_tissue)) "all" else filter_tissue))

# Define grid
all_tissues <- c("deepseq_488B", "deepseq_489", "lowseq_488B", "lowseq_489")
tilesizes <- c(500, 5000)
binarize_options <- c(FALSE, TRUE)

# Load all available CSVs
log_msg("step", "Loading UMAP embedding CSVs...")

all_data <- list()

tissues_to_process <- if (!is.null(filter_tissue)) filter_tissue else all_tissues

for (tissue in tissues_to_process) {
  for (ts in tilesizes) {
    for (bin in binarize_options) {
      csv_file <- file.path(umap_scores_dir,
                           sprintf("%s_%dbp_binarize%s_umap_embeddings.csv",
                                   tissue, ts, bin))

      if (file.exists(csv_file)) {
        df <- read.csv(csv_file, stringsAsFactors = FALSE)
        df$variant <- sprintf("%dbp_%s", ts, bin)
        df$tissue <- tissue
        all_data[[length(all_data) + 1]] <- df
        log_msg("step", sprintf("Loaded: %s (%dbp, %s)", tissue, ts, bin))
      } else {
        log_msg("warn", sprintf("CSV not found: %s", csv_file))
      }
    }
  }
}

if (length(all_data) == 0) {
  log_msg("error", "No UMAP CSVs found. Run build_archr_variant_project.R first.")
  stop("No data to compare")
}

# Combine all data
combined_df <- do.call(rbind, all_data)
rownames(combined_df) <- NULL

log_msg("step", sprintf("Combined data from %d variants, %d total cells",
                        length(all_data), nrow(combined_df)))

# Identify marker genes (columns that aren't metadata)
metadata_cols <- c("cellID", "Clusters", "UMAP_1", "UMAP_2", "variant", "tissue")
marker_cols <- setdiff(colnames(combined_df), metadata_cols)

log_msg("step", sprintf("Found %d marker genes in combined data", length(marker_cols)))

# Generate comparison plots
log_msg("step", "Generating comparison PDF...")

output_pdf <- file.path(comparison_dir, "archr_umap_genescore_comparison.pdf")

pdf(output_pdf, width = 16, height = 12, onefile = TRUE)

# Page 1: UMAP colored by Clusters, faceted by variant
p_clusters <- ggplot(combined_df, aes(x = UMAP_1, y = UMAP_2, color = Clusters)) +
  geom_point(size = 0.8, alpha = 0.7) +
  facet_wrap(~variant, nrow = 2) +
  labs(title = sprintf("UMAP Clustering Comparison: %s",
                       if (!is.null(filter_tissue)) filter_tissue else "All tissues"),
       x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

print(p_clusters)

# Pages 2+: One gene per page, faceted by variant
genes_per_page <- 4
n_pages <- ceiling(length(marker_cols) / genes_per_page)

for (page in 1:n_pages) {
  start_idx <- (page - 1) * genes_per_page + 1
  end_idx <- min(page * genes_per_page, length(marker_cols))
  page_genes <- marker_cols[start_idx:end_idx]

  for (gene in page_genes) {
    # Skip if all NA
    if (all(is.na(combined_df[[gene]]))) {
      next
    }

    p <- ggplot(combined_df, aes(x = UMAP_1, y = UMAP_2, color = !!sym(gene))) +
      geom_point(size = 0.8, alpha = 0.7) +
      facet_wrap(~variant, nrow = 2) +
      scale_color_viridis_c(name = gene, option = "plasma") +
      labs(title = sprintf("Gene Score: %s", gene),
           x = "UMAP 1", y = "UMAP 2") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = "right"
      )

    print(p)
  }
}

dev.off()

log_msg("step", sprintf("Saved comparison PDF: %s", output_pdf))

# Summary report
log_msg("step", "Generating summary report...")

output_summary <- file.path(comparison_dir, "umap_comparison_summary.txt")

sink(output_summary)

cat("====================================\n")
cat("ArchR UMAP Variant Comparison Summary\n")
cat("====================================\n\n")

cat("Generated:", format(Sys.time(), "%F %T"), "\n\n")

cat("Variants analyzed:\n")
for (variant in unique(combined_df$variant)) {
  n_cells <- nrow(combined_df[combined_df$variant == variant, ])
  n_clusters <- length(unique(combined_df[combined_df$variant == variant, "Clusters"]))
  cat(sprintf("  %s: %d cells, %d clusters\n", variant, n_cells, n_clusters))
}

cat("\nMarker genes evaluated:\n")
cat(sprintf("  %s\n", paste(marker_cols, collapse = ", ")))

cat("\nComparison plots:\n")
cat("  Page 1: UMAP colored by Clusters (visual assessment of cluster separation)\n")
cat("  Pages 2+: One marker gene per page, faceted by variant\n")
cat("           Check for consistency of marker expression patterns across variants\n")

cat("\nKey observations to check:\n")
cat("  1. Cluster separation: Are clusters similarly separated across all variants?\n")
cat("  2. Marker specificity: Do markers show expected tissue/cell-type patterns?\n")
cat("  3. Consistency: Does 500bp vs 5000bp or binarize vs not-binarize show major differences?\n")
cat("  4. Artifacts: Any unusual/noisy embedding patterns in particular variant(s)?\n")

sink()

log_msg("step", sprintf("Saved summary report: %s", output_summary))

log_msg("done", "ArchR UMAP comparison completed successfully")
