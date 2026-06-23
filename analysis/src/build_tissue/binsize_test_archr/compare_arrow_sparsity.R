#!/usr/bin/env Rscript
# ============================================================================
# compare_arrow_sparsity.R
#
# SAFE aggregator: reads pre-computed metrics from text files (created by
# create_arrow_variants.R), aggregates them into CSV, and generates visualizations.
#
# Does NOT touch arrow files, does NOT call ArchR — eliminates risk of
# accidentally corrupting previously-computed arrows.
#
# Usage:
#   Rscript compare_arrow_sparsity.R [tissue_name]
#
# If tissue_name is provided, only analyze that tissue
# If not provided, analyze all tissues
#
# Output:
#   1. Sparsity comparison table (CSV)
#   2. Sparsity plots (PDF)
#   3. Performance metrics summary (TXT)
#
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(tidyverse)
})

# Logging function
log_msg <- function(tag, msg) {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), tag, msg))
}

# Setup paths and configuration
project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
analysis_output_dir <- file.path(project_root, "analysis", "binsize_comparison")

dir.create(analysis_output_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("start", "===== Arrow Sparsity Comparison Analysis =====")

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
filter_tissue <- if (length(args) > 0) args[1] else NULL

# Define all tissues
all_tissues <- c("deepseq_488B", "deepseq_489", "lowseq_488B", "lowseq_489")
tilesizes <- c(500, 5000)
binarize_options <- c(FALSE, TRUE)

log_msg("info", sprintf("Filter tissue: %s", if (is.null(filter_tissue)) "all" else filter_tissue))

# Function to parse metrics from a text file
parse_metrics_file <- function(metrics_file) {
  tryCatch({
    if (!file.exists(metrics_file)) {
      return(NULL)
    }

    lines <- readLines(metrics_file)
    parsed <- list()

    for (line in lines) {
      if (nchar(line) == 0) next
      parts <- strsplit(line, ": ")[[1]]
      if (length(parts) == 2) {
        key <- trimws(parts[1])
        value <- trimws(parts[2])
        parsed[[key]] <- value
      }
    }

    # Convert numeric fields
    # Extract tilesize (remove " bp" suffix)
    tilesize_str <- sub(" bp", "", parsed[["Tilesize"]])
    # Extract sparsity (remove percentage in parentheses)
    sparsity_str <- sub(" \\(.*", "", parsed[["Sparsity"]])

    result <- list(
      tissue = parsed[["Tissue"]],
      tilesize = as.integer(tilesize_str),
      binarize = as.logical(parsed[["Binarize"]]),
      n_cells = as.integer(parsed[["Cells"]]),
      n_tiles = as.integer(parsed[["Tiles"]]),
      n_nnz = as.integer(parsed[["Non-zero elements"]]),
      density = as.numeric(parsed[["Density"]]),
      sparsity = as.numeric(sparsity_str),
      mean_cell_coverage = as.numeric(parsed[["Mean tiles per cell"]]),
      median_cell_coverage = as.numeric(parsed[["Median tiles per cell"]]),
      mean_tile_coverage = as.numeric(parsed[["Mean cells per tile"]]),
      median_tile_coverage = as.numeric(parsed[["Median cells per tile"]])
    )

    return(as.data.frame(result))
  }, error = function(e) {
    log_msg("warn", sprintf("Error parsing metrics file %s: %s", metrics_file, e$message))
    return(NULL)
  })
}

# Collect metrics for all combinations
log_msg("step", "Collecting sparsity metrics from text files...")

metrics_list <- list()
idx <- 0

tissues_to_process <- if (!is.null(filter_tissue)) filter_tissue else all_tissues

for (tissue in tissues_to_process) {
  for (ts in tilesizes) {
    for (bin in binarize_options) {
      metrics_file <- file.path(analysis_output_dir,
                                sprintf("%s_%dbp_binarize%s_metrics.txt", tissue, ts, bin))

      if (file.exists(metrics_file)) {
        metrics <- parse_metrics_file(metrics_file)
        if (!is.null(metrics)) {
          idx <- idx + 1
          metrics_list[[idx]] <- metrics
          log_msg("step", sprintf("Loaded: %s (%dbp, binarize=%s)", tissue, ts, bin))
        }
      } else {
        log_msg("warn", sprintf("Metrics file not found: %s", metrics_file))
      }
    }
  }
}

if (length(metrics_list) == 0) {
  log_msg("error", "No metrics files found. Run create_arrow_variants.R first.")
  stop("No data to analyze")
}

# Convert to data frame
metrics_df <- do.call(rbind, lapply(metrics_list, as.data.frame))
rownames(metrics_df) <- NULL
metrics_df <- as.data.table(metrics_df)

log_msg("step", sprintf("Collected metrics for %d combinations", nrow(metrics_df)))

# Save metrics table
output_table <- file.path(analysis_output_dir, "sparsity_comparison_table.csv")
fwrite(metrics_df, output_table)
log_msg("step", sprintf("Saved metrics table to: %s", output_table))

# Print summary
log_msg("step", "Sparsity Metrics Summary:")
print(metrics_df)

# Create comparison plots
log_msg("step", "Creating comparison plots...")

output_pdf <- file.path(analysis_output_dir, "sparsity_comparison_plots.pdf")

pdf(output_pdf, width = 12, height = 10)

# Plot 1: Sparsity comparison by tissue and tilesize
p1 <- ggplot(metrics_df, aes(x = factor(tilesize), y = sparsity * 100, fill = factor(binarize))) +
  geom_col(position = "dodge") +
  facet_wrap(~tissue) +
  labs(title = "Sparsity by Tile Size and Binarization",
       x = "Tile Size (bp)", y = "Sparsity (%)",
       fill = "Binarize") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)

# Plot 2: Density comparison
p2 <- ggplot(metrics_df, aes(x = factor(tilesize), y = density, fill = factor(binarize))) +
  geom_col(position = "dodge") +
  facet_wrap(~tissue) +
  labs(title = "Density by Tile Size and Binarization",
       x = "Tile Size (bp)", y = "Density",
       fill = "Binarize") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)

# Plot 3: Cell coverage comparison
p3 <- ggplot(metrics_df, aes(x = factor(tilesize), y = mean_cell_coverage, fill = factor(binarize))) +
  geom_col(position = "dodge") +
  facet_wrap(~tissue) +
  labs(title = "Mean Tiles per Cell by Tile Size and Binarization",
       x = "Tile Size (bp)", y = "Mean Tiles per Cell",
       fill = "Binarize") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p3)

# Plot 4: Comparison across tissues
p4 <- ggplot(metrics_df, aes(x = tissue, y = sparsity * 100,
                             fill = interaction(factor(tilesize), factor(binarize)))) +
  geom_col(position = "dodge") +
  labs(title = "Sparsity Comparison Across All Conditions",
       x = "Tissue", y = "Sparsity (%)",
       fill = "Tilesize×Binarize") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

print(p4)

dev.off()

log_msg("step", sprintf("Saved plots to: %s", output_pdf))

# Generate recommendations
log_msg("step", "Generating performance recommendations...")

recommendations <- list()

# Find best (highest density / lowest sparsity) for each tissue
for (tissue in unique(metrics_df$tissue)) {
  tissue_data <- metrics_df[tissue == tissue]
  best_idx <- which.min(tissue_data$sparsity)
  best_row <- tissue_data[best_idx]

  recommendations[[tissue]] <- sprintf(
    "%s: Best density=%.6f with %d bp tiles, binarize=%s (mean tiles/cell=%.2f)",
    tissue, best_row$density, best_row$tilesize, best_row$binarize, best_row$mean_cell_coverage
  )
}

# Find overall best
overall_best_idx <- which.min(metrics_df$sparsity)
overall_best <- metrics_df[overall_best_idx]

output_summary <- file.path(analysis_output_dir, "sparsity_recommendations.txt")
sink(output_summary)

cat("====================================\n")
cat("Arrow Sparsity Comparison Summary\n")
cat("====================================\n\n")

cat("Generated:", format(Sys.time(), "%F %T"), "\n\n")

cat("Overall Best Density Configuration:\n")
cat(sprintf("  Tissue: %s\n", overall_best$tissue))
cat(sprintf("  Tile Size: %d bp\n", overall_best$tilesize))
cat(sprintf("  Binarize: %s\n", overall_best$binarize))
cat(sprintf("  Density: %.6f\n", overall_best$density))
cat(sprintf("  Sparsity: %.4f (%.2f%%)\n", overall_best$sparsity, overall_best$sparsity * 100))
cat(sprintf("  Mean tiles per cell: %.2f\n", overall_best$mean_cell_coverage))
cat(sprintf("  Median tiles per cell: %.1f\n\n", overall_best$median_cell_coverage))

cat("Per-Tissue Recommendations:\n")
for (tissue in unique(metrics_df$tissue)) {
  cat(sprintf("  %s\n", recommendations[[tissue]]))
}

cat("\nInterpretation Guide:\n")
cat("  - Lower sparsity = more non-zero elements = more information retained\n")
cat("  - Higher density = more non-zero values, potentially better for downstream analysis\n")
cat("  - Tiles per cell = information density at the single-barcode level\n")
cat("    * For spatial ATAC with ~1-10 cells/spot, target mean >= 100 tiles/cell for robust CNV calling\n")
cat("    * If mean < 100, consider smaller tile size (higher resolution) to capture more signal\n\n")

cat("Full Metrics Table:\n")
sink()

metrics_wide <- metrics_df %>%
  select(tissue, tilesize, binarize, density, sparsity, mean_cell_coverage) %>%
  pivot_wider(names_from = c(tilesize, binarize),
              values_from = c(density, sparsity, mean_cell_coverage),
              names_sep = "_")

print(metrics_wide, file = output_summary, append = TRUE)

log_msg("step", sprintf("Saved recommendations to: %s", output_summary))
log_msg("done", "Arrow sparsity comparison completed successfully")

# Print to console as well
cat("\n\n===== SUMMARY =====\n")
print(readLines(output_summary))
