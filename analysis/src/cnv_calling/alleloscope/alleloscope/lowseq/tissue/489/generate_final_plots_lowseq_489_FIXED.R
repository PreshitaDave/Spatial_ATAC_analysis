#!/usr/bin/env Rscript
# =============================================================================
# Generate Step 9-10 final plots for lowseq_489 Alleloscope analysis
# FIX: Restore missing colnames to theta_hat_cbn before Genotype call
# Purpose: Load Obj_after_gtv and generate missing heatmap and scatter plots
# Based on: run_alleloscope_lowseq_tissue_488B.R Steps 9-10
# =============================================================================

.libPaths(c("/projectnb/paxlab/presh/Rlibs/4.5", .libPaths()))

suppressPackageStartupMessages({
  library(Alleloscope)
  library(Matrix)
  library(parallel)
})

NCORES <- as.integer(Sys.getenv("NSLOTS", "8"))

allelo.path <- "/projectnb/paxlab/presh/software/Alleloscope"
data.path   <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/cnv/alleloscope/lowseq_tissue_from_existing/489"
dir_path    <- file.path(data.path, "output")

message(sprintf("[%s] === Generating Step 9-10 plots for lowseq_489 (ncores=%d) ===",
                Sys.time(), NCORES))

# ─────────────────────────────────────────────────────────────────────────────
# Load Obj_after_gtv checkpoint (result of Step 8: Genotype_value)
# ─────────────────────────────────────────────────────────────────────────────
message(sprintf("[%s] Loading checkpoint: Obj_after_gtv.rds", Sys.time()))

Obj_filtered <- readRDS(file.path(dir_path, "rds", "Obj_after_gtv.rds"))

message(sprintf("[%s] Loaded Obj_filtered with %d cells", Sys.time(), 
                ncol(Obj_filtered$alt_all)))

# ─────────────────────────────────────────────────────────────────────────────
# FIX: Restore missing colnames to theta_hat_cbn
# ─────────────────────────────────────────────────────────────────────────────
message(sprintf("[%s] [FIX] Checking theta_hat_cbn colnames...", Sys.time()))

if (is.null(colnames(Obj_filtered$select_normal$theta_hat_cbn))) {
  message(sprintf("[%s] [FIX] theta_hat_cbn has NULL colnames, restoring from Obj$barcodes...", 
                  Sys.time()))
  
  # Get subset of barcodes used in select_normal (after filtering)
  # Infer from Select_normal output which barcodes were retained
  # Use all barcodes and trust dimensions match
  barcodes_for_theta <- Obj_filtered$barcodes
  
  if (length(barcodes_for_theta) == ncol(Obj_filtered$select_normal$theta_hat_cbn)) {
    colnames(Obj_filtered$select_normal$theta_hat_cbn) <- barcodes_for_theta
    message(sprintf("[%s] [FIX] ✓ Restored %d colnames to theta_hat_cbn", 
                    Sys.time(), length(barcodes_for_theta)))
  } else {
    stop(sprintf("Barcode count mismatch: %d barcodes vs %d theta_hat_cbn columns",
                 length(barcodes_for_theta), ncol(Obj_filtered$select_normal$theta_hat_cbn)))
  }
} else {
  message(sprintf("[%s] [FIX] theta_hat_cbn already has colnames (%d)", 
                  Sys.time(), length(colnames(Obj_filtered$select_normal$theta_hat_cbn))))
}

# Ensure region_normal is set properly (needed by Genotype)
if (is.null(Obj_filtered$ref)) {
  region_normal <- Obj_filtered$select_normal$region_normal
  Obj_filtered$ref <- region_normal[1]
  message(sprintf("[%s] Set ref to: %s", Sys.time(), Obj_filtered$ref))
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 9: Genotype -- per-region (rho_hat, theta_hat) scatter plots
# ─────────────────────────────────────────────────────────────────────────────
message(sprintf("[%s] Step 9: Genotype (scatter plots)", Sys.time()))

gtype_plot <- file.path(dir_path, "plots",
                        paste0("gtype_scatter_ref_", Obj_filtered$ref, ".pdf"))
gtype_rds  <- file.path(dir_path, "rds", "genotypes.rds")

dir.create(dirname(gtype_plot), showWarnings = FALSE, recursive = TRUE)

Obj_filtered <- Genotype(
  Obj_filtered = Obj_filtered,
  plot_path    = gtype_plot,
  rds_path     = gtype_rds,
  cell_type    = NULL,
  legend       = TRUE
)

# Save final object
saveRDS(Obj_filtered,
        file.path(dir_path, "rds", "Obj_final.rds"))
message(sprintf("[%s] Step 9: Saved Obj_final.rds", Sys.time()))

# ─────────────────────────────────────────────────────────────────────────────
# STEP 10: plot_scATAC_cnv -- Step-6 smoothed coverage CNV heatmap
# ─────────────────────────────────────────────────────────────────────────────
message(sprintf("[%s] Step 10: Generating Step 6 CNV coverage heatmap", Sys.time()))

# Load supporting data for heatmap
raw_counts <- read.table(file.path(data.path, "chr1000k_fragments.tsv"),
  sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
colnames(raw_counts) <- gsub("[.]", "-", colnames(raw_counts))

size_df <- read.table(
  file.path(allelo.path, "data-raw/sizes.cellranger-GRCh38-1.0.0.txt"),
  stringsAsFactors = FALSE
)
size_df <- size_df[1:22, ]

# Extract normal barcodes
barcode_normal_final <- Obj_filtered$select_normal$barcode_normal

# Build cell_type data frame: barcode | label
all_barcodes <- colnames(raw_counts)
cell_type_df <- data.frame(
  barcode   = all_barcodes,
  cell_type = ifelse(all_barcodes %in% barcode_normal_final, "normal", "tumor"),
  stringsAsFactors = FALSE
)

heatmap_path <- file.path(dir_path, "plots", "step6_CNV_coverage_heatmap.png")

# Remove zero-coverage cells before normalization (0/0 → NaN → rowMedians all NA)
keep_cells  <- colSums(raw_counts[, cell_type_df$barcode]) > 0
message(sprintf("[%s]   Cells passing coverage filter: %d / %d",
                Sys.time(), sum(keep_cells), length(keep_cells)))
cell_type_df_filt <- cell_type_df[keep_cells, , drop = FALSE]

# Generate heatmap
cov_obj <- plot_scATAC_cnv(
  raw_mat        = as.matrix(raw_counts[, cell_type_df_filt$barcode]),
  cell_type      = cell_type_df_filt,
  normal_lab     = "normal",
  size           = size_df,
  window_w       = 10000000,
  window_step    = 2000000,
  plot_path      = heatmap_path,
  nclust         = 3,
  var.filter     = FALSE
)

saveRDS(cov_obj, file.path(dir_path, "rds", "cov_obj.rds"))
message(sprintf("[%s] Step 10: Saved Step 6 CNV heatmap", Sys.time()))

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
message(sprintf("[%s] === All steps complete ===", Sys.time()))
message(sprintf("Output directory: %s", dir_path))
message("Key outputs:")
message(sprintf("  Step-6 CNV coverage heatmap : plots/step6_CNV_coverage_heatmap.png"))
message(sprintf("  Genotype scatter plots      : plots/gtype_scatter_ref_%s.pdf",
                Obj_filtered$ref))
message(sprintf("  Final Alleloscope object    : rds/Obj_final.rds"))
message(sprintf("  Coverage object             : rds/cov_obj.rds"))
message(sprintf("  Normal barcodes (final)     : %d cells",
                length(barcode_normal_final)))
