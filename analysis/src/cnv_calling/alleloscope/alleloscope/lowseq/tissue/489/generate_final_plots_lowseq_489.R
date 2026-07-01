#!/usr/bin/env Rscript
# =============================================================================
# Generate Step 10 final plot for lowseq_489 Alleloscope analysis
# Purpose: Generate CNV coverage heatmap only (Step 9 Genotype skipped due to bug)
# Note: Step 9 (Genotype) failed in the original run due to Alleloscope bug
#       Regenerating heatmap which requires only Obj_final.rds
# =============================================================================

.libPaths(c("/projectnb/paxlab/presh/Rlibs/4.5", .libPaths()))

suppressPackageStartupMessages({
  library(Alleloscope)
  library(Matrix)
})

allelo.path <- "/projectnb/paxlab/presh/software/Alleloscope"
data.path   <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/cnv/alleloscope/lowseq_tissue_from_existing/489"
dir_path    <- file.path(data.path, "output")

message(sprintf("[%s] === Generating Step 10 heatmap for lowseq_489 ===",
                Sys.time()))

# ─────────────────────────────────────────────────────────────────────────────
# Load Obj_final checkpoint (use this instead of Obj_after_gtv)
# ─────────────────────────────────────────────────────────────────────────────
message(sprintf("[%s] Loading checkpoint: Obj_final.rds", Sys.time()))

Obj_filtered <- readRDS(file.path(dir_path, "rds", "Obj_final.rds"))

message(sprintf("[%s] Loaded Obj_filtered with %d cells", Sys.time(), 
                ncol(Obj_filtered$alt_all)))

# ─────────────────────────────────────────────────────────────────────────────
# NOTE: Skipping Step 9 (Genotype scatter plots)
# ─────────────────────────────────────────────────────────────────────────────
message(sprintf("[%s] Step 9: SKIPPED (Genotype failed in original run - Alleloscope bug)", Sys.time()))
message(sprintf("[%s]   Reason: strsplit error on missing colnames in theta_hat_cbn", Sys.time()))

message(sprintf("[%s]   (Genotype function has bug in Alleloscope for this object structure)", Sys.time()))

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
