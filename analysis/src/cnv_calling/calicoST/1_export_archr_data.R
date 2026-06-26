#!/usr/bin/env Rscript
# ============================================================================
# 1_export_archr_data.R
#
# Export TileMatrix (raw tile counts) and spatial coordinates from an ArchR
# project for downstream CalicoST input preparation (script 2).
#
# Usage:
#   Rscript 1_export_archr_data.R <tissue> [bin_size]
#   Example: Rscript 1_export_archr_data.R lowseq_489 5000
#
# Outputs (written to Data/04_analysis/cnv/calicoST/<tissue>/intermediate/):
#   archr_tilematrix.mtx  — sparse tile count matrix (tiles × cells), Market Exchange format
#   tile_ranges.csv       — genomic coordinates per tile row (chr, start, end)
#   barcodes.csv          — cell barcodes in matrix column order
#   spatial_coords.csv    — per-cell barcode, x_spatial, y_spatial, array_row, array_col
# ============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(Matrix)
})

log_msg <- function(tag, msg) cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), tag, msg))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 1_export_archr_data.R <tissue> [bin_size]")

tissue   <- args[1]
bin_size <- if (length(args) >= 2) as.integer(args[2]) else 5000L

log_msg("start", sprintf("Exporting ArchR data for tissue: %s (bin_size=%d)", tissue, bin_size))

project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
proj_dir <- file.path(project_root, "analysis/binsize_comparison/archr_projects",
                      sprintf("%s_%dbp_binarizeFALSE", tissue, bin_size))
spatial_coord_file <- file.path(project_root, "Data/01_inputs/spatial/tissue_positions_list.csv")
out_dir <- file.path(project_root, "Data/04_analysis/cnv/calicoST", tissue, "intermediate")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(proj_dir)) stop(sprintf("ArchR project not found: %s", proj_dir))

set.seed(42)
addArchRGenome("hg38")
addArchRThreads(threads = 4)

log_msg("step", sprintf("Loading ArchR project: %s", proj_dir))
proj <- loadArchRProject(path = proj_dir, force = TRUE)
log_msg("step", sprintf("Project loaded: %d cells", nrow(proj@cellColData)))

# ============================================================================
# Attach spatial coordinates (order-preserving — see compare_spatial_smoothing_methods.R)
# ============================================================================
tissue_name_map <- list(
  "lowseq_489"  = "Lowseq_489",
  "lowseq_488B" = "Lowseq_488B",
  "deepseq_488B" = "Deepseq_488B",
  "deepseq_489"  = "Deepseq_489"
)
sample_prefix <- tissue_name_map[[tissue]]
if (is.null(sample_prefix)) stop(sprintf("Unknown tissue: %s", tissue))

tissue_locs <- read.csv(spatial_coord_file)
tissue_locs <- tissue_locs[tissue_locs$in_tissue == 1, ]
tissue_locs$cellName <- paste0(sample_prefix, "#", tissue_locs$barcode, "-1")

# Use match() to preserve ArchR cell order (never merge/sort)
match_idx <- match(rownames(proj@cellColData), tissue_locs$cellName)
proj@cellColData$x_spatial  <- tissue_locs$x_spatial[match_idx]
proj@cellColData$y_spatial  <- tissue_locs$y_spatial[match_idx]
proj@cellColData$array_row  <- tissue_locs$array_row[match_idx]
proj@cellColData$array_col  <- tissue_locs$array_col[match_idx]

n_matched <- sum(!is.na(proj@cellColData$x_spatial))
log_msg("step", sprintf("Spatial coords matched: %d / %d cells", n_matched, nrow(proj@cellColData)))

# ============================================================================
# Export TileMatrix
# ============================================================================
log_msg("step", "Extracting TileMatrix...")

tile_se <- getMatrixFromProject(proj, useMatrix = "TileMatrix", binarize = FALSE)

# tile_se is a SummarizedExperiment: rows = tiles, cols = cells
tile_mat <- assay(tile_se, "TileMatrix")  # sparse Matrix (tiles × cells)

# Genomic ranges from rowData (rowRanges is NULL for TileMatrix in this ArchR version)
# rowData has columns: seqnames, idx, start; end = start + bin_size
rd <- as.data.frame(rowData(tile_se))
tile_ranges_df <- data.frame(
  chr   = as.character(rd$seqnames),
  start = as.integer(rd$start),
  end   = as.integer(rd$start) + bin_size
)

# Barcodes in column order (same as ArchR cellColData row order)
barcodes_vec <- colnames(tile_mat)

log_msg("step", sprintf("TileMatrix: %d tiles × %d cells", nrow(tile_mat), ncol(tile_mat)))

# ============================================================================
# Export spatial coordinates per cell
# ============================================================================
meta_df <- as.data.frame(proj@cellColData)

# Bare barcode (strip "SamplePrefix#" and "-1" suffix) for joining with numbat
bare_barcode <- sub("^[^#]+#", "", rownames(meta_df))
bare_barcode <- sub("-1$", "", bare_barcode)

spatial_df <- data.frame(
  archr_barcode = rownames(meta_df),
  bare_barcode  = bare_barcode,
  x_spatial     = meta_df$x_spatial,
  y_spatial     = meta_df$y_spatial,
  array_row     = meta_df$array_row,
  array_col     = meta_df$array_col
)

# ============================================================================
# Write outputs
# ============================================================================
log_msg("step", "Writing outputs...")

writeMM(tile_mat, file.path(out_dir, "archr_tilematrix.mtx"))
write.csv(tile_ranges_df, file.path(out_dir, "tile_ranges.csv"), row.names = FALSE)
write.csv(data.frame(barcode = barcodes_vec), file.path(out_dir, "barcodes.csv"), row.names = FALSE)
write.csv(spatial_df, file.path(out_dir, "spatial_coords.csv"), row.names = FALSE)

log_msg("done", sprintf("Outputs written to: %s", out_dir))
log_msg("done", sprintf("  archr_tilematrix.mtx : %d tiles x %d cells", nrow(tile_mat), ncol(tile_mat)))
log_msg("done", sprintf("  tile_ranges.csv      : %d rows", nrow(tile_ranges_df)))
log_msg("done", sprintf("  barcodes.csv         : %d barcodes", length(barcodes_vec)))
log_msg("done", sprintf("  spatial_coords.csv   : %d rows (%d with spatial coords)",
                        nrow(spatial_df), sum(!is.na(spatial_df$x_spatial))))
