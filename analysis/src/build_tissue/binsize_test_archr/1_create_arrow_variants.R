#!/usr/bin/env Rscript
# ============================================================================
# create_arrow_variants_v2.R
#
# Fixed version of create_arrow_variants.R: the original never applied the
# no_edge_effect barcode whitelist (validBarcodes was defined but unused) and
# used a stricter arrow-level QC (minTSS=3, minFrags=1000) than the main
# pipeline (minTSS=2, minFrags=100), which combined to drop the deepseq_488B
# 5000bp non-binarized cell count from ~11k to 7,841.
#
# This version creates the arrow file with the SAME permissive QC as
# 0_create_archr_qc_cluster.R (minTSS=2, minFrags=100, maxFrags=Inf) so the
# real filtering happens downstream via the no_edge_effect barcode list +
# high-nFrags doublet removal, in build_archr_variant_project_v2.R.
#
# Usage:
#   Rscript create_arrow_variants_v2.R <tissue_name> <tilesize> <binarize>
#   Example: Rscript create_arrow_variants_v2.R deepseq_488B 5000 FALSE
#
# Output:
#   /Data/01_inputs/arrow/arrow_not_binarize_v2/{tissue}_{tilesize}bp.arrow (or arrow_binarize_v2/)
# ============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(data.table)
})

log_msg <- function(tag, msg) {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), tag, msg))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript create_arrow_variants_v2.R <tissue_name> <tilesize> <binarize>")
}

tissue_name <- args[1]
tilesize <- as.integer(args[2])
binarize_arg <- args[3]
binarize <- if (tolower(binarize_arg) %in% c("true", "t", "1")) TRUE else if (tolower(binarize_arg) %in% c("false", "f", "0")) FALSE else as.logical(binarize_arg)

log_msg("info", sprintf("Parameters: tissue=%s, tilesize=%d, binarize=%s", tissue_name, tilesize, binarize))

valid_tissues <- c("deepseq_488B", "deepseq_489", "lowseq_488B", "lowseq_489", "deepseq_combined", "lowseq_combined")
if (!(tissue_name %in% valid_tissues)) stop(sprintf("Invalid tissue_name: %s", tissue_name))
if (!(tilesize %in% c(500, 5000))) stop("tilesize must be 500 or 5000")
if (!is.logical(binarize)) stop("binarize must be TRUE or FALSE")

set.seed(42)
addArchRGenome("hg38")
addArchRThreads(threads = as.integer(Sys.getenv("NSLOTS", "8")))

project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
data_dir <- file.path(project_root, "Data", "01_inputs")
arrow_output_subdir <- if (binarize) "arrow_binarize_v2" else "arrow_not_binarize_v2"
arrow_output_dir <- file.path(data_dir, "arrow", arrow_output_subdir)
dir.create(arrow_output_dir, recursive = TRUE, showWarnings = FALSE)

tissue_metadata <- list(
  deepseq_488B = list(
    fragments = file.path(data_dir, "fragments", "deepseq_488B", "deepseq_488B.fragments.sort.filtered.bed.gz"),
    sample_name = "Deepseq_488B",
    barcodes = file.path(data_dir, "barcodes", "tissue_barcodes", "deepseq_488B", "deepseq_488B.no_edge_effect.barcodes.tsv")
  ),
  deepseq_489 = list(
    fragments = file.path(data_dir, "fragments", "deepseq_489", "deepseq_489.fragments.sort.filtered.bed.gz"),
    sample_name = "Deepseq_489",
    barcodes = file.path(data_dir, "barcodes", "tissue_barcodes", "deepseq_489", "deepseq_489.no_edge_effect.barcodes.tsv")
  ),
  lowseq_488B = list(
    fragments = file.path(data_dir, "fragments", "lowseq_488B", "lowseq_488B.fragments.sort.filtered.bed.gz"),
    sample_name = "Lowseq_488B",
    barcodes = file.path(data_dir, "barcodes", "tissue_barcodes", "lowseq_488B", "lowseq_488B.no_edge_effect.barcodes.tsv")
  ),
  lowseq_489 = list(
    fragments = file.path(data_dir, "fragments", "lowseq_489", "lowseq_489.fragments.sort.filtered.bed.gz"),
    sample_name = "Lowseq_489",
    barcodes = file.path(data_dir, "barcodes", "tissue_barcodes", "lowseq_489", "lowseq_489.no_edge_effect.barcodes.tsv")
  ),
  deepseq_combined = list(
    fragments = file.path(data_dir, "fragments", "deepseq_combined", "deepseq_combined.fragments.sort.filtered.bed.gz"),
    sample_name = "Deepseq_combined",
    barcodes = NULL
  ),
  lowseq_combined = list(
    fragments = file.path(data_dir, "fragments", "lowseq_combined", "lowseq_combined.fragments.sort.filtered.bed.gz"),
    sample_name = "Lowseq_combined",
    barcodes = NULL
  )
)

metadata <- tissue_metadata[[tissue_name]]

valid_barcodes <- NULL
if (!is.null(metadata$barcodes) && file.exists(metadata$barcodes)) {
  bcs <- readLines(metadata$barcodes)
  bcs <- bcs[nzchar(bcs)]
  # Fragment file barcodes carry a "-1" suffix; the no_edge_effect list is bare
  bcs <- ifelse(grepl("-1$", bcs), bcs, paste0(bcs, "-1"))
  valid_barcodes <- unique(bcs)
  log_msg("step", sprintf("Restricting arrow creation to %d no_edge_effect barcodes (huge speedup vs scanning the whole barcode universe)", length(valid_barcodes)))
} else {
  log_msg("warn", sprintf("No barcode whitelist found for %s - arrow will be created over the full barcode universe (slow)", tissue_name))
}
output_arrow <- file.path(arrow_output_dir, sprintf("%s_%dbp.arrow", tissue_name, tilesize))

log_msg("step", sprintf("Fragment file: %s", metadata$fragments))
log_msg("step", sprintf("Output arrow: %s", output_arrow))

if (!file.exists(metadata$fragments)) stop(sprintf("Fragment file not found: %s", metadata$fragments))

if (file.exists(output_arrow) && file.size(output_arrow) > 1e8) {
  log_msg("warn", sprintf("Output arrow already exists, skipping creation: %s", output_arrow))
  quit(status = 0)
}

tryCatch({
  temp_dir <- tempdir()
  temp_arrow_dir <- file.path(temp_dir, sprintf("arrow_creation_v2_%s_%d_%s", tissue_name, tilesize, binarize))
  dir.create(temp_arrow_dir, recursive = TRUE, showWarnings = FALSE)

  current_dir <- getwd()
  setwd(temp_arrow_dir)

  tile_mat_params <- list(tileSize = tilesize, binarize = binarize)

  log_msg("step", sprintf(
    "Creating arrow with permissive QC (minTSS=2, minFrags=100, maxFrags=Inf), TileMatrix tilesize=%d, binarize=%s",
    tilesize, binarize
  ))

  validBarcodes_arg <- NULL
  if (!is.null(valid_barcodes)) {
    validBarcodes_arg <- setNames(list(valid_barcodes), metadata$sample_name)
  }

  arrow_files <- createArrowFiles(
    inputFiles = metadata$fragments,
    sampleNames = metadata$sample_name,
    outputNames = sprintf("%s_%dbp", tissue_name, tilesize),
    validBarcodes = validBarcodes_arg,
    minTSS = 2,
    minFrags = 100,
    maxFrags = Inf,
    addTileMat = TRUE,
    TileMatParams = tile_mat_params,
    addGeneScoreMat = TRUE,
    force = TRUE
  )

  setwd(current_dir)

  arrow_pattern <- sprintf("%s_%dbp\\.arrow$", tissue_name, tilesize)
  created_arrows <- list.files(temp_arrow_dir, pattern = arrow_pattern, full.names = TRUE)

  if (length(created_arrows) > 0) {
    file.copy(created_arrows[1], output_arrow, overwrite = TRUE)
  } else if (length(arrow_files) > 0 && file.exists(arrow_files[1])) {
    file.copy(arrow_files[1], output_arrow, overwrite = TRUE)
  } else {
    stop("createArrowFiles() did not produce output")
  }

  if (!file.exists(output_arrow)) stop("Arrow file not created at output location")

  log_msg("done", sprintf("Successfully created arrow file: %s (%.2f GB)",
                          output_arrow, file.size(output_arrow) / 1e9))
}, error = function(e) {
  log_msg("error", sprintf("Failed to create arrow variant: %s", e$message))
  stop(e$message)
})

log_msg("done", "Arrow variant (v2) creation completed successfully")
