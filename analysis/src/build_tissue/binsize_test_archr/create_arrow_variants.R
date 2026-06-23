#!/usr/bin/env Rscript
# ============================================================================
# create_arrow_variants.R
#
# Create arrow files with different tilesize and binarization options
# directly from fragment files for comparison of sparsity and ArchR performance
#
# Usage:
#   Rscript create_arrow_variants.R <tissue_name> <tilesize> <binarize>
#
# Arguments:
#   tissue_name: One of {deepseq_488B, deepseq_489, lowseq_488B, lowseq_489, deepseq_combined, lowseq_combined}
#   tilesize: 500 or 5000 (in bp)
#   binarize: TRUE or FALSE (logical)
#
# Example:
#   Rscript create_arrow_variants.R deepseq_488B 500 FALSE
#   Rscript create_arrow_variants.R lowseq_489 5000 TRUE
#
# Output:
#   Arrow files in:
#   - /Data/01_inputs/arrow/arrow_binarize/{tissue}_{tilesize}bp.arrow (if binarize=TRUE)
#   - /Data/01_inputs/arrow/arrow_not_binarize/{tissue}_{tilesize}bp.arrow (if binarize=FALSE)
#
# Creates arrow files from fragment files using createArrowFiles() with
# specified TileMatrix parameters (tilesize and binarization).
# This ensures different arrow variants are genuinely different at storage level.
#
# ============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(data.table)
  library(SummarizedExperiment)
})

# Logging function
log_msg <- function(tag, msg) {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), tag, msg))
}

# Helper function to compute and save metrics from an arrow file
compute_and_save_metrics <- function(arrow_file, tissue_name, tilesize, binarize, output_metrics) {
  tryCatch({
    log_msg("step", "Computing tile matrix statistics...")

    temp_dir <- tempdir()
    temp_proj_dir <- file.path(temp_dir, sprintf("metrics_proj_%s_%d_%s", tissue_name, tilesize, binarize))

    proj <- ArchRProject(
      ArrowFiles = arrow_file,
      outputDirectory = temp_proj_dir,
      copyArrows = FALSE
    )

    log_msg("step", sprintf("Project has %d cells", ncol(proj)))

    tile_matrix <- getMatrixFromArrow(
      ArrowFile = getArrowFiles(proj)[1],
      useMatrix = "TileMatrix",
      verbose = FALSE,
      binarize = binarize
    )

    if (is.null(tile_matrix)) {
      log_msg("error", "Tile matrix is NULL")
      return(FALSE)
    }

    # Extract the sparse matrix from SummarizedExperiment
    tile_mat <- assay(tile_matrix)

    n_tiles <- as.numeric(nrow(tile_mat))
    n_cells <- as.numeric(ncol(tile_mat))
    n_nnz <- as.numeric(Matrix::nnzero(tile_mat))
    total_elements <- n_tiles * n_cells
    density <- n_nnz / total_elements
    sparsity <- 1.0 - density

    # Compute per-cell and per-tile coverage
    cell_coverage <- Matrix::colSums(tile_mat)
    tile_coverage <- Matrix::rowSums(tile_mat)

    mean_cell_coverage <- mean(cell_coverage)
    median_cell_coverage <- median(cell_coverage)
    mean_tile_coverage <- mean(tile_coverage)
    median_tile_coverage <- median(tile_coverage)

    log_msg("step", sprintf("Tile matrix: %d tiles × %d cells", n_tiles, n_cells))
    log_msg("step", sprintf("Non-zero elements: %d", n_nnz))
    log_msg("step", sprintf("Density: %.6f, Sparsity: %.4f (%.2f%%)", density, sparsity, 100*sparsity))
    log_msg("step", sprintf("Cell coverage: mean=%.2f, median=%.1f", mean_cell_coverage, median_cell_coverage))
    log_msg("step", sprintf("Tile coverage: mean=%.4f, median=%.1f", mean_tile_coverage, median_tile_coverage))

    # Write metrics to file
    log_msg("step", sprintf("Writing metrics to: %s", output_metrics))
    metrics_text <- sprintf(
      "Tissue: %s\nTilesize: %d bp\nBinarize: %s\nCells: %d\nTiles: %d\nNon-zero elements: %d\nDensity: %.6f\nSparsity: %.4f (%.2f%%)\nMean tiles per cell: %.2f\nMedian tiles per cell: %.1f\nMean cells per tile: %.4f\nMedian cells per tile: %.1f\n",
      tissue_name, tilesize, binarize, n_cells, n_tiles, n_nnz, density, sparsity, 100*sparsity,
      mean_cell_coverage, median_cell_coverage, mean_tile_coverage, median_tile_coverage
    )

    writeLines(metrics_text, con = output_metrics)
    log_msg("done", sprintf("Metrics saved to: %s", output_metrics))

    return(TRUE)
  }, error = function(e) {
    log_msg("error", sprintf("Failed to compute metrics: %s", e$message))
    return(FALSE)
  })
}

# Setup paths and configuration
project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
threads <- as.integer(Sys.getenv("NSLOTS", "8"))
min_tss <- 3
min_frags <- 1000

log_msg("start", "===== Arrow File Creation with Variant Parameters =====")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript create_arrow_variants.R <tissue_name> <tilesize> <binarize>")
}

tissue_name <- args[1]
tilesize <- as.integer(args[2])
binarize_arg <- args[3]
binarize <- if (tolower(binarize_arg) %in% c("true", "t", "1")) TRUE else if (tolower(binarize_arg) %in% c("false", "f", "0")) FALSE else as.logical(binarize_arg)

log_msg("info", sprintf("Parameters: tissue=%s, tilesize=%d, binarize=%s", tissue_name, tilesize, binarize))

# Validate parameters
valid_tissues <- c("deepseq_488B", "deepseq_489", "lowseq_488B", "lowseq_489", "deepseq_combined", "lowseq_combined")
if (!(tissue_name %in% valid_tissues)) {
  stop(sprintf("Invalid tissue_name: %s. Must be one of: %s", tissue_name, paste(valid_tissues, collapse = ", ")))
}
if (!(tilesize %in% c(500, 5000))) {
  stop("tilesize must be 500 or 5000")
}
if (!is.logical(binarize)) {
  stop("binarize must be TRUE or FALSE")
}

# Configure ArchR
set.seed(42)
addArchRGenome("hg38")
addArchRThreads(threads = threads)

# Setup paths
data_dir <- file.path(project_root, "Data", "01_inputs")
arrow_input_dir <- file.path(data_dir, "arrow")
arrow_output_subdir <- if (binarize) "arrow_binarize" else "arrow_not_binarize"
arrow_output_dir <- file.path(arrow_input_dir, arrow_output_subdir)
analysis_output_dir <- file.path(project_root, "analysis", "binsize_comparison")

# Create output directories
dir.create(arrow_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(analysis_output_dir, recursive = TRUE, showWarnings = FALSE)

# Define tissue metadata - points to fragment files (will create arrows from these)
tissue_metadata <- list(
  deepseq_488B = list(
    fragments = file.path(data_dir, "fragments", "deepseq_488B", "deepseq_488B.fragments.sort.filtered.bed.gz"),
    barcodes = file.path(data_dir, "barcodes", "tissue_barcodes", "deepseq_488B", "deepseq_488B.no_edge_effect.barcodes.tsv"),
    sample_name = "Deepseq_488B"
  ),
  deepseq_489 = list(
    fragments = file.path(data_dir, "fragments", "deepseq_489", "deepseq_489.fragments.sort.filtered.bed.gz"),
    barcodes = file.path(data_dir, "barcodes", "tissue_barcodes", "deepseq_489", "deepseq_489.no_edge_effect.barcodes.tsv"),
    sample_name = "Deepseq_489"
  ),
  lowseq_488B = list(
    fragments = file.path(data_dir, "fragments", "lowseq_488B", "lowseq_488B.fragments.sort.filtered.bed.gz"),
    barcodes = file.path(data_dir, "barcodes", "tissue_barcodes", "lowseq_488B", "lowseq_488B.no_edge_effect.barcodes.tsv"),
    sample_name = "Lowseq_488B"
  ),
  lowseq_489 = list(
    fragments = file.path(data_dir, "fragments", "lowseq_489", "lowseq_489.fragments.sort.filtered.bed.gz"),
    barcodes = file.path(data_dir, "barcodes", "tissue_barcodes", "lowseq_489", "lowseq_489.no_edge_effect.barcodes.tsv"),
    sample_name = "Lowseq_489"
  ),
  deepseq_combined = list(
    fragments = file.path(data_dir, "fragments", "deepseq_combined", "deepseq_combined.fragments.sort.filtered.bed.gz"),
    barcodes = NULL,
    sample_name = "Deepseq_combined"
  ),
  lowseq_combined = list(
    fragments = file.path(data_dir, "fragments", "lowseq_combined", "lowseq_combined.fragments.sort.filtered.bed.gz"),
    barcodes = NULL,
    sample_name = "Lowseq_combined"
  )
)

metadata <- tissue_metadata[[tissue_name]]
if (is.null(metadata)) {
  stop(sprintf("Tissue %s not found in metadata", tissue_name))
}

# Construct output filenames
output_arrow <- file.path(arrow_output_dir, sprintf("%s_%dbp.arrow", tissue_name, tilesize))
output_metrics <- file.path(analysis_output_dir, sprintf("%s_%dbp_binarize%s_metrics.txt", tissue_name, tilesize, binarize))

log_msg("step", sprintf("Fragment file: %s", metadata$fragments))
log_msg("step", sprintf("Output arrow: %s", output_arrow))
log_msg("step", sprintf("Metrics output: %s", output_metrics))
log_msg("step", sprintf("Tilesize: %d bp, Binarize: %s", tilesize, binarize))

# Check if fragment file exists
if (!file.exists(metadata$fragments)) {
  log_msg("error", sprintf("Fragment file not found: %s", metadata$fragments))
  stop("Fragment file missing")
}

# Check if output already exists
if (file.exists(output_arrow) && file.size(output_arrow) > 1e8) {  # > 100MB means it's real
  log_msg("warn", sprintf("Output arrow already exists, skipping creation: %s", output_arrow))

  # Always recompute metrics (they may have been skipped or failed in previous run)
  log_msg("step", "Computing metrics for existing arrow file...")
  success <- compute_and_save_metrics(output_arrow, tissue_name, tilesize, binarize, output_metrics)

  if (success) {
    log_msg("done", "Metrics computation succeeded")
    quit(status = 0)
  } else {
    log_msg("error", "Metrics computation failed")
    stop("Could not compute metrics for existing arrow file")
  }
}

tryCatch({
  log_msg("step", "Creating arrow file from fragments using createArrowFiles()...")
  
  # Setup temporary directory for arrow creation
  temp_dir <- tempdir()
  temp_arrow_dir <- file.path(temp_dir, sprintf("arrow_creation_%s_%d_%s", tissue_name, tilesize, binarize))
  dir.create(temp_arrow_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Change to temp directory for createArrowFiles (it creates files in working directory)
  current_dir <- getwd()
  setwd(temp_arrow_dir)
  
  # Define TileMatrix parameters with specified tilesize and binarization
  # binarize=TRUE: stores binary (0/1) in arrow
  # binarize=FALSE: stores count values in arrow
  tile_mat_params <- list(
    tileSize = tilesize,
    binarize = binarize
  )
  
  log_msg("step", sprintf("Creating arrow with TileMatrix params: tilesize=%d, binarize=%s", tilesize, binarize))
  
  # Create arrow file from fragment file
  # This is the key difference: arrow is created with specified binarization at storage level
  arrow_files <- createArrowFiles(
    inputFiles = metadata$fragments,
    sampleNames = metadata$sample_name,
    outputNames = sprintf("%s_%dbp", tissue_name, tilesize),
    minTSS = min_tss,
    minFrags = min_frags,
    addTileMat = TRUE,
    TileMatParams = tile_mat_params,  # Binarization is stored in arrow!
    addGeneScoreMat = TRUE,
    force = TRUE
  )
  
  # Return to original directory
  setwd(current_dir)
  
  log_msg("step", sprintf("Arrow creation completed. Files: %s", paste(arrow_files, collapse=", ")))
  
  # Find and copy the created arrow file to output directory
  # createArrowFiles creates files in the working directory, so look in temp_arrow_dir
  arrow_pattern <- sprintf("%s_%dbp\\.arrow$", tissue_name, tilesize)
  created_arrows <- list.files(temp_arrow_dir, pattern = arrow_pattern, full.names = TRUE)
  
  if (length(created_arrows) > 0) {
    arrow_file_source <- created_arrows[1]
    file.copy(arrow_file_source, output_arrow, overwrite = TRUE)
    log_msg("step", sprintf("Copied arrow file to output location: %s", output_arrow))
  } else if (length(arrow_files) > 0 && file.exists(arrow_files[1])) {
    file.copy(arrow_files[1], output_arrow, overwrite = TRUE)
    log_msg("step", sprintf("Copied arrow file to output location: %s", output_arrow))
  } else {
    log_msg("error", "Arrow file creation failed - no output files found")
    stop("createArrowFiles() did not produce output")
  }

  # Verify output file exists
  if (!file.exists(output_arrow)) {
    log_msg("error", sprintf("Output arrow file not found: %s", output_arrow))
    stop("Arrow file not created at output location")
  }

  output_size <- file.size(output_arrow)
  log_msg("done", sprintf("Successfully created arrow file: %s (%.2f GB)",
                          output_arrow, output_size / 1e9))

  # Now compute metrics
  log_msg("step", "Computing metrics for newly created arrow file...")
  success <- compute_and_save_metrics(output_arrow, tissue_name, tilesize, binarize, output_metrics)

  if (success) {
    log_msg("done", sprintf("Successfully created and analyzed arrow variant: %s (tilesize=%d, binarize=%s)",
                            tissue_name, tilesize, binarize))
  } else {
    log_msg("error", "Metrics computation failed after arrow creation")
    stop("Could not compute metrics for created arrow file")
  }
  
}, error = function(e) {
  log_msg("error", sprintf("Failed to create arrow variant: %s", e$message))
  stop(e$message)
})

log_msg("done", "Arrow variant creation completed successfully")
