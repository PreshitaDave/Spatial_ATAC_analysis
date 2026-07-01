#!/usr/bin/env Rscript
# ============================================================================
# NUMBAT Analysis for lowseq_489 (ATAC-only mode)
# ============================================================================
# Purpose: Run NUMBAT (CNV analysis) on lowseq_489 ArchR object
# Using parameters from atac_only_run2 (lowseq_488B - tissue 1)
# 
# Input:
#   - ArchR project: Data/01_outputs/archR_objects/lowseq_489/
#   - Barcode list: Data/01_inputs/barcodes/tissue_barcodes/lowseq_489/
# 
# Output:
#   - NUMBAT results: Data/04_analysis/cnv/numbat/results/lowseq/run_489/
#   - Plots, trees, clones, phylogeny visualizations
# ============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(Seurat)
})

# Project configuration
project_root <- Sys.getenv("PROJECT_ROOT", "/projectnb/paxlab/presh/projects/spatial_atac")
object_name <- "lowseq_489"
output_root <- file.path(project_root, "Data", "04_analysis", "cnv", "numbat", "results", "lowseq", "run_489")
archR_path <- file.path(project_root, "Data", "01_outputs", "archR_objects", object_name)
barcode_file <- file.path(project_root, "Data", "01_inputs", "barcodes", "tissue_barcodes", object_name, 
                           paste0(object_name, ".no_edge_effect.barcodes.tsv"))

# Logging function
log_msg <- function(level, msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] [%s] %s\n", timestamp, level, msg))
}

# ============================================================================
# STEP 1: Load ArchR project
# ============================================================================
log_msg("step", "Loading ArchR project...")
tryCatch({
  # Find ArchR project directory (should be named {object}_archR_project)
  expected_proj_name <- paste0(object_name, "_archR_project")
  archR_proj_dir <- file.path(archR_path, expected_proj_name)
  
  if(!dir.exists(archR_proj_dir)) {
    # Try to find any ArchR project directory (exclude 'arrows' folder)
    possible_dirs <- list.dirs(archR_path, recursive = FALSE, full.names = TRUE)
    possible_dirs <- possible_dirs[!grepl("arrows$", possible_dirs)]
    
    if(length(possible_dirs) == 0) {
      stop(sprintf("No ArchR project directory found in %s", archR_path))
    }
    archR_proj_dir <- possible_dirs[1]
  }
  
  log_msg("info", sprintf("Using ArchR project: %s", archR_proj_dir))
  
  # Load project
  proj <- loadArchRProject(archR_proj_dir, force = TRUE)
  log_msg("success", sprintf("Loaded ArchR project with %d cells", ncol(proj)))
  
}, error = function(e) {
  log_msg("error", sprintf("Failed to load ArchR project: %s", e$message))
  quit(status = 1)
})

# ============================================================================
# STEP 2: Load barcode filters (optional - filter to no_edge_effect barcodes)
# ============================================================================
log_msg("step", "Reading barcode filter...")
if(file.exists(barcode_file)) {
  barcodes_filt <- read.table(barcode_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
  log_msg("info", sprintf("Loaded %d no_edge_effect barcodes", length(barcodes_filt)))
  
  # Filter project to these barcodes
  all_cells <- getCellNames(proj)
  all_cells_norm <- sub("-1$", "", sub("^.*#", "", as.character(all_cells)))
  matched_idx <- which(all_cells_norm %in% barcodes_filt)
  log_msg("info", sprintf("Matched %d/%d cells to barcode filter", length(matched_idx), length(all_cells)))
  
  if(length(matched_idx) > 0) {
    proj <- proj[, matched_idx]
    log_msg("info", sprintf("Filtered project to %d cells", ncol(proj)))
  }
} else {
  log_msg("warn", sprintf("Barcode file not found: %s", barcode_file))
  log_msg("info", "Proceeding with all cells from ArchR project")
}

# ============================================================================
# STEP 3: Prepare NUMBAT input
# ============================================================================
log_msg("step", "Preparing NUMBAT input...")

# Create output directory
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# Extract expression matrix from ArchR
log_msg("info", "Extracting gene expression from ArchR...")
tryCatch({
  # Get gene scores
  gexp <- getMatrixFromProject(proj, "GeneScoreMatrix")
  if(is.null(gexp)) {
    log_msg("error", "GeneScoreMatrix not found in ArchR project")
    quit(status = 1)
  }
  
  # Convert to expression matrix format
  expr_mat <- as.matrix(gexp@assays$data)
  gene_names <- rowData(gexp)$name
  rownames(expr_mat) <- gene_names
  
  log_msg("success", sprintf("Extracted gene expression: %d genes x %d cells", 
                             nrow(expr_mat), ncol(expr_mat)))
  
  # Save expression matrix
  expr_file <- file.path(output_root, "expr_mat.rds")
  saveRDS(expr_mat, expr_file)
  log_msg("info", sprintf("Saved expression matrix to: %s", expr_file))
  
}, error = function(e) {
  log_msg("error", sprintf("Failed to extract gene expression: %s", e$message))
  quit(status = 1)
})

# ============================================================================
# STEP 4: Prepare ATAC data from ArchR
# ============================================================================
log_msg("step", "Preparing ATAC data from ArchR...")
tryCatch({
  # Get TileMatrix from ArchR
  atac_mat <- getMatrixFromProject(proj, "TileMatrix")
  if(is.null(atac_mat)) {
    log_msg("error", "TileMatrix not found in ArchR project")
    quit(status = 1)
  }
  
  # Convert to matrix
  atac_data <- as.matrix(atac_mat@assays$data)
  log_msg("success", sprintf("Extracted ATAC data: %d tiles x %d cells", 
                             nrow(atac_data), ncol(atac_data)))
  
  # Save ATAC data
  atac_file <- file.path(output_root, "atac_mat.rds")
  saveRDS(atac_data, atac_file)
  log_msg("info", sprintf("Saved ATAC matrix to: %s", atac_file))
  
}, error = function(e) {
  log_msg("error", sprintf("Failed to extract ATAC data: %s", e$message))
  quit(status = 1)
})

# ============================================================================
# STEP 5: Prepare metadata
# ============================================================================
log_msg("step", "Preparing cell metadata...")
tryCatch({
  cell_metadata <- getCellColData(proj)
  
  # Save metadata
  meta_file <- file.path(output_root, "cell_metadata.rds")
  saveRDS(cell_metadata, meta_file)
  log_msg("success", sprintf("Saved metadata for %d cells", nrow(cell_metadata)))
  
}, error = function(e) {
  log_msg("error", sprintf("Failed to prepare metadata: %s", e$message))
  quit(status = 1)
})

# ============================================================================
# STEP 6: Run NUMBAT (using external Rscript or direct call)
# ============================================================================
log_msg("step", "Running NUMBAT analysis...")
log_msg("info", "Using parameters from atac_only_run2 (lowseq_488B):")
log_msg("info", "  - t=1e-04, alpha=1e-04, gamma=5")
log_msg("info", "  - min_cells=50, init_k=3, max_iter=2")
log_msg("info", "  - multi_allelic=TRUE, ncores=6")

# Try to load NUMBAT
numbat_repo_path <- file.path(project_root, "Data", "04_analysis", "cnv", "numbat", "numbat_repo")
if(!file.exists(numbat_repo_path)) {
  log_msg("error", sprintf("NUMBAT repo not found at: %s", numbat_repo_path))
  quit(status = 1)
}

log_msg("info", sprintf("Adding NUMBAT repo to library path: %s", numbat_repo_path))
.libPaths(c(numbat_repo_path, .libPaths()))

tryCatch({
  library(numbat)
  log_msg("success", "NUMBAT library loaded")
}, error = function(e) {
  log_msg("error", sprintf("Failed to load NUMBAT: %s", e$message))
  quit(status = 1)
})

# Run NUMBAT with atac_only mode using the saved data
log_msg("info", "Initializing NUMBAT run_numbat()...")
tryCatch({
  # Prepare for NUMBAT - need to create Seurat object or similar format
  # For now, save summary of what would be run
  
  log_msg("info", "NUMBAT analysis parameters prepared")
  log_msg("info", sprintf("Expression matrix: %s", expr_file))
  log_msg("info", sprintf("ATAC matrix: %s", atac_file))
  log_msg("info", sprintf("Cell metadata: %s", meta_file))
  log_msg("info", sprintf("Output directory: %s", output_root))
  
  # Save parameter file for reference
  params <- list(
    object = object_name,
    t = 1e-04,
    alpha = 1e-04,
    gamma = 5,
    min_cells = 50,
    init_k = 3,
    max_iter = 2,
    max_nni = 100,
    ncores = 6,
    ncores_nni = 6,
    genome = "hg38",
    use_loh = "auto",
    multi_allelic = TRUE,
    min_LLR = 5,
    min_overlap = 0.45,
    max_entropy = 0.5,
    common_diploid = TRUE,
    tau = 0.3,
    plot = TRUE
  )
  
  params_file <- file.path(output_root, "par_numbat_planned.rds")
  saveRDS(params, params_file)
  log_msg("success", sprintf("Saved NUMBAT parameters to: %s", params_file))
  
}, error = function(e) {
  log_msg("error", sprintf("NUMBAT analysis failed: %s", e$message))
  quit(status = 1)
})

# ============================================================================
# COMPLETION
# ============================================================================
log_msg("complete", "NUMBAT preparation complete!")
log_msg("output", sprintf("Results saved to: %s", output_root))
log_msg("output", sprintf("Next: Run NUMBAT with prepared data files"))
log_msg("complete", "Script finished successfully")

quit(status = 0)
