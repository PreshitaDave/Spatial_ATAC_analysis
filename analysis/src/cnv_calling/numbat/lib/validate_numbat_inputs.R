#!/usr/bin/env Rscript
#
# validate_numbat_inputs.R
# Purpose: Validate barcode consistency between NUMBAT input files
# Ensures ATAC matrix and allele counts use the same cell barcodes
# 
# Usage: Rscript validate_numbat_inputs.R <tissue> [--fix]
#        Rscript validate_numbat_inputs.R lowseq_489 --fix
#

library(data.table)
library(Matrix)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript validate_numbat_inputs.R <tissue> [--fix]")
}

tissue <- args[1]
fix_if_needed <- "--fix" %in% args

# Project root
project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"

message(sprintf("[%s] NUMBAT INPUT VALIDATION: %s", Sys.time(), tissue))
message(sprintf("[%s] ════════════════════════════════════════════", Sys.time()))

# Paths
atac_file <- file.path(project_root, sprintf("Data/04_analysis/cnv/numbat/inputs/%s/atac_bin/%s_atac_bin.rds", tissue, tissue))
allele_file <- file.path(project_root, sprintf("Data/04_analysis/cnv/numbat/inputs/%s/alleles/%s_atac_allele_counts.tsv.gz", tissue, tissue))
barcode_for_pileup <- file.path(project_root, sprintf("Data/04_analysis/cnv/numbat/inputs/%s/barcodes/%s_atac_barcodes_for_pileup.tsv", tissue, tissue))

# ════════════════════════════════════════════════════════════════════════════════
# STEP 1: Load ATAC matrix and extract cell names
# ════════════════════════════════════════════════════════════════════════════════
message(sprintf("[%s] STEP 1: Loading ATAC bin matrix...", Sys.time()))
if (!file.exists(atac_file)) {
  stop(sprintf("ATAC file not found: %s", atac_file))
}

count_mat <- readRDS(atac_file)
atac_cells <- colnames(count_mat)
n_atac_cells <- length(atac_cells)
message(sprintf("[%s]   - ATAC cells: %d", Sys.time(), n_atac_cells))
message(sprintf("[%s]   - ATAC bins: %d", Sys.time(), nrow(count_mat)))

# ════════════════════════════════════════════════════════════════════════════════
# STEP 2: Load allele counts and extract unique cells
# ════════════════════════════════════════════════════════════════════════════════
message(sprintf("[%s] STEP 2: Loading allele counts file...", Sys.time()))
if (!file.exists(allele_file)) {
  stop(sprintf("Allele file not found: %s", allele_file))
}

df_allele <- data.table::fread(allele_file)
allele_cells <- unique(df_allele$cell)
n_allele_cells <- length(allele_cells)
n_allele_rows <- nrow(df_allele)
message(sprintf("[%s]   - Allele cells: %d", Sys.time(), n_allele_cells))
message(sprintf("[%s]   - Allele rows (variants): %d", Sys.time(), n_allele_rows))

# ════════════════════════════════════════════════════════════════════════════════
# STEP 3: Check barcode overlap
# ════════════════════════════════════════════════════════════════════════════════
message(sprintf("[%s] STEP 3: Checking barcode overlap...", Sys.time()))

overlap <- intersect(atac_cells, allele_cells)
n_overlap <- length(overlap)

atac_only <- setdiff(atac_cells, allele_cells)
n_atac_only <- length(atac_only)

allele_only <- setdiff(allele_cells, atac_cells)
n_allele_only <- length(allele_only)

message(sprintf("[%s]   - Overlap (cells in both): %d", Sys.time(), n_overlap))
message(sprintf("[%s]   - ATAC-only cells: %d", Sys.time(), n_atac_only))
message(sprintf("[%s]   - Allele-only cells: %d", Sys.time(), n_allele_only))

# ════════════════════════════════════════════════════════════════════════════════
# STEP 4: Check pileup barcode file
# ════════════════════════════════════════════════════════════════════════════════
message(sprintf("[%s] STEP 4: Checking pileup barcode file...", Sys.time()))
if (!file.exists(barcode_for_pileup)) {
  message(sprintf("[%s]   ⚠ Pileup barcode file not found (optional check)", Sys.time()))
  pileup_barcodes <- NULL
} else {
  pileup_barcodes <- data.table::fread(barcode_for_pileup, header = FALSE)[[1]]
  n_pileup <- length(pileup_barcodes)
  message(sprintf("[%s]   - Pileup barcodes: %d", Sys.time(), n_pileup))
  
  # Check overlap with allele cells
  pileup_allele_overlap <- intersect(pileup_barcodes, allele_cells)
  message(sprintf("[%s]   - Overlap (pileup ∩ allele): %d", Sys.time(), length(pileup_allele_overlap)))
}

# ════════════════════════════════════════════════════════════════════════════════
# STEP 5: Verdict
# ════════════════════════════════════════════════════════════════════════════════
message(sprintf("[%s] STEP 5: Validation Verdict", Sys.time()))
message(sprintf("[%s] ────────────────────────────────────────────", Sys.time()))

is_consistent <- (n_atac_cells == n_allele_cells) && (n_overlap == n_atac_cells) && (n_atac_only == 0) && (n_allele_only == 0)

if (is_consistent) {
  message(sprintf("[%s] ✓ PASS: Barcode files are CONSISTENT", Sys.time()))
  message(sprintf("[%s] All %d cells match between ATAC and allele counts", Sys.time(), n_atac_cells))
  quit(status = 0)
} else {
  message(sprintf("[%s] ✗ FAIL: Barcode files are INCONSISTENT", Sys.time()))
  message(sprintf("[%s] Issue: Different barcode sets used in ATAC vs pileup stages", Sys.time()))
  
  if (n_allele_only > 0) {
    message(sprintf("[%s]   → %d cells in allele file but NOT in ATAC matrix", Sys.time(), n_allele_only))
    message(sprintf("[%s]   → These cells will be FILTERED OUT by NUMBAT (0 coverage)", Sys.time()))
  }
  if (n_atac_only > 0) {
    message(sprintf("[%s]   → %d cells in ATAC matrix but NOT in allele file", Sys.time(), n_atac_only))
    message(sprintf("[%s]   → These cells will have NO VARIANT DATA in NUMBAT", Sys.time()))
  }
  
  if (!fix_if_needed) {
    message(sprintf("[%s] ────────────────────────────────────────────", Sys.time()))
    message(sprintf("[%s] RECOMMENDED FIX:", Sys.time()))
    message(sprintf("[%s] Regenerate ATAC matrix with pileup barcode file:", Sys.time()))
    message(sprintf("[%s]   tissue <- '%s'", Sys.time(), tissue))
    message(sprintf("[%s]   Rscript analysis/src/cnv_calling/numbat/get_binned_atac_fixed.R \\", Sys.time()))
    message(sprintf("[%s]     --CB Data/04_analysis/cnv/numbat/inputs/{tissue}/barcodes/{tissue}_atac_barcodes_for_pileup.tsv \\", Sys.time()))
    message(sprintf("[%s]     --frag Data/01_inputs/fragments/{tissue}/ \\", Sys.time()))
    message(sprintf("[%s]     --outFile Data/04_analysis/cnv/numbat/inputs/{tissue}/atac_bin/{tissue}_atac_bin_FIXED.rds", Sys.time()))
    message(sprintf("[%s]", Sys.time()))
    message(sprintf("[%s] Then update references and rerun NUMBAT analysis", Sys.time()))
    quit(status = 1)
  } else {
    message(sprintf("[%s] ────────────────────────────────────────────", Sys.time()))
    message(sprintf("[%s] --fix flag provided. Would regenerate, but not implemented in validation script.", Sys.time()))
    message(sprintf("[%s] Please run the regeneration command above manually.", Sys.time()))
    quit(status = 1)
  }
}
