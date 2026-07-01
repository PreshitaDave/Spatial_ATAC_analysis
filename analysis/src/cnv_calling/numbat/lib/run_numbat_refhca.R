#!/usr/bin/env Rscript
#
# run_numbat_refhca.R
# Purpose: Run NUMBAT CNV analysis with refhca reference
# Usage: Rscript run_numbat_refhca.R <tissue> <atac_file> <allele_file> <output_dir> <ncores>
#

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript run_numbat_refhca.R <tissue> <atac_file> <allele_file> <output_dir> <ncores>")
}

tissue <- args[1]
atac_file <- args[2]
allele_file <- args[3]
out_dir <- args[4]
ncores <- as.numeric(args[5])

# Timestamp helper
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

# Startup messages
message(sprintf("[%s] ============================================================", ts()))
message(sprintf("[%s] NUMBAT refhca Analysis", ts()))
message(sprintf("[%s] ============================================================", ts()))
message(sprintf("[%s] R version: %s", ts(), R.version$version.string))
message(sprintf("[%s] Tissue: %s", ts(), tissue))
message(sprintf("[%s] NCORES: %d", ts(), ncores))
message(sprintf("[%s] Hostname: %s", ts(), system("hostname", intern=TRUE)))
message(sprintf("[%s]", ts()))

# ============================================================================
# STEP 1: Load Libraries
# ============================================================================
message(sprintf("[%s] STEP 1: Loading libraries...", ts()))
tryCatch({
  suppressPackageStartupMessages({
    library(numbat)
    library(data.table)
    library(Matrix)
  })
  message(sprintf("[%s] ✓ Libraries loaded successfully", ts()))
}, error = function(e) {
  message(sprintf("[%s] ✗ ERROR loading libraries: %s", ts(), e$message))
  quit(status = 1)
})

# ============================================================================
# STEP 2: Verify Input Files
# ============================================================================
message(sprintf("[%s] STEP 2: Verifying input files...", ts()))

if (!file.exists(atac_file)) {
  message(sprintf("[%s] ✗ ERROR: ATAC file not found: %s", ts(), atac_file))
  quit(status = 1)
}
atac_size <- file.size(atac_file) / (1024^2)
message(sprintf("[%s] ✓ ATAC file: %s (%.1f MB)", ts(), basename(atac_file), atac_size))

if (!file.exists(allele_file)) {
  message(sprintf("[%s] ✗ ERROR: Allele file not found: %s", ts(), allele_file))
  quit(status = 1)
}
allele_size <- file.size(allele_file) / (1024^2)
message(sprintf("[%s] ✓ Allele file: %s (%.1f MB)", ts(), basename(allele_file), allele_size))

# Create output directory
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(out_dir)) {
  message(sprintf("[%s] ✗ ERROR: Failed to create output directory: %s", ts(), out_dir))
  quit(status = 1)
}
message(sprintf("[%s] ✓ Output directory created: %s", ts(), out_dir))
message(sprintf("[%s]", ts()))

# ============================================================================
# STEP 3: Load Input Data
# ============================================================================
message(sprintf("[%s] STEP 3: Loading input data...", ts()))

message(sprintf("[%s] Loading ATAC bin matrix...", ts()))
tryCatch({
  count_mat <- readRDS(atac_file)
  n_bins <- nrow(count_mat)
  n_cells <- ncol(count_mat)
  message(sprintf("[%s] ✓ ATAC loaded: %d bins × %d cells", ts(), n_bins, n_cells))
}, error = function(e) {
  message(sprintf("[%s] ✗ ERROR loading ATAC: %s", ts(), e$message))
  quit(status = 1)
})

message(sprintf("[%s] Loading allele counts...", ts()))
tryCatch({
  df_allele <- data.table::fread(allele_file, showProgress = FALSE)
  n_rows <- nrow(df_allele)
  n_cols <- ncol(df_allele)
  message(sprintf("[%s] ✓ Allele counts loaded: %d rows × %d columns", ts(), n_rows, n_cols))
}, error = function(e) {
  message(sprintf("[%s] ✗ ERROR loading alleles: %s", ts(), e$message))
  quit(status = 1)
})

message(sprintf("[%s] Loading reference...", ts()))
tryCatch({
  data(ref_hca)
  message(sprintf("[%s] ✓ ref_hca reference loaded", ts()))
}, error = function(e) {
  message(sprintf("[%s] ✗ ERROR loading reference: %s", ts(), e$message))
  quit(status = 1)
})

message(sprintf("[%s]", ts()))

# ============================================================================
# STEP 4: Run NUMBAT Analysis
# ============================================================================
message(sprintf("[%s] STEP 4: Running NUMBAT analysis...", ts()))
message(sprintf("[%s] Parameters: t=1e-5, genome=hg38, ncores=%d", ts(), ncores))

tryCatch({
  numbat_obj <- run_numbat(
    count_mat = count_mat,
    lambdas_ref = ref_hca,
    df_allele = df_allele,
    genome = "hg38",
    t = 1e-5,
    ncores = ncores,
    plot = TRUE,
    out_dir = out_dir,
    verbose = TRUE
  )
  message(sprintf("[%s] ✓ NUMBAT analysis completed successfully", ts()))
}, error = function(e) {
  message(sprintf("[%s] ✗ ERROR during NUMBAT analysis: %s", ts(), e$message))
  message(sprintf("[%s] Call: %s", ts(), paste(deparse(sys.call(-1)), collapse=" ")))
  quit(status = 1)
})

# ============================================================================
# STEP 5: Verify Outputs
# ============================================================================
message(sprintf("[%s] STEP 5: Verifying outputs...", ts()))

output_files <- list.files(out_dir, full.names = FALSE)
if (length(output_files) == 0) {
  message(sprintf("[%s] ✗ WARNING: No output files found in %s", ts(), out_dir))
} else {
  message(sprintf("[%s] ✓ Output files created:", ts()))
  for (f in output_files) {
    size <- file.size(file.path(out_dir, f))
    size_str <- if (size > 1024^2) sprintf("%.1f MB", size/(1024^2)) else sprintf("%.1f KB", size/1024)
    message(sprintf("[%s]   - %s (%s)", ts(), f, size_str))
  }
}

message(sprintf("[%s]", ts()))

# ============================================================================
# COMPLETION
# ============================================================================
message(sprintf("[%s] ============================================================", ts()))
message(sprintf("[%s] ✓ NUMBAT Analysis Complete!", ts()))
message(sprintf("[%s] Results saved to: %s", ts(), out_dir))
message(sprintf("[%s] ============================================================", ts()))

quit(status = 0)
