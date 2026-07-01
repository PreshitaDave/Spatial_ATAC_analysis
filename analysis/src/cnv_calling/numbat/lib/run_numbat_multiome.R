#!/usr/bin/env Rscript
#
# run_numbat_multiome.R
# Purpose: Run NUMBAT CNV analysis in ATAC-bin mode with command-line arguments
# Usage: Rscript run_numbat_multiome.R \
#          --countmat <path> \
#          --alleledf <path> \
#          --out_dir <path> \
#          --ref <path> \
#          --gtf <path> \
#          --parL <path>
#

# Suppress initial warnings for cleaner output
suppressWarnings(suppressMessages({
  library(optparse)
  library(data.table)
  library(Matrix)
  library(numbat)
}))

ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

message(sprintf("[%s] ════════════════════════════════════════════════════════════════", ts()))
message(sprintf("[%s] NUMBAT multiome Analysis (ATAC-bin mode)", ts()))
message(sprintf("[%s] ════════════════════════════════════════════════════════════════", ts()))
message(sprintf("[%s]", ts()))

# ============================================================================
# PARSE COMMAND-LINE ARGUMENTS
# ============================================================================

option_list <- list(
  make_option(c("--countmat"), type="character", help="ATAC bin matrix RDS file"),
  make_option(c("--alleledf"), type="character", help="Allele counts TSV.gz file"),
  make_option(c("--out_dir"), type="character", help="Output directory"),
  make_option(c("--ref"), type="character", help="Lambda reference RDS file"),
  make_option(c("--gtf"), type="character", help="Genomic bins (var220kb.rds)"),
  make_option(c("--parL"), type="character", help="NUMBAT parameters RDS file")
)

parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser, positional_arguments = 0)

# Extract arguments
countmat_file <- opts$countmat
alleledf_file <- opts$alleledf
out_dir <- opts$out_dir
ref_file <- opts$ref
gtf_file <- opts$gtf
parL_file <- opts$parL

# ============================================================================
# VALIDATE INPUTS
# ============================================================================

message(sprintf("[%s] STEP 1: Validating inputs...", ts()))

inputs <- list(
  "Count matrix" = countmat_file,
  "Allele DF" = alleledf_file,
  "Reference" = ref_file,
  "GTF/Bins" = gtf_file,
  "Parameters" = parL_file
)

for (name in names(inputs)) {
  path <- inputs[[name]]
  if (!file.exists(path)) {
    message(sprintf("[%s] ✗ ERROR: %s not found: %s", ts(), name, path))
    quit(status = 1)
  }
  size_mb <- file.size(path) / (1024^2)
  message(sprintf("[%s] ✓ %s: %s (%.1f MB)", ts(), name, basename(path), size_mb))
}

# Create output directory
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(out_dir)) {
  message(sprintf("[%s] ✗ ERROR: Failed to create output directory: %s", ts(), out_dir))
  quit(status = 1)
}
message(sprintf("[%s] ✓ Output directory: %s", ts(), out_dir))
message(sprintf("[%s]", ts()))

# ============================================================================
# LOAD INPUT DATA
# ============================================================================

message(sprintf("[%s] STEP 2: Loading input data...", ts()))

# Load count matrix (ATAC)
message(sprintf("[%s] Loading ATAC count matrix...", ts()))
tryCatch({
  count_mat <- readRDS(countmat_file)
  n_bins <- nrow(count_mat)
  n_cells <- ncol(count_mat)
  message(sprintf("[%s] ✓ ATAC matrix: %d bins × %d cells", ts(), n_bins, n_cells))
}, error = function(e) {
  message(sprintf("[%s] ✗ ERROR loading count matrix: %s", ts(), e$message))
  quit(status = 1)
})

# Load allele counts
message(sprintf("[%s] Loading allele counts...", ts()))
tryCatch({
  df_allele <- data.table::fread(alleledf_file, showProgress = FALSE)
  n_variants <- nrow(df_allele)
  message(sprintf("[%s] ✓ Allele counts: %d variant observations", ts(), n_variants))
}, error = function(e) {
  message(sprintf("[%s] ✗ ERROR loading allele DF: %s", ts(), e$message))
  quit(status = 1)
})

# Load reference lambda
message(sprintf("[%s] Loading lambda reference...", ts()))
tryCatch({
  lambdas_ref <- readRDS(ref_file)
  message(sprintf("[%s] ✓ Lambda reference loaded", ts()))
}, error = function(e) {
  message(sprintf("[%s] ✗ ERROR loading reference: %s", ts(), e$message))
  quit(status = 1)
})

# Load GTF/bins (not strictly needed for run_numbat but verify it exists)
message(sprintf("[%s] Verifying GTF/bins file...", ts()))
tryCatch({
  gtf_bins <- readRDS(gtf_file)
  message(sprintf("[%s] ✓ GTF bins loaded", ts()))
}, error = function(e) {
  message(sprintf("[%s] ✗ ERROR loading GTF: %s", ts(), e$message))
  quit(status = 1)
})

# Load parameters (optional - verify it exists)
message(sprintf("[%s] Verifying parameters file...", ts()))
if (file.exists(parL_file)) {
  tryCatch({
    params <- readRDS(parL_file)
    message(sprintf("[%s] ✓ Parameters file loaded", ts()))
  }, error = function(e) {
    message(sprintf("[%s] ✗ WARNING: Could not load parameters: %s", ts(), e$message))
  })
}

message(sprintf("[%s]", ts()))

# ============================================================================
# STEP 3: RUN NUMBAT
# ============================================================================

message(sprintf("[%s] STEP 3: Running NUMBAT analysis...", ts()))
message(sprintf("[%s] Parameters:", ts()))
message(sprintf("[%s]   - genome: hg38", ts()))
message(sprintf("[%s]   - t: 1e-5", ts()))
message(sprintf("[%s]   - ncores: 8", ts()))
message(sprintf("[%s]   - plot: TRUE", ts()))
message(sprintf("[%s]", ts()))

tryCatch({
  message(sprintf("[%s] Calling run_numbat()...", ts()))
  
  numbat_obj <- run_numbat(
    count_mat = count_mat,
    lambdas_ref = lambdas_ref,
    df_allele = df_allele,
    genome = "hg38",
    t = 1e-5,
    ncores = 8,
    plot = TRUE,
    out_dir = out_dir,
    verbose = TRUE
  )
  
  message(sprintf("[%s] ✓ NUMBAT analysis completed successfully!", ts()))
  
}, error = function(e) {
  message(sprintf("[%s] ✗ ERROR during NUMBAT analysis:", ts()))
  message(sprintf("[%s]   Message: %s", ts(), e$message))
  message(sprintf("[%s]   Call: %s", ts(), paste(deparse(sys.call(-1)), collapse=" ")))
  quit(status = 1)
}, finally = {
  message(sprintf("[%s]", ts()))
})

# ============================================================================
# STEP 4: VERIFY OUTPUTS
# ============================================================================

message(sprintf("[%s] STEP 4: Verifying outputs...", ts()))

output_files <- list.files(out_dir, full.names = FALSE)
if (length(output_files) == 0) {
  message(sprintf("[%s] ✗ WARNING: No output files found in %s", ts(), out_dir))
} else {
  message(sprintf("[%s] ✓ Output files generated:", ts()))
  for (f in sort(output_files)) {
    size <- file.size(file.path(out_dir, f))
    size_str <- if (size > 1024^2) sprintf("%.1f MB", size/(1024^2)) else sprintf("%.1f KB", size/1024)
    message(sprintf("[%s]   - %s (%s)", ts(), f, size_str))
  }
}

message(sprintf("[%s]", ts()))

# ============================================================================
# COMPLETION
# ============================================================================

message(sprintf("[%s] ════════════════════════════════════════════════════════════════", ts()))
message(sprintf("[%s] ✓ NUMBAT Analysis Complete!", ts()))
message(sprintf("[%s] Results saved to: %s", ts(), out_dir))
message(sprintf("[%s] ════════════════════════════════════════════════════════════════", ts()))

quit(status = 0)
