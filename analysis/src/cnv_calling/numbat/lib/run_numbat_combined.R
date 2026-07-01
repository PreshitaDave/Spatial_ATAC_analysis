#!/usr/bin/env Rscript
# NUMBAT analysis on combined lowseq tissues (488B + 489)
# Purpose: CNV calling and allele phasing on combined spatial ATAC

# Suppress warnings for cleaner output
options(warn = -1)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_numbat_combined.R <dataset> <tissue>")
}
DATASET <- args[1]  # "lowseq"
TISSUE <- args[2]   # "combined"

cat(sprintf("[START] NUMBAT analysis: %s_%s\n", DATASET, TISSUE))

# Load libraries
cat("[STEP 1] Loading required libraries...\n")
library(ArchR, quietly = TRUE)
library(numbat, quietly = TRUE)
set.seed(1)

# Define paths
PROJECT_ROOT <- "/projectnb/paxlab/presh/projects/spatial_atac"
DATASET_TISSUE <- sprintf("%s_%s", DATASET, TISSUE)

# ArchR project path
ARCHR_PROJECT_PATH <- file.path(
  PROJECT_ROOT, "Data/01_outputs/archR_objects",
  DATASET_TISSUE, sprintf("%s_archR_project_final", DATASET_TISSUE)
)

# Input directory for NUMBAT
NUMBAT_INPUT_DIR <- file.path(
  PROJECT_ROOT, "Data/04_analysis/cnv/numbat/inputs",
  DATASET_TISSUE
)

# Output directory for NUMBAT results
NUMBAT_OUTPUT_DIR <- file.path(
  PROJECT_ROOT, "Data/04_analysis/cnv/numbat",
  DATASET_TISSUE
)

# Verify directories exist
cat("[STEP 2] Verifying directories...\n")
if (!dir.exists(ARCHR_PROJECT_PATH)) {
  stop(sprintf("ERROR: ArchR project not found at %s", ARCHR_PROJECT_PATH))
}
if (!dir.exists(NUMBAT_INPUT_DIR)) {
  stop(sprintf("ERROR: NUMBAT input directory not found at %s", NUMBAT_INPUT_DIR))
}
cat(sprintf("✓ ArchR project: %s\n", ARCHR_PROJECT_PATH))
cat(sprintf("✓ NUMBAT input dir: %s\n", NUMBAT_INPUT_DIR))

# Create output directory
dir.create(NUMBAT_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# STEP 3: Load ArchR project and extract metadata
cat("[STEP 3] Loading ArchR project...\n")
proj <- loadArchRProject(ARCHR_PROJECT_PATH, showLogo = FALSE)
cat(sprintf("✓ Loaded %d cells\n", ncol(proj)))

# Get cell barcodes
cell_barcodes <- colnames(proj)
metadata <- as.data.frame(getCellColData(proj))
cat(sprintf("✓ Metadata loaded: %d cells\n", nrow(metadata)))

# STEP 4: Load NUMBAT input files
cat("[STEP 4] Loading NUMBAT input files...\n")

# Look for barcode file
BARCODE_FILE <- file.path(NUMBAT_INPUT_DIR, sprintf("barcode_%s_%s.txt", DATASET, TISSUE))
if (!file.exists(BARCODE_FILE)) {
  # Try alternative naming
  BARCODE_FILE <- list.files(NUMBAT_INPUT_DIR, pattern = "barcode.*\\.txt$", full.names = TRUE)[1]
  if (is.na(BARCODE_FILE)) {
    stop(sprintf("ERROR: Could not find barcode file in %s", NUMBAT_INPUT_DIR))
  }
}
cat(sprintf("✓ Barcode file: %s\n", basename(BARCODE_FILE)))

# Look for bin matrix file
BIN_MATRIX_FILE <- list.files(NUMBAT_INPUT_DIR, pattern = "bincounts.*\\.h5$", full.names = TRUE)[1]
if (is.na(BIN_MATRIX_FILE)) {
  stop(sprintf("ERROR: Could not find bin matrix file in %s", NUMBAT_INPUT_DIR))
}
cat(sprintf("✓ Bin matrix file: %s\n", basename(BIN_MATRIX_FILE)))

# Look for allele file
ALLELE_FILE <- list.files(NUMBAT_INPUT_DIR, pattern = "allele_*.h5", full.names = TRUE)[1]
if (is.na(ALLELE_FILE)) {
  stop(sprintf("ERROR: Could not find allele file in %s", NUMBAT_INPUT_DIR))
}
cat(sprintf("✓ Allele file: %s\n", basename(ALLELE_FILE)))

# STEP 5: Run NUMBAT
cat("[STEP 5] Running NUMBAT analysis...\n")

# Load reference files
VAR220_FILE <- file.path(PROJECT_ROOT, "Data/02_references/var220kb.rds")
LAMBDAS_FILE <- file.path(PROJECT_ROOT, "Data/02_references/lambdas_ATAC_bincnt.rds")

if (!file.exists(VAR220_FILE)) {
  stop(sprintf("ERROR: Reference file not found: %s", VAR220_FILE))
}
if (!file.exists(LAMBDAS_FILE)) {
  stop(sprintf("ERROR: Reference file not found: %s", LAMBDAS_FILE))
}
cat(sprintf("✓ Reference files loaded\n"))

# Initialize NUMBAT
cat("[STEP 6] Initializing NUMBAT object...\n")
numbat_obj <- initialize_numbat(
  ref_var = VAR220_FILE,
  lambdas_pooled = LAMBDAS_FILE,
  tumor_cells = readLines(BARCODE_FILE),
  bincounts_file = BIN_MATRIX_FILE,
  allele_file = ALLELE_FILE,
  seurat_obj = NULL,
  ncores = 8
)
cat(sprintf("✓ NUMBAT object initialized with %d cells\n", nrow(numbat_obj@allele)))

# Run bulk mode
cat("[STEP 7] Running bulk allele calling...\n")
numbat_obj <- run_bulk_allele(numbat_obj)

# Run clustering
cat("[STEP 8] Running CNV clustering...\n")
numbat_obj <- run_clonealign(numbat_obj, verbose = TRUE)

# Generate plots
cat("[STEP 9] Generating plots...\n")
plot_dir <- file.path(NUMBAT_OUTPUT_DIR, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

pdf(file.path(plot_dir, "heatmap.pdf"), width = 12, height = 8)
plot_phylo_heatmap(numbat_obj)
dev.off()

pdf(file.path(plot_dir, "umap.pdf"), width = 10, height = 8)
plot_dimred(numbat_obj)
dev.off()

# STEP 10: Extract and save results
cat("[STEP 10] Extracting and saving results...\n")

# Save NUMBAT object
FINAL_OBJ_PATH <- file.path(NUMBAT_OUTPUT_DIR, "final_obj.rds")
saveRDS(numbat_obj, FINAL_OBJ_PATH)
cat(sprintf("✓ NUMBAT object saved: %s\n", FINAL_OBJ_PATH))

# Extract allele posteriors
allele_post <- as.data.frame(numbat_obj@allele)
for (col in colnames(allele_post)) {
  if (grepl("post_", col)) {
    post_file <- file.path(NUMBAT_OUTPUT_DIR, sprintf("allele_%s.tsv", col))
    write.table(allele_post[, c("cell_barcode", col)], 
                post_file, 
                sep = "\t", 
                quote = FALSE, 
                row.names = FALSE)
  }
}
cat(sprintf("✓ Allele posteriors saved\n"))

# Extract CNV calls
cnv_calls <- as.data.frame(numbat_obj@cnv)
CNV_FILE <- file.path(NUMBAT_OUTPUT_DIR, "cnv_calls.tsv")
write.table(cnv_calls, CNV_FILE, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("✓ CNV calls saved: %s\n", CNV_FILE))

# STEP 11: Summary statistics
cat("\n[SUMMARY]\n")
cat(sprintf("NUMBAT Analysis: %s_%s\n", DATASET, TISSUE))
cat(sprintf("  - Cells analyzed: %d\n", nrow(numbat_obj@allele)))
cat(sprintf("  - CNV clones detected: %d\n", length(unique(numbat_obj@clone))))
cat(sprintf("  - Output directory: %s\n", NUMBAT_OUTPUT_DIR))
cat("\n✓ NUMBAT analysis complete!\n")
