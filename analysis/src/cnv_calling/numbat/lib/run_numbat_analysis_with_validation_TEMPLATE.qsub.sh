#!/bin/bash
#
# run_numbat_analysis_with_validation.qsub.sh - TEMPLATE FOR ANY TISSUE
# 
# Purpose: NUMBAT CNV analysis with CRITICAL validation checks
#          Ensures ATAC matrix and allele counts use SAME cell barcodes
#
# Usage: Copy and customize for each tissue:
#   cp run_numbat_analysis_with_validation.qsub.sh run_numbat_analysis_{TISSUE}_refhca.qsub.sh
#   Edit: tissue and dataset variables
#   qsub run_numbat_analysis_{TISSUE}_refhca.qsub.sh
#

#$ -N numbat_analysis_TEMPLATE
#$ -l h_rt=08:00:00
#$ -pe omp 8
#$ -P paxlab
#$ -l mem_per_core=8G
#$ -o analysis/qsub_logs/numbat_analysis_TEMPLATE_$JOB_ID.out
#$ -e analysis/qsub_logs/numbat_analysis_TEMPLATE_$JOB_ID.err

set -eo pipefail

# Initialize module system
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

module load R
which Rscript && Rscript --version

project_root="/projectnb/paxlab/presh/projects/spatial_atac"
cd "$project_root"

# ════════════════════════════════════════════════════════════════════════════════
# CONFIGURATION - CUSTOMIZE FOR YOUR TISSUE
# ════════════════════════════════════════════════════════════════════════════════
tissue="lowseq_489"      # CHANGE THIS FOR YOUR TISSUE
dataset="lowseq"         # CHANGE THIS FOR YOUR DATASET
output_suffix="refhca"   # CHANGE THIS IF USING DIFFERENT REFERENCE

echo "════════════════════════════════════════════════════════════════"
echo "NUMBAT ANALYSIS WITH VALIDATION"
echo "════════════════════════════════════════════════════════════════"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Tissue: $tissue"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Dataset: $dataset"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Job ID: $JOB_ID"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Hostname: $(hostname)"
echo ""

# ════════════════════════════════════════════════════════════════════════════════
# STEP 0: CRITICAL VALIDATION - CHECK INPUT CONSISTENCY
# ════════════════════════════════════════════════════════════════════════════════
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ════════════════════════════════════════════"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP 0: INPUT VALIDATION"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ════════════════════════════════════════════"

# Run validation script
Rscript analysis/src/cnv_calling/numbat/validate_numbat_inputs.R "$tissue"
VALIDATION_EXIT=$?

if [[ $VALIDATION_EXIT -ne 0 ]]; then
  echo ""
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ VALIDATION FAILED!"
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] The ATAC matrix and allele counts use DIFFERENT barcode files."
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] This will cause NUMBAT to fail with 'No matching cell names' error."
  echo ""
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] REQUIRED ACTION:"
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] Regenerate ATAC matrix with correct barcode file:"
  echo "[$(date +'%Y-%m-%d %H:%M:%S')]   qsub analysis/src/cnv_calling/numbat/regenerate_atac_with_correct_barcodes.qsub.sh"
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] Then resubmit this analysis job."
  exit 1
fi

echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ VALIDATION PASSED - Inputs are consistent"
echo ""

# ════════════════════════════════════════════════════════════════════════════════
# STEP 1: PREPARE PATHS AND VERIFY FILES
# ════════════════════════════════════════════════════════════════════════════════
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ════════════════════════════════════════════"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP 1: PREPARE INPUTS"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ════════════════════════════════════════════"

# Input paths
atac_bin="${project_root}/Data/04_analysis/cnv/numbat/inputs/${tissue}/atac_bin/${tissue}_atac_bin.rds"
allele_df="${project_root}/Data/04_analysis/cnv/numbat/inputs/${tissue}/alleles/${tissue}_atac_allele_counts.tsv.gz"

# Output directory
out_dir="${project_root}/Data/04_analysis/cnv/numbat/results/${tissue}/"
mkdir -p "$out_dir"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] ATAC input: $atac_bin"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Allele input: $allele_df"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Output dir: $out_dir"

if [[ ! -f "$atac_bin" ]]; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ ERROR: ATAC file not found"
  exit 1
fi

if [[ ! -f "$allele_df" ]]; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ ERROR: Allele file not found"
  exit 1
fi

echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ All inputs verified"
echo ""

# ════════════════════════════════════════════════════════════════════════════════
# STEP 2: RUN NUMBAT ANALYSIS
# ════════════════════════════════════════════════════════════════════════════════
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ════════════════════════════════════════════"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP 2: RUNNING NUMBAT CNV ANALYSIS"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ════════════════════════════════════════════"
echo ""

Rscript - <<'EOFR'
# NUMBAT Analysis with ref_hca (built-in reference)

library(numbat)
library(data.table)
library(Matrix)

message(sprintf("[%s] Loading inputs...", Sys.time()))

# Parse command-line-like variables (set via environment or function)
tissue <- Sys.getenv("tissue", "lowseq_489")
dataset <- Sys.getenv("dataset", "lowseq")
project_root <- Sys.getenv("project_root", "/projectnb/paxlab/presh/projects/spatial_atac")
out_dir <- Sys.getenv("out_dir", file.path(project_root, sprintf("Data/04_analysis/cnv/numbat/results/%s/", tissue)))

# Use function args if available (set by parent script)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) tissue <- args[1]
if (length(args) >= 2) dataset <- args[2]
if (length(args) >= 3) project_root <- args[3]
if (length(args) >= 4) out_dir <- args[4]

message(sprintf("[%s] Tissue: %s", Sys.time(), tissue))
message(sprintf("[%s] Dataset: %s", Sys.time(), dataset))
message(sprintf("[%s] Project root: %s", Sys.time(), project_root))
message(sprintf("[%s] Output dir: %s", Sys.time(), out_dir))

# Input paths
atac_bin <- file.path(project_root, "Data/04_analysis/cnv/numbat/inputs", tissue, "atac_bin", 
                      sprintf("%s_atac_bin.rds", tissue))
allele_file <- file.path(project_root, "Data/04_analysis/cnv/numbat/inputs", tissue, "alleles",
                         sprintf("%s_atac_allele_counts.tsv.gz", tissue))

# Load data
message(sprintf("[%s] Loading ATAC bin matrix from: %s", Sys.time(), atac_bin))
count_mat <- readRDS(atac_bin)
message(sprintf("[%s]   Dimensions: %d bins × %d cells", Sys.time(), nrow(count_mat), ncol(count_mat)))

message(sprintf("[%s] Loading allele counts from: %s", Sys.time(), allele_file))
df_allele <- data.table::fread(allele_file)
message(sprintf("[%s]   Loaded: %d variants × %d columns", Sys.time(), nrow(df_allele), ncol(df_allele)))

# Load reference (ref_hca built-in)
message(sprintf("[%s] Loading built-in reference (ref_hca)...", Sys.time()))
data(ref_hca)
message(sprintf("[%s]   Reference loaded", Sys.time()))

# Run NUMBAT
message(sprintf("[%s] Starting NUMBAT analysis...", Sys.time()))
message(sprintf("[%s]   ncores = 8, t = 1e-5", Sys.time()))

numbat_obj <- run_numbat(
  count_mat = count_mat,
  lambdas_ref = ref_hca,
  df_allele = df_allele,
  genome = "hg38",
  t = 1e-5,
  ncores = 8,
  plot = TRUE,
  out_dir = out_dir,
  verbose = TRUE
)

message(sprintf("[%s] ✓ NUMBAT analysis completed successfully", Sys.time()))
message(sprintf("[%s] Results saved to: %s", Sys.time(), out_dir))

# Summary statistics
message(sprintf("[%s] Results Summary:", Sys.time()))
message(sprintf("[%s]   Output files in: %s", Sys.time(), out_dir))
list_files <- list.files(out_dir)
for (f in list_files) {
  fpath <- file.path(out_dir, f)
  fsize <- file.size(fpath)
  message(sprintf("[%s]   - %s (%.1f MB)", Sys.time(), f, fsize / 1e6))
}
EOFR

ANALYSIS_EXIT=$?

echo ""
if [[ $ANALYSIS_EXIT -eq 0 ]]; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ NUMBAT analysis completed successfully"
else
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ NUMBAT analysis failed with exit code $ANALYSIS_EXIT"
  exit $ANALYSIS_EXIT
fi

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ ANALYSIS COMPLETE"
echo "════════════════════════════════════════════════════════════════"
