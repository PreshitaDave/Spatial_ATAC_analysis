#!/bin/bash
################################################################################
# run_numbat_analysis_lowseq_488B_refhca.qsub.sh
#
# PURPOSE: Run NUMBAT CNV analysis for lowseq_488B using built-in ref_hca
#          WITH CRITICAL VALIDATION of input consistency
#
# NEW FEATURE: Validates that ATAC matrix and allele counts use same cell barcodes
#              Prevents silent failures from barcode file mismatches
#
################################################################################

#$ -N numbat_analysis_lowseq_488B_refhca
#$ -l h_rt=08:00:00
#$ -pe omp 8
#$ -P paxlab
#$ -l mem_per_core=8G
#$ -o analysis/qsub_logs/numbat_analysis_lowseq_488B_refhca_$JOB_ID.out
#$ -e analysis/qsub_logs/numbat_analysis_lowseq_488B_refhca_$JOB_ID.err

set -eo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
TISSUE="lowseq_488B"
NCORES=8

cd "$PROJECT_ROOT"

echo "════════════════════════════════════════════════════════════════════════════"
echo "NUMBAT ANALYSIS: $TISSUE with VALIDATION"
echo "════════════════════════════════════════════════════════════════════════════"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Job ID: $JOB_ID"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Hostname: $(hostname)"
echo ""

# Initialize module system
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

# Load R
module load R 2>/dev/null || true
which Rscript && Rscript --version

echo ""
echo "════════════════════════════════════════════════════════════════════════════"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP 0: CRITICAL VALIDATION"
echo "════════════════════════════════════════════════════════════════════════════"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Checking input consistency..."
echo ""

# Run validation
Rscript analysis/src/cnv_calling/numbat/validate_numbat_inputs.R "$TISSUE"
VALIDATION_EXIT=$?

if [[ $VALIDATION_EXIT -ne 0 ]]; then
  echo ""
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗✗✗ VALIDATION FAILED ✗✗✗"
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] "
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] PROBLEM: ATAC matrix and allele counts use DIFFERENT barcode files"
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] This causes NUMBAT to fail with:"
  echo "[$(date +'%Y-%m-%d %H:%M:%S')]   'No matching cell names between count_mat and df_allele'"
  echo "[$(date +'%Y-%m-%d %H:%M:%S')]   'Filtering out X cells with 0 coverage'"
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] "
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] SOLUTION: Regenerate ATAC matrix with correct (pileup) barcode file"
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] "
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] NOTE: For 488B, first ensure prepare_atac_488B.qsub.sh completes successfully"
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] "
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] Then resubmit this analysis job after regeneration completes."
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] "
  exit 1
fi

echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ VALIDATION PASSED"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ATAC matrix and allele counts have consistent cell barcodes"
echo ""

echo "════════════════════════════════════════════════════════════════════════════"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP 1: PREPARE INPUTS"
echo "════════════════════════════════════════════════════════════════════════════"

ATAC_BIN="Data/04_analysis/cnv/numbat/inputs/${TISSUE}/atac_bin/${TISSUE}_atac_bin.rds"
ALLELE_DF="Data/04_analysis/cnv/numbat/inputs/${TISSUE}/alleles/${TISSUE}_atac_allele_counts.tsv.gz"
OUTPUT_DIR="Data/04_analysis/cnv/numbat/results/${TISSUE}/ref_hca"

mkdir -p "$OUTPUT_DIR"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] ATAC matrix: $ATAC_BIN"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Allele counts: $ALLELE_DF"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Output dir: $OUTPUT_DIR"

if [[ ! -f "$ATAC_BIN" ]]; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ ERROR: ATAC file not found: $ATAC_BIN"
  exit 1
fi

if [[ ! -f "$ALLELE_DF" ]]; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ ERROR: Allele file not found: $ALLELE_DF"
  exit 1
fi

echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ All input files verified"
echo ""

echo "════════════════════════════════════════════════════════════════════════════"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP 2: RUNNING NUMBAT ANALYSIS"
echo "════════════════════════════════════════════════════════════════════════════"
echo ""

# Run NUMBAT analysis in R
Rscript - "$TISSUE" "$ATAC_BIN" "$ALLELE_DF" "$OUTPUT_DIR" <<'EOFR'
args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
atac_file <- args[2]
allele_file <- args[3]
out_dir <- args[4]

library(numbat)
library(data.table)
library(Matrix)

message(sprintf("[%s] Tissue: %s", Sys.time(), tissue))
message(sprintf("[%s] Loading ATAC matrix...", Sys.time()))
count_mat <- readRDS(atac_file)
message(sprintf("[%s] Dimensions: %d bins × %d cells", Sys.time(), nrow(count_mat), ncol(count_mat)))

message(sprintf("[%s] Loading allele counts...", Sys.time()))
df_allele <- data.table::fread(allele_file)
message(sprintf("[%s] Loaded: %d rows × %d columns", Sys.time(), nrow(df_allele), ncol(df_allele)))

message(sprintf("[%s] Loading ref_hca reference...", Sys.time()))
data(ref_hca)

message(sprintf("[%s] Running NUMBAT analysis (ncores=%d, t=1e-5)...", Sys.time(), 8))
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

message(sprintf("[%s] ✓ NUMBAT analysis complete", Sys.time()))
message(sprintf("[%s] Results saved to: %s", Sys.time(), out_dir))
EOFR

ANALYSIS_EXIT=$?

echo ""
if [[ $ANALYSIS_EXIT -eq 0 ]]; then
  echo "════════════════════════════════════════════════════════════════════════════"
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ ANALYSIS COMPLETE"
  echo "════════════════════════════════════════════════════════════════════════════"
  echo ""
  echo "Output files:"
  ls -lh "$OUTPUT_DIR"
else
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ Analysis failed with exit code $ANALYSIS_EXIT"
  exit $ANALYSIS_EXIT
fi
