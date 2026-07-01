#!/bin/bash
#
# run_numbat_analysis_lowseq_488B_refhca_fixed.qsub.sh
# Purpose: NUMBAT refhca analysis for lowseq_488B with CORRECTED barcode matching
# Prerequisites: Regenerated ATAC matrix using correct pileup barcodes
#

#$ -N numbat_refhca_488B_fixed
#$ -l h_rt=08:00:00
#$ -pe omp 8
#$ -P paxlab
#$ -l mem_per_core=8G
#$ -o analysis/qsub_logs/numbat_refhca_488B_fixed_$JOB_ID.log

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

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
cd "$PROJECT_ROOT"

echo "════════════════════════════════════════════════════════════════════════════"
echo "NUMBAT refhca Analysis: lowseq_488B (with corrected barcodes)"
echo "════════════════════════════════════════════════════════════════════════════"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Job ID: $JOB_ID"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Hostname: $(hostname)"
echo ""

# Configuration
TISSUE="lowseq_488B"
NCORES=8
ATAC_BIN="Data/04_analysis/cnv/numbat/inputs/${TISSUE}/atac_bin/${TISSUE}_atac_bin.rds"
ALLELE_DF="Data/04_analysis/cnv/numbat/inputs/${TISSUE}/alleles/${TISSUE}_atac_allele_counts.tsv.gz"
OUTPUT_DIR="Data/04_analysis/cnv/numbat/results/${TISSUE}/refhca_run_20260519_fixed/"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] TISSUE: $TISSUE"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] NCORES: $NCORES"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ATAC: $ATAC_BIN"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Alleles: $ALLELE_DF"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Output: $OUTPUT_DIR"
echo ""

# Verify inputs exist
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Verifying input files..."
if [[ ! -f "$ATAC_BIN" ]]; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ ERROR: ATAC file not found: $ATAC_BIN"
  exit 1
fi
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ ATAC file: $(ls -lh "$ATAC_BIN" | awk '{print $5}')"

if [[ ! -f "$ALLELE_DF" ]]; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ ERROR: Allele file not found: $ALLELE_DF"
  exit 1
fi
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ Allele file: $(ls -lh "$ALLELE_DF" | awk '{print $5}')"

# Create output directory
mkdir -p "$OUTPUT_DIR"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ Output directory created"
echo ""

# Run NUMBAT analysis
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Starting NUMBAT analysis..."
echo ""

Rscript --vanilla - "$TISSUE" "$ATAC_BIN" "$ALLELE_DF" "$OUTPUT_DIR" "$NCORES" <<'EOFR'
args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
atac_file <- args[2]
allele_file <- args[3]
out_dir <- args[4]
ncores <- as.numeric(args[5])

message(sprintf("[%s] Loading libraries...", Sys.time()))
suppressPackageStartupMessages({
  library(numbat)
  library(data.table)
  library(Matrix)
})

message(sprintf("[%s] Loading ATAC matrix...", Sys.time()))
count_mat <- readRDS(atac_file)
message(sprintf("[%s] ATAC dimensions: %d bins x %d cells", Sys.time(), nrow(count_mat), ncol(count_mat)))

message(sprintf("[%s] Loading allele counts...", Sys.time()))
df_allele <- data.table::fread(allele_file)
message(sprintf("[%s] Allele counts: %d rows (variants x cells)", Sys.time(), nrow(df_allele)))

message(sprintf("[%s] Loading reference...", Sys.time()))
data(ref_hca)

message(sprintf("[%s] Running NUMBAT analysis...", Sys.time()))
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

message(sprintf("[%s] NUMBAT analysis complete!", Sys.time()))
message(sprintf("[%s] Results saved to: %s", Sys.time(), out_dir))
EOFR

echo ""
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Analysis completed successfully"
echo "════════════════════════════════════════════════════════════════════════════"
