#!/bin/bash
################################################################################
# run_numbat_analysis_lowseq_489_refhca_simple.qsub.sh
#
# PURPOSE: Run NUMBAT CNV analysis for lowseq_489 using ref_hca
#          Simplified version WITHOUT validation wrapper (validation causes silent failure)
#
################################################################################

#$ -N numbat_refhca_489_simple
#$ -l h_rt=08:00:00
#$ -pe omp 8
#$ -P paxlab
#$ -l mem_per_core=8G
#$ -o analysis/qsub_logs/numbat_refhca_489_simple_$JOB_ID.log

set -eo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
TISSUE="lowseq_489"
NCORES=8

cd "$PROJECT_ROOT"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] ════════════════════════════════════════════"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] NUMBAT refhca: $TISSUE"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Job ID: $JOB_ID"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Hostname: $(hostname)"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ════════════════════════════════════════════"

# Initialize modules
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

module load R 2>/dev/null || true
echo "[$(date +'%Y-%m-%d %H:%M:%S')] R: $(which Rscript)"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Version: $(Rscript --version 2>&1)"

# Prepare paths
ATAC_BIN="$PROJECT_ROOT/Data/04_analysis/cnv/numbat/inputs/${TISSUE}/atac_bin/${TISSUE}_atac_bin.rds"
ALLELE_DF="$PROJECT_ROOT/Data/04_analysis/cnv/numbat/inputs/${TISSUE}/alleles/${TISSUE}_atac_allele_counts.tsv.gz"
OUTPUT_DIR="$PROJECT_ROOT/Data/04_analysis/cnv/numbat/results/${TISSUE}/refhca_run_20260519/"

mkdir -p "$OUTPUT_DIR"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] ATAC: $ATAC_BIN"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Alleles: $ALLELE_DF"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Output: $OUTPUT_DIR"

# Verify files
if [[ ! -f "$ATAC_BIN" ]]; then
  echo "[ERROR] ATAC file not found: $ATAC_BIN" >&2
  exit 1
fi

if [[ ! -f "$ALLELE_DF" ]]; then
  echo "[ERROR] Allele file not found: $ALLELE_DF" >&2
  exit 1
fi

echo "[$(date +'%Y-%m-%d %H:%M:%S')] Starting NUMBAT analysis..."
echo ""

# Run NUMBAT - inline R script
Rscript - "$TISSUE" "$ATAC_BIN" "$ALLELE_DF" "$OUTPUT_DIR" <<'EOFR'
args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
atac_file <- args[2]
allele_file <- args[3]
out_dir <- args[4]

message(sprintf("[%s] Loading libraries...", Sys.time()))
library(numbat)
library(data.table)
library(Matrix)

message(sprintf("[%s] Tissue: %s", Sys.time(), tissue))
message(sprintf("[%s] Loading ATAC matrix: %s", Sys.time(), atac_file))
count_mat <- readRDS(atac_file)
message(sprintf("[%s]   Dimensions: %d bins × %d cells", Sys.time(), nrow(count_mat), ncol(count_mat)))

message(sprintf("[%s] Loading allele counts: %s", Sys.time(), allele_file))
df_allele <- data.table::fread(allele_file)
message(sprintf("[%s]   Loaded: %d rows", Sys.time(), nrow(df_allele)))

message(sprintf("[%s] Loading ref_hca...", Sys.time()))
data(ref_hca)

message(sprintf("[%s] Running NUMBAT analysis with ref_hca...", Sys.time()))
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

message(sprintf("[%s] ✓ Analysis complete", Sys.time()))
message(sprintf("[%s] Results in: %s", Sys.time(), out_dir))
EOFR

EXIT_CODE=$?

echo ""
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Job exit code: $EXIT_CODE"

if [[ $EXIT_CODE -eq 0 ]]; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ SUCCESS"
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] Output files:"
  ls -lh "$OUTPUT_DIR" | head -20
else
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ FAILED with exit code $EXIT_CODE"
  exit $EXIT_CODE
fi
