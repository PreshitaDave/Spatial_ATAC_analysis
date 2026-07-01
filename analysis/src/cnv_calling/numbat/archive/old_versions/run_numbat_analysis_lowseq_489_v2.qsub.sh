#!/bin/bash
#
# run_numbat_analysis_lowseq_489_v2.qsub.sh
# Purpose: NUMBAT refhca analysis for lowseq_489 (FIXED VERSION)
# Calls standalone R script to avoid heredoc issues
#

#$ -N numbat_489_v2
#$ -l h_rt=08:00:00
#$ -pe omp 8
#$ -P paxlab
#$ -l mem_per_core=8G
#$ -o analysis/qsub_logs/numbat_489_v2_$JOB_ID.log

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
echo "NUMBAT refhca Analysis: lowseq_489 (FIXED VERSION WITH CORRECT BARCODES)"
echo "════════════════════════════════════════════════════════════════════════════"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Job ID: $JOB_ID"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Hostname: $(hostname)"
echo ""

# Configuration
TISSUE="lowseq_489"
NCORES=8
ATAC_BIN="Data/04_analysis/cnv/numbat/inputs/${TISSUE}/atac_bin/${TISSUE}_atac_bin.rds"
ALLELE_DF="Data/04_analysis/cnv/numbat/inputs/${TISSUE}/alleles/${TISSUE}_atac_allele_counts.tsv.gz"
OUTPUT_DIR="Data/04_analysis/cnv/numbat/results/${TISSUE}/refhca_run_20260519_v2/"
R_SCRIPT="analysis/src/cnv_calling/numbat/run_numbat_refhca.R"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] Configuration:"
echo "[$(date +'%Y-%m-%d %H:%M:%S')]   Tissue: $TISSUE"
echo "[$(date +'%Y-%m-%d %H:%M:%S')]   NCORES: $NCORES"
echo "[$(date +'%Y-%m-%d %H:%M:%S')]   ATAC: $ATAC_BIN"
echo "[$(date +'%Y-%m-%d %H:%M:%S')]   Alleles: $ALLELE_DF"
echo "[$(date +'%Y-%m-%d %H:%M:%S')]   Output: $OUTPUT_DIR"
echo "[$(date +'%Y-%m-%d %H:%M:%S')]   R Script: $R_SCRIPT"
echo ""

# Verify inputs exist
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Pre-flight checks..."
for f in "$ATAC_BIN" "$ALLELE_DF" "$R_SCRIPT"; do
  if [[ ! -f "$f" ]]; then
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ ERROR: File not found: $f"
    exit 1
  fi
done
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ All input files verified"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ Output directory created: $OUTPUT_DIR"
echo ""

# Run NUMBAT via R script
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ──────────────────────────────────────────────────────"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Calling: Rscript $R_SCRIPT"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ──────────────────────────────────────────────────────"
echo ""

Rscript "$R_SCRIPT" "$TISSUE" "$ATAC_BIN" "$ALLELE_DF" "$OUTPUT_DIR" "$NCORES"

EXIT_CODE=$?

echo ""
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ──────────────────────────────────────────────────────"
if [[ $EXIT_CODE -eq 0 ]]; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ Job completed successfully (exit code: 0)"
else
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ Job failed with exit code: $EXIT_CODE"
fi
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ──────────────────────────────────────────────────────"
echo ""

exit $EXIT_CODE
