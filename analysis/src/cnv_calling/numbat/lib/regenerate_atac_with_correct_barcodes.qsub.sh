#!/bin/bash
#
# regenerate_atac_with_correct_barcodes.qsub.sh
# Purpose: Regenerate ATAC bin matrix using the CORRECT barcode file 
#          (the one used for pileup/allele counts)
#
# Usage: qsub regenerate_atac_with_correct_barcodes.qsub.sh
#

#$ -N regenerate_atac_lowseq_489
#$ -l h_rt=04:00:00
#$ -pe omp 8
#$ -P paxlab
#$ -l mem_per_core=8G
#$ -o analysis/qsub_logs/regenerate_atac_lowseq_489_$JOB_ID.out
#$ -e analysis/qsub_logs/regenerate_atac_lowseq_489_$JOB_ID.err

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

echo "════════════════════════════════════════════════════════════════"
echo "REGENERATING ATAC MATRIX WITH CORRECT BARCODES"
echo "════════════════════════════════════════════════════════════════"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Job ID: $JOB_ID"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Hostname: $(hostname)"
echo ""

# ════════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ════════════════════════════════════════════════════════════════════════════════
tissue="lowseq_489"
dataset="lowseq"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] Tissue: $tissue"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Dataset: $dataset"
echo ""

# ════════════════════════════════════════════════════════════════════════════════
# PATHS
# ════════════════════════════════════════════════════════════════════════════════
echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP 1: Verifying all input paths..."

FRAG_DIR="Data/01_inputs/fragments/${tissue}/"
FRAG_FILE="${FRAG_DIR}${tissue}.fragments.sort.filtered.bed.gz"
CORRECT_BC_FILE="Data/04_analysis/cnv/numbat/inputs/${tissue}/barcodes/${tissue}_atac_barcodes_for_pileup.tsv"
BINS_FILE="Data/04_analysis/cnv/numbat/reference/var220kb.rds"
OUT_DIR="Data/04_analysis/cnv/numbat/inputs/${tissue}/atac_bin/"
OUT_FILE="${OUT_DIR}${tissue}_atac_bin_CORRECTED.rds"
BACKUP_FILE="${OUT_DIR}${tissue}_atac_bin_ORIGINAL_BACKUP.rds"

# Verify fragment file exists
if [[ ! -f "$FRAG_FILE" ]]; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ ERROR: Fragment file not found: $FRAG_FILE"
  exit 1
fi
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ Fragment file found: $FRAG_FILE"

# Verify barcode file exists
if [[ ! -f "$CORRECT_BC_FILE" ]]; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ ERROR: Barcode file not found: $CORRECT_BC_FILE"
  exit 1
fi
BC_COUNT=$(wc -l < "$CORRECT_BC_FILE")
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ Barcode file found with $BC_COUNT cells: $CORRECT_BC_FILE"

# Verify bins file exists
if [[ ! -f "$BINS_FILE" ]]; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ ERROR: Bins file not found: $BINS_FILE"
  exit 1
fi
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ Genomic bins file found: $BINS_FILE"

# Create output directory
mkdir -p "$OUT_DIR"

# Backup original ATAC if it exists
if [[ -f "${OUT_DIR}${tissue}_atac_bin.rds" ]]; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] Backing up original ATAC matrix..."
  cp "${OUT_DIR}${tissue}_atac_bin.rds" "$BACKUP_FILE"
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ Backup saved: $BACKUP_FILE"
fi

echo ""
echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP 2: Regenerating ATAC matrix with corrected barcodes..."
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Running: Rscript get_binned_atac_fixed.R"
echo ""

# ════════════════════════════════════════════════════════════════════════════════
# RUN REGENERATION
# ════════════════════════════════════════════════════════════════════════════════
Rscript analysis/src/cnv_calling/numbat/get_binned_atac_fixed.R \
  --CB "$CORRECT_BC_FILE" \
  --frag "$FRAG_FILE" \
  --binGR "$BINS_FILE" \
  --outFile "$OUT_FILE"

echo ""
echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP 3: Verifying corrected ATAC matrix..."

if [[ ! -f "$OUT_FILE" ]]; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ ERROR: Output file not created: $OUT_FILE"
  exit 1
fi

FILE_SIZE=$(ls -lh "$OUT_FILE" | awk '{print $5}')
FILE_SIZE_BYTES=$(stat --format=%s "$OUT_FILE")
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ Output file created: $OUT_FILE"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ File size: $FILE_SIZE ($FILE_SIZE_BYTES bytes)"

# Verify file is valid (at least 10MB)
if [[ $FILE_SIZE_BYTES -lt 10000000 ]]; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ⚠ WARNING: File seems small. Check if regeneration worked correctly."
  ls -lh "$OUT_FILE"
fi

echo ""
echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP 4: Replacing original with corrected matrix..."
mv "$OUT_FILE" "${OUT_DIR}${tissue}_atac_bin.rds"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ Corrected matrix now active"
echo ""

echo "════════════════════════════════════════════════════════════════"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ REGENERATION COMPLETE"
echo "════════════════════════════════════════════════════════════════"
echo ""
echo "Next steps:"
echo "  1. Run validation: Rscript analyze_numbat_inputs.R lowseq_489"
echo "  2. Then resubmit NUMBAT analysis: qsub analysis/src/cnv_calling/numbat/run_numbat_analysis_lowseq_489_refhca.qsub.sh"
