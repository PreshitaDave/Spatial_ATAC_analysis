#!/bin/bash
set -eo pipefail

################################################################################
# Script: regenerate_atac_full_barcodes.qsub.sh
# Purpose: Regenerate ATAC bin-by-cell matrix using FULL barcode list
#          (same as pileup stage) to ensure barcode consistency with allele counts
#
# CRITICAL FIX FOR BARCODE MISMATCH:
# ───────────────────────────────────────────────────────────────────────────
# Problem:  NUMBAT requires exact cell-by-cell correspondence between ATAC 
#           matrix and allele counts. If different barcode files are used:
#           - Pileup stage (allele generation) used: {tissue}.barcodes.tsv (FULL)
#           - ATAC stage previously used: {tissue}.no_edge_effect.barcodes.tsv (FILTERED)
#           Result: 460+ cells mismatch → NUMBAT filters all to 0 coverage → No CNV calls
#
# Solution: Regenerate ATAC matrix with EXACT same barcode file as pileup stage
#           This ensures 100% barcode match between inputs
#
# Usage: qsub regenerate_atac_full_barcodes.qsub.sh tissue_name
#        Example: qsub regenerate_atac_full_barcodes.qsub.sh lowseq_489
#
# Inputs (automatically located for any tissue):
#   - Fragment file: Data/01_inputs/fragments/{tissue}/{tissue}.fragments.sort.filtered.bed.gz
#   - Barcode file:  Data/01_inputs/barcodes/tissue_barcodes/{tissue}/{tissue}.barcodes.tsv
#   - Genomic bins:  Data/04_analysis/cnv/numbat/reference/var220kb.rds (80K, 220kb windows)
#
# Outputs:
#   - New ATAC matrix: Data/04_analysis/cnv/numbat/inputs/{tissue}/atac_bin/{tissue}_atac_bin.rds
#   - Backup of old:   Data/04_analysis/cnv/numbat/inputs/{tissue}/atac_bin/{tissue}_atac_bin_ORIGINAL_BACKUP.rds
#
# Reproducibility Notes:
#   1. This script uses FULL barcode file ({tissue}.barcodes.tsv), not filtered version
#   2. Barcode counts must match pileup stage counts (verify before submission)
#   3. Fragment file must have "-1" or "-2" suffix (10X Cell Ranger format)
#   4. Output size should be >50M (if <1M, barcode matching failed)
#   5. After successful completion, run validation:
#      Rscript analysis/src/cnv_calling/numbat/validate_numbat_inputs.R {tissue}
#   6. Then submit NUMBAT analysis:
#      qsub analysis/src/cnv_calling/numbat/run_numbat_analysis_{tissue}_refhca.qsub.sh
#
# Error Handling:
#   - Backs up original ATAC matrix before regeneration
#   - Verifies output file size (sanity check against corruption)
#   - Checks all input paths exist before starting
#   - Progress indicators for monitoring long-running job
################################################################################

#$ -N regenerate_atac
#$ -cwd
#$ -o analysis/qsub_logs/regenerate_atac_\$JOB_ID.out
#$ -e analysis/qsub_logs/regenerate_atac_\$JOB_ID.err
#$ -pe omp 8
#$ -l h_rt=04:00:00
#$ -l mem_per_core=8G
#$ -P paxlab

echo "════════════════════════════════════════════════════════════════"
echo "ATAC MATRIX REGENERATION WITH FULL BARCODE LIST"
echo "════════════════════════════════════════════════════════════════"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Job ID: $JOB_ID"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Hostname: $(hostname)"

# ────────────────────────────────────────────────────────────────────
# STEP 0: Parse arguments and validate
# ────────────────────────────────────────────────────────────────────
TISSUE="${1:?Error: Must provide tissue name. Usage: qsub script.sh tissue_name}"

# Determine dataset type (lowseq vs deepseq)
if [[ "$TISSUE" =~ ^lowseq ]]; then
  DATASET="lowseq"
elif [[ "$TISSUE" =~ ^deepseq ]]; then
  DATASET="deepseq"
else
  echo "[ERROR] Unknown tissue type. Must start with 'lowseq' or 'deepseq'"
  exit 1
fi

echo "[$(date +'%Y-%m-%d %H:%M:%S')] Tissue: $TISSUE"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Dataset: $DATASET"

# ────────────────────────────────────────────────────────────────────
# STEP 1: Initialize module system (CRITICAL for SGE jobs)
# ────────────────────────────────────────────────────────────────────
echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP 1: Initializing module system..."

set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

module load R
which Rscript && Rscript --version || {
  echo "[ERROR] Failed to load R module"
  exit 1
}

# ────────────────────────────────────────────────────────────────────
# STEP 2: Define paths and verify all inputs exist
# ────────────────────────────────────────────────────────────────────
echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP 2: Verifying input paths..."

FRAG_DIR="Data/01_inputs/fragments/${TISSUE}"
FRAG_FILE="${FRAG_DIR}/${TISSUE}.fragments.sort.filtered.bed.gz"
BARCODE_FILE="Data/01_inputs/barcodes/tissue_barcodes/${TISSUE}/${TISSUE}.barcodes.tsv"
BIN_FILE="Data/04_analysis/cnv/numbat/reference/var220kb.rds"

# Output directories
ATAC_BIN_DIR="Data/04_analysis/cnv/numbat/inputs/${TISSUE}/atac_bin"
OUTPUT_FILE="${ATAC_BIN_DIR}/${TISSUE}_atac_bin.rds"
OUTPUT_FILE_TEMP="${ATAC_BIN_DIR}/${TISSUE}_atac_bin_FULL_BARCODES.rds"
BACKUP_FILE="${ATAC_BIN_DIR}/${TISSUE}_atac_bin_ORIGINAL_BACKUP.rds"

# Verify input files
if [[ ! -f "$FRAG_FILE" ]]; then
  echo "[ERROR] Fragment file not found: $FRAG_FILE"
  exit 1
fi
echo "[✓] Fragment file: $FRAG_FILE ($(stat -c%s "$FRAG_FILE" | numfmt --to=iec))"

if [[ ! -f "$BARCODE_FILE" ]]; then
  echo "[ERROR] Barcode file not found: $BARCODE_FILE"
  exit 1
fi
BARCODE_COUNT=$(wc -l < "$BARCODE_FILE")
echo "[✓] Barcode file: $BARCODE_FILE ($BARCODE_COUNT cells)"

if [[ ! -f "$BIN_FILE" ]]; then
  echo "[ERROR] Genomic bins file not found: $BIN_FILE"
  exit 1
fi
echo "[✓] Genomic bins file: $BIN_FILE ($(stat -c%s "$BIN_FILE" | numfmt --to=iec))"

# Create output directory if it doesn't exist
mkdir -p "$ATAC_BIN_DIR"
echo "[✓] Output directory ready: $ATAC_BIN_DIR"

# ────────────────────────────────────────────────────────────────────
# STEP 3: Backup original ATAC matrix
# ────────────────────────────────────────────────────────────────────
echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP 3: Backing up original ATAC matrix..."

if [[ -f "$OUTPUT_FILE" ]]; then
  cp "$OUTPUT_FILE" "$BACKUP_FILE"
  echo "[✓] Backup saved: $BACKUP_FILE ($(stat -c%s "$BACKUP_FILE" | numfmt --to=iec))"
else
  echo "[⚠] No existing ATAC matrix found at $OUTPUT_FILE (first run)"
fi

# ────────────────────────────────────────────────────────────────────
# STEP 4: Regenerate ATAC bin-by-cell matrix with correct barcodes
# ────────────────────────────────────────────────────────────────────
echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP 4: Regenerating ATAC matrix..."
echo "  - Using barcode file: $BARCODE_FILE ($BARCODE_COUNT cells)"
echo "  - Using fragment file: $FRAG_FILE"
echo "  - Using genomic bins: $BIN_FILE"
echo ""

Rscript analysis/src/cnv_calling/numbat/get_binned_atac_fixed.R \
  --CB "$BARCODE_FILE" \
  --frag "$FRAG_FILE" \
  --binGR "$BIN_FILE" \
  --outFile "$OUTPUT_FILE_TEMP"

# ────────────────────────────────────────────────────────────────────
# STEP 5: Verify regenerated ATAC matrix
# ────────────────────────────────────────────────────────────────────
echo ""
echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP 5: Verifying regenerated ATAC matrix..."

if [[ ! -f "$OUTPUT_FILE_TEMP" ]]; then
  echo "[ERROR] Regeneration failed - output file not created"
  exit 1
fi

OUTPUT_SIZE=$(stat -c%s "$OUTPUT_FILE_TEMP")
OUTPUT_SIZE_MB=$((OUTPUT_SIZE / 1024 / 1024))

echo "[✓] Output file created: $OUTPUT_FILE_TEMP"
echo "[✓] File size: ${OUTPUT_SIZE_MB}M (${OUTPUT_SIZE} bytes)"

# Sanity check: output should be >50M (corrupted outputs are <1M)
if [[ $OUTPUT_SIZE -lt 52428800 ]]; then  # 50M threshold
  echo "[ERROR] Output file suspiciously small (${OUTPUT_SIZE_MB}M < 50M)"
  echo "[ERROR] This likely indicates barcode matching failed"
  echo "[ERROR] Check fragment vs barcode compatibility"
  exit 1
fi

# ────────────────────────────────────────────────────────────────────
# STEP 6: Replace active ATAC matrix with regenerated version
# ────────────────────────────────────────────────────────────────────
echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP 6: Replacing active ATAC matrix..."

mv "$OUTPUT_FILE_TEMP" "$OUTPUT_FILE"
echo "[✓] New ATAC matrix is now active: $OUTPUT_FILE"

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ REGENERATION COMPLETE"
echo "════════════════════════════════════════════════════════════════"

echo ""
echo "NEXT STEPS:"
echo "─────────────────────────────────────────────────────────────────"
echo "1. Validate barcode consistency:"
echo "   Rscript analysis/src/cnv_calling/numbat/validate_numbat_inputs.R $TISSUE"
echo ""
echo "2. If validation PASSES (expect: 100% overlap with ~${BARCODE_COUNT} cells):"
echo "   qsub analysis/src/cnv_calling/numbat/run_numbat_analysis_${TISSUE}_refhca.qsub.sh"
echo ""
echo "3. If validation FAILS:"
echo "   - Check which barcode file was used for pileup stage"
echo "   - Regenerate again with correct barcode file"
echo "   - See: .github/copilot-instructions.md CRITICAL ISSUE 2 section"
echo ""

exit 0
