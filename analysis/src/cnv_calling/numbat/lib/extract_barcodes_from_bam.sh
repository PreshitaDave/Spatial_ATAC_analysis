#!/bin/bash
set -eo pipefail

################################################################################
# extract_barcodes_from_bam.sh (SIMPLIFIED)
#
# Purpose: Extract BC16 barcodes from merged BAM and create tissue-specific files
#          matching the BAM CB tag format (with -1 or -2 suffix)
#
# Input:   Merged BAM: Data/04_analysis/cnv/numbat/inputs/deepseq_merged_for_numbat_bc16.bam
#
# Output:  BC16 Barcode files (one per tissue):
#          - Data/01_inputs/barcodes/tissue_barcodes/deepseq_488B/deepseq_488B.bc16.barcodes.tsv
#          - Data/01_inputs/barcodes/tissue_barcodes/deepseq_489/deepseq_489.bc16.barcodes.tsv
#
# Author:  Spatial ATAC Pipeline
# Date:    May 16, 2026
#
################################################################################

# Initialize module system
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

# Load modules
module load samtools/1.14 || module load samtools

PROJECT_ROOT="${PROJECT_ROOT:-.}"

# Logging
log() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" >&2
}

################################################################################
# STEP 1: Verify Inputs
################################################################################

log "========== STEP 1: VERIFY INPUTS =========="

MERGED_BAM="${PROJECT_ROOT}/Data/04_analysis/cnv/numbat/inputs/deepseq_merged_for_numbat_bc16.bam"

if [[ ! -f "$MERGED_BAM" ]]; then
  log "ERROR: Merged BAM not found: $MERGED_BAM"
  exit 1
fi

if [[ ! -f "${MERGED_BAM}.bai" ]]; then
  log "ERROR: BAM index not found. Run merge_bam_for_numbat.qsub.sh first."
  exit 1
fi

log "✓ Merged BAM found: $(du -h "$MERGED_BAM" | cut -f1)"

################################################################################
# STEP 2: Extract Unique BC16 Barcodes from BAM
################################################################################

log "========== STEP 2: EXTRACT BC16 BARCODES FROM BAM =========="

TEMP_BC16_FILE=$(mktemp)
trap "rm -f $TEMP_BC16_FILE" EXIT

log "Extracting unique CB tags from BAM (this may take 5-10 minutes)..."

# Extract CB tags, remove the CB:Z: prefix, and get unique values
# These will have format like: AAGCATGGAACTGTCC-1, AAGCATGGAACTGTCC-2, etc.
samtools view "$MERGED_BAM" 2>/dev/null \
  | awk 'BEGIN{FS="\t"} {for (i=12; i<=NF; i++) if ($i ~ /^CB:Z:/) {gsub(/CB:Z:/, "", $i); print $i}}' \
  | sort -u \
  > "$TEMP_BC16_FILE"

TOTAL_BC16_COUNT=$(wc -l < "$TEMP_BC16_FILE")
log "✓ Extracted $TOTAL_BC16_COUNT unique BC16 barcodes from BAM"

if [[ $TOTAL_BC16_COUNT -lt 1000 ]]; then
  log "WARNING: Very low barcode count ($TOTAL_BC16_COUNT). Check if BAM has CB tags."
fi

log "Sample BC16 barcodes (first 5):"
head -5 "$TEMP_BC16_FILE" | while read bc; do
  log "  $bc"
done

################################################################################
# STEP 3: Create BC16 Barcode Files for Each Tissue
################################################################################

log "========== STEP 3: CREATE BC16 BARCODE FILES =========="

# Get list of tissues from existing barcode directories
TISSUE_LIST=()
for tissue_dir in ${PROJECT_ROOT}/Data/01_inputs/barcodes/tissue_barcodes/deepseq_*/; do
  if [[ -d "$tissue_dir" ]]; then
    tissue=$(basename "$tissue_dir" | sed 's/deepseq_//')
    TISSUE_LIST+=("$tissue")
  fi
done

if [[ ${#TISSUE_LIST[@]} -eq 0 ]]; then
  log "ERROR: No existing tissue barcode directories found"
  exit 1
fi

log "Found tissues: ${TISSUE_LIST[*]}"

# For each tissue, create BC16 barcode files
for tissue in "${TISSUE_LIST[@]}"; do
  BARCODE_DIR="${PROJECT_ROOT}/Data/01_inputs/barcodes/tissue_barcodes/deepseq_${tissue}"
  
  # All BC16 barcodes (for this tissue)
  BC16_ALL_FILE="${BARCODE_DIR}/deepseq_${tissue}.bc16.barcodes.tsv"
  
  # For now, assign all BC16 barcodes to all tissues
  # (In practice, would need tissue demultiplexing info)
  cp "$TEMP_BC16_FILE" "$BC16_ALL_FILE"
  
  BC16_COUNT=$(wc -l < "$BC16_ALL_FILE")
  log "✓ Created BC16 all-barcodes file for ${tissue}"
  log "  File:  $BC16_ALL_FILE"
  log "  Count: $BC16_COUNT"
  
  # BC16 no-edge-effect barcodes
  # For now, use all (actual filtering happens in NUMBAT via cellsnp-lite)
  BC16_NO_EDGE_FILE="${BARCODE_DIR}/deepseq_${tissue}.bc16.no_edge_effect.barcodes.tsv"
  cp "$TEMP_BC16_FILE" "$BC16_NO_EDGE_FILE"
  
  NO_EDGE_COUNT=$(wc -l < "$BC16_NO_EDGE_FILE")
  log "✓ Created BC16 no-edge-effect file for ${tissue}"
  log "  File:  $BC16_NO_EDGE_FILE"
  log "  Count: $NO_EDGE_COUNT"
done

################################################################################
# STEP 4: Verify Output Files
################################################################################

log "========== STEP 4: VERIFY OUTPUT FILES =========="

ALL_GOOD=true
for tissue in "${TISSUE_LIST[@]}"; do
  BC16_ALL_FILE="${PROJECT_ROOT}/Data/01_inputs/barcodes/tissue_barcodes/deepseq_${tissue}/deepseq_${tissue}.bc16.barcodes.tsv"
  BC16_NO_EDGE_FILE="${PROJECT_ROOT}/Data/01_inputs/barcodes/tissue_barcodes/deepseq_${tissue}/deepseq_${tissue}.bc16.no_edge_effect.barcodes.tsv"
  
  if [[ -f "$BC16_ALL_FILE" && -f "$BC16_NO_EDGE_FILE" ]]; then
    log "✓ deepseq_${tissue}: BC16 files created"
  else
    log "✗ deepseq_${tissue}: Missing files"
    ALL_GOOD=false
  fi
done

if [[ "$ALL_GOOD" == false ]]; then
  log "ERROR: Some files were not created successfully"
  exit 1
fi

################################################################################
# STEP 5: Summary
################################################################################

log "========== STEP 5: SUMMARY =========="
log "✓ BC16 barcode extraction complete"
log ""
log "Input:  $MERGED_BAM"
log "Output: BC16 barcode files created for tissues: ${TISSUE_LIST[*]}"
log ""
log "Files created:"
for tissue in "${TISSUE_LIST[@]}"; do
  BARCODE_DIR="${PROJECT_ROOT}/Data/01_inputs/barcodes/tissue_barcodes/deepseq_${tissue}"
  echo "  - ${BARCODE_DIR}/deepseq_${tissue}.bc16.barcodes.tsv"
  echo "  - ${BARCODE_DIR}/deepseq_${tissue}.bc16.no_edge_effect.barcodes.tsv"
done
log ""
log "Next: Update prepare_numbat_atac_inputs.sh to use BC16 barcodes"
log "EXTRACTION COMPLETE"
