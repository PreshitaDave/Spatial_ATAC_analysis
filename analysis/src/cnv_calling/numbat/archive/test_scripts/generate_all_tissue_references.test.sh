#!/bin/bash
################################################################################
# Generate NUMBAT References for All Tissues
# 
# Purpose: Generate lambdas_ATAC_bincnt.rds reference files for each tissue
#          using the fixed barcode-compatible binning script
#
# Tissues: lowseq_488B, lowseq_489, deepseq_488B, deepseq_489
# 
# Each reference is generated from that tissue's ATAC data using:
# get_binned_atac_fixed.R --generateAggRef flag
#
# This creates: lambdas_ATAC_bincnt.rds (aggregated reference for CNV detection)
################################################################################

set -eo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
NCORES=8

# Reference directory
REF_DIR="$PROJECT_ROOT/Data/04_analysis/cnv/numbat/reference"
mkdir -p "$REF_DIR"

# Paths to common resources
NUMBAT_EXTDATA="/projectnb/paxlab/presh/Rlibs/4.5/numbat/extdata"
BINGR="$NUMBAT_EXTDATA/var220kb.rds"

cd "$PROJECT_ROOT"

echo "==============================================="
echo "NUMBAT Reference Generation for All Tissues"
echo "==============================================="
echo "Start Time: $(date)"
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

# Array of tissues to process
declare -a DATASETS=("lowseq" "lowseq" "deepseq" "deepseq")
declare -a TISSUES=("488B" "489" "488B" "489")

# Process each tissue
for i in "${!TISSUES[@]}"; do
  DATASET="${DATASETS[$i]}"
  TISSUE="${TISSUES[$i]}"
  DATASET_TISSUE="${DATASET}_${TISSUE}"
  
  echo ""
  echo "=========== Processing: $DATASET_TISSUE ==========="
  
  # Input paths
  BARCODE_FILE="$PROJECT_ROOT/Data/01_inputs/barcodes/tissue_barcodes/${DATASET_TISSUE}/${DATASET_TISSUE}.barcodes.tsv"
  FRAGMENT_FILE="$PROJECT_ROOT/Data/01_inputs/fragments/${DATASET_TISSUE}/${DATASET_TISSUE}.fragments.sort.filtered.bed.gz"
  LAMBDA_OUTPUT="$REF_DIR/lambdas_${DATASET_TISSUE}_ATAC_bincnt.rds"
  
  # Check inputs
  if [[ ! -f "$BARCODE_FILE" ]]; then
    echo "⚠ Skipping $DATASET_TISSUE: barcode file not found"
    continue
  fi
  
  if [[ ! -f "$FRAGMENT_FILE" ]]; then
    echo "⚠ Skipping $DATASET_TISSUE: fragment file not found"
    continue
  fi
  
  echo "[STEP 1] Verifying inputs for $DATASET_TISSUE..."
  echo "  ✓ Barcodes: $(basename $BARCODE_FILE)"
  echo "  ✓ Fragments: $(basename $FRAGMENT_FILE)"
  echo "  ✓ Bins: $(basename $BINGR)"
  
  echo ""
  echo "[STEP 2] Generating aggregated reference for $DATASET_TISSUE..."
  echo "  Output: $(basename $LAMBDA_OUTPUT)"
  echo ""
  
  # Run binning with --generateAggRef flag
  Rscript "$PROJECT_ROOT/analysis/src/cnv_calling/numbat/get_binned_atac_fixed.R" \
    --CB "$BARCODE_FILE" \
    --frag "$FRAGMENT_FILE" \
    --binGR "$BINGR" \
    --outFile "$LAMBDA_OUTPUT" \
    --generateAggRef \
    2>&1 | tee "$REF_DIR/${DATASET_TISSUE}_reference_generation.log"
  
  # Verify output
  if [[ -f "$LAMBDA_OUTPUT" ]]; then
    SIZE=$(ls -lh "$LAMBDA_OUTPUT" | awk '{print $5}')
    echo "  ✓ Reference created: $SIZE"
  else
    echo "  ✗ FAILED to create reference"
  fi
  
  echo ""
done

echo "==============================================="
echo "✓ Reference generation complete!"
echo "  Completed: $(date)"
echo "==============================================="
echo ""
echo "Generated references:"
ls -lh "$REF_DIR"/lambdas_*ATAC_bincnt.rds 2>/dev/null || echo "No references found"
