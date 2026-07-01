#!/bin/bash
################################################################################
# Generate NUMBAT lambda reference file (lambdas_ATAC_bincnt.rds)
# 
# This creates the aggregated pooled reference for NUMBAT ATAC-bin mode
# using the fixed binning script that handles barcode format normalization
################################################################################

set -eo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
DATASET="lowseq"
TISSUE="489"

# Paths
NUMBAT_EXTDATA="/projectnb/paxlab/presh/Rlibs/4.5/numbat/extdata"

# Input paths
BARCODE_FILE="$PROJECT_ROOT/Data/01_inputs/barcodes/tissue_barcodes/${DATASET}_${TISSUE}/${DATASET}_${TISSUE}.barcodes.tsv"
FRAGMENT_FILE="$PROJECT_ROOT/Data/01_inputs/fragments/${DATASET}_${TISSUE}/${DATASET}_${TISSUE}.fragments.sort.filtered.bed.gz"
BINGR="$NUMBAT_EXTDATA/var220kb.rds"

# Output paths
REF_DIR="$PROJECT_ROOT/Data/02_references"
LAMBDA_OUTPUT="$REF_DIR/lambdas_ATAC_bincnt.rds"

cd "$PROJECT_ROOT"

echo "==============================================="
echo "Generating NUMBAT Lambda Reference File"
echo "==============================================="
echo "Dataset: $DATASET"
echo "Tissue: $TISSUE"
echo "Output: $LAMBDA_OUTPUT"
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

# Load R module
module load R 2>/dev/null || true

# Verify inputs
echo "[STEP 1] Verifying inputs..."
test -f "$BARCODE_FILE" && echo "  ✓ Barcodes: $(basename $BARCODE_FILE)" || { echo "  ✗ Missing: $BARCODE_FILE"; exit 1; }
test -f "$FRAGMENT_FILE" && echo "  ✓ Fragments: $(basename $FRAGMENT_FILE)" || { echo "  ✗ Missing: $FRAGMENT_FILE"; exit 1; }
test -f "$BINGR" && echo "  ✓ Bins: $(basename $BINGR)" || { echo "  ✗ Missing: $BINGR"; exit 1; }

# Create reference directory if needed
mkdir -p "$REF_DIR"
test -w "$REF_DIR" && echo "  ✓ Reference dir writable" || { echo "  ✗ Not writable: $REF_DIR"; exit 1; }

echo ""
echo "[STEP 2] Generating aggregated ATAC reference..."
echo "  Running: get_binned_atac_fixed.R --generateAggRef"
echo ""

Rscript "$PROJECT_ROOT/analysis/src/cnv_calling/numbat/get_binned_atac_fixed.R" \
  --CB "$BARCODE_FILE" \
  --frag "$FRAGMENT_FILE" \
  --binGR "$BINGR" \
  --outFile "$LAMBDA_OUTPUT" \
  --generateAggRef \
  2>&1 | tee "$REF_DIR/lambda_generation.log"

echo ""
echo "[STEP 3] Verification..."
if [[ -f "$LAMBDA_OUTPUT" ]]; then
  SIZE=$(ls -lh "$LAMBDA_OUTPUT" | awk '{print $5}')
  echo "  ✓ Lambda file created: $SIZE"
else
  echo "  ✗ Lambda file NOT created"
  exit 1
fi

echo ""
echo "==============================================="
echo "✓ Lambda reference generation complete!"
echo "==============================================="
