#!/bin/bash
################################################################################
# run_numbat_analysis_atac.sh
#
# PURPOSE: Run NUMBAT analysis in ATAC-bin mode
# REQUIRES: Prepared inputs from prepare_numbat_inputs.sh
#
# INPUTS:
#   - ATAC bin matrix (RDS)
#   - Allele counts (TSV.GZ)
#   - ATAC reference (RDS)
#
# USAGE: bash run_numbat_analysis_atac.sh <dataset> <tissue>
#
################################################################################

set -eo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
DATASET="${1:-lowseq}"
TISSUE="${2:-488B}"
NCORES=8

NUMBAT_BIN="/projectnb/paxlab/presh/Rlibs/4.5/numbat/bin"
NUMBAT_EXTDATA="$NUMBAT_BIN/../extdata"

# Input paths
ATAC_BIN="Data/04_analysis/cnv/numbat/inputs/${DATASET}_${TISSUE}_atac_bin.rds"
ALLELE_DF="Data/04_analysis/cnv/numbat/inputs/${DATASET}_${TISSUE}_comb_allele_counts.tsv.gz"
ATAC_REF="Data/02_references/genome/hg38_resources/numbat/lambdas_ATAC_bincnt.rds"

# Output path
OUTPUT_DIR="Data/04_analysis/cnv/numbat/results/${DATASET}_${TISSUE}_atac"

# Parameters file (optional)
PAR_FILE="$NUMBAT_EXTDATA/par_numbatm.rds"
BINGR="$NUMBAT_EXTDATA/var220kb.rds"

cd "$PROJECT_ROOT"

echo "==============================================="
echo "NUMBAT ATAC-BIN ANALYSIS: ${DATASET} ${TISSUE}"
echo "==============================================="
echo ""

# Verify inputs
echo "[STEP 1] Verifying inputs..."
for input in "$ATAC_BIN" "$ALLELE_DF" "$ATAC_REF" "$BINGR"; do
  if [[ ! -f "$input" ]]; then
    echo "[ERROR] Missing: $input"
    exit 1
  fi
  echo "  ✓ $(basename $input)"
done

mkdir -p "$OUTPUT_DIR"

# Create parameters R script
echo "[STEP 2] Creating parameters..."
cat > /tmp/par_numbatm_local.R <<'EOF'
# NUMBAT parameters for ATAC-bin mode
par_numbatm <- list(
  min_cells = 5,           # Minimum cells per clone
  max_iters = 100,         # Max EM iterations
  ncores = 8,              # Parallel cores
  diploid_chrs = "Y",      # Known diploid chromosomes
  verbose = TRUE
)
EOF

echo "[STEP 3] Running NUMBAT analysis in ATAC-bin mode..."
module load R 2>/dev/null || true

Rscript "$NUMBAT_BIN/run_numbat_multiome.R" \
  --countmat "$ATAC_BIN" \
  --alleledf "$ALLELE_DF" \
  --out_dir "$OUTPUT_DIR" \
  --ref "$ATAC_REF" \
  --gtf "$BINGR" \
  --mode ATAC-bin \
  --ncores $NCORES \
  2>&1 | tee "$OUTPUT_DIR/numbat_${DATASET}_${TISSUE}.log"

echo ""
echo "[STEP 4] Checking results..."
if [[ -f "$OUTPUT_DIR/cnv_calls.rds" ]]; then
  echo "  ✓ CNV calls saved"
else
  echo "  ⚠ CNV calls not found (check log above)"
fi

if [[ -f "$OUTPUT_DIR/phylogeny.png" ]]; then
  echo "  ✓ Phylogeny plot saved"
fi

echo ""
echo "[✓] Analysis complete: $OUTPUT_DIR"
