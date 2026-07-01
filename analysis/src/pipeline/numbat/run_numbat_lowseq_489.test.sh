#!/bin/bash
# ============================================================================
# TEST: NUMBAT on lowseq_489 - Compute Node Verification
# ============================================================================
# Purpose: Quick test to verify all prerequisites for NUMBAT before full run
# Run on: COMPUTE NODE (use qrsh + tmux, NOT on login node)
# Usage: ./run_numbat_lowseq_489.test.sh
# ============================================================================

set -eo pipefail

echo "[$(date +'%F %T')] [TEST] NUMBAT Prerequisites Check for lowseq_489"
echo "[$(date +'%F %T')] [TEST] Hostname: $(hostname)"

# Must NOT be on login node
if [[ $(hostname) == scc1* ]]; then
  echo "[ERROR] You are on login node! Use qrsh to get compute node first"
  exit 1
fi

# Initialize modules
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

# Load R
echo "[$(date +'%F %T')] [TEST] Loading R module..."
module load R
which Rscript && Rscript --version

# Set project paths
PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
export PROJECT_ROOT

echo ""
echo "[$(date +'%F %T')] [TEST] Checking input files..."

# Check ArchR project
ARCHR_PATH="$PROJECT_ROOT/Data/01_outputs/archR_objects/lowseq_489"
if [[ ! -d "$ARCHR_PATH" ]]; then
  echo "[ERROR] ArchR directory not found: $ARCHR_PATH"
  exit 1
fi
echo "[✓] ArchR project directory exists: $ARCHR_PATH"

# Check barcode file
BARCODE_FILE="$PROJECT_ROOT/Data/01_inputs/barcodes/tissue_barcodes/lowseq_489/lowseq_489.no_edge_effect.barcodes.tsv"
if [[ ! -f "$BARCODE_FILE" ]]; then
  echo "[ERROR] Barcode file not found: $BARCODE_FILE"
  exit 1
fi
BARCODE_COUNT=$(wc -l < "$BARCODE_FILE")
echo "[✓] Barcode file exists: $BARCODE_FILE ($BARCODE_COUNT barcodes)"

# Check output directory writable
OUTPUT_ROOT="$PROJECT_ROOT/Data/04_analysis/cnv/numbat/results/lowseq/run_489"
mkdir -p "$OUTPUT_ROOT"
if ! touch "$OUTPUT_ROOT/.test" 2>/dev/null; then
  echo "[ERROR] Cannot write to output directory: $OUTPUT_ROOT"
  exit 1
fi
rm -f "$OUTPUT_ROOT/.test"
echo "[✓] Output directory writable: $OUTPUT_ROOT"

# Test R loading of ArchR
echo ""
echo "[$(date +'%F %T')] [TEST] Testing ArchR library..."
Rscript -e "
suppressPackageStartupMessages(library(ArchR))
cat('ArchR version:', as.character(packageVersion('ArchR')), '\n')
cat('✓ ArchR loaded successfully\n')
" 2>&1 || { echo "[ERROR] Failed to load ArchR"; exit 1; }

echo ""
echo "[$(date +'%F %T')] [TEST] ✓ All prerequisites verified!"
echo "[$(date +'%F %T')] [TEST] Ready to run NUMBAT pipeline"
echo ""
echo "Next steps:"
echo "  1. Review parameters in copilot-instructions for NUMBAT settings"
echo "  2. Run full NUMBAT: Rscript analysis/src/pipeline/numbat/run_numbat_lowseq_489.R"
echo "  3. Monitor progress and output files in: Data/04_analysis/cnv/numbat/results/lowseq/run_489/"
echo ""

exit 0
