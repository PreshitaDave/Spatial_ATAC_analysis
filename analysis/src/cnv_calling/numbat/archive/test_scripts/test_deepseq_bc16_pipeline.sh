#!/bin/bash
set -eo pipefail

################################################################################
# test_deepseq_bc16_pipeline.sh
#
# Quick validation test to verify all setup before running actual NUMBAT jobs
# Run this on a compute node (via tmux + qrsh) BEFORE submitting merge/extract jobs
#
# Usage:
#   1. tmux a -t spatial_atac_work (attach to existing tmux session)
#   2. qrsh -l h_rt=4:00:00 -pe omp 4 -P paxlab -l mem_per_core=8G
#   3. cd /projectnb/paxlab/presh/projects/spatial_atac
#   4. bash analysis/src/cnv_calling/numbat/test_deepseq_bc16_pipeline.sh
#
################################################################################

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log_pass() {
  echo -e "${GREEN}✓${NC} $*"
}

log_fail() {
  echo -e "${RED}✗${NC} $*"
}

log_warn() {
  echo -e "${YELLOW}⚠${NC} $*"
}

log_info() {
  echo "[INFO] $*"
}

################################################################################
# STEP 1: Verify Compute Node
################################################################################

log_info "========== STEP 1: VERIFY COMPUTE NODE =========="

HOSTNAME=$(hostname)
if [[ "$HOSTNAME" == scc1* ]]; then
  log_fail "Running on login node ($HOSTNAME)"
  log_info "Use: qrsh -l h_rt=4:00:00 -pe omp 4 -P paxlab -l mem_per_core=8G"
  exit 1
fi
log_pass "Running on compute node: $HOSTNAME"

################################################################################
# STEP 2: Check Module System
################################################################################

log_info "========== STEP 2: CHECK MODULE SYSTEM =========="

if ! command -v module >/dev/null 2>&1; then
  log_fail "Module system not available"
  exit 1
fi
log_pass "Module system available"

# Try loading samtools
if module load samtools/1.14 2>/dev/null || module load samtools 2>/dev/null; then
  log_pass "samtools module loaded"
  which samtools || log_warn "samtools not in PATH after module load"
else
  log_fail "Cannot load samtools module"
  exit 1
fi

################################################################################
# STEP 3: Verify Source BAMs
################################################################################

log_info "========== STEP 3: VERIFY SOURCE BAMS =========="

BAM_SOURCE_DIR="${PROJECT_ROOT}/Data/04_analysis/03_intermediate/variant_calling/monopogen/variant_calling/deepseq/Bam"

if [[ ! -d "$BAM_SOURCE_DIR" ]]; then
  log_fail "BAM source directory not found: $BAM_SOURCE_DIR"
  exit 1
fi
log_pass "BAM source directory exists: $BAM_SOURCE_DIR"

# Count BAMs
BAM_COUNT=$(find "$BAM_SOURCE_DIR" -name "deepseq_chr*.filter.bam" -type f | wc -l)
log_info "Found $BAM_COUNT chromosome BAMs"

if [[ $BAM_COUNT -eq 22 ]]; then
  log_pass "All 22 chromosome BAMs present"
else
  log_fail "Expected 22 BAMs, found $BAM_COUNT"
  exit 1
fi

# Check sizes
log_info "Sample BAM sizes:"
for chr in 1 10 22; do
  if [[ -f "$BAM_SOURCE_DIR/deepseq_chr${chr}.filter.bam" ]]; then
    size=$(du -h "$BAM_SOURCE_DIR/deepseq_chr${chr}.filter.bam" | cut -f1)
    log_info "  chr${chr}: $size"
  fi
done

################################################################################
# STEP 4: Verify Output Directories
################################################################################

log_info "========== STEP 4: VERIFY OUTPUT DIRECTORIES =========="

OUTPUT_DIR="${PROJECT_ROOT}/Data/04_analysis/cnv/numbat/inputs"
if [[ ! -d "$OUTPUT_DIR" ]]; then
  log_fail "Output directory doesn't exist: $OUTPUT_DIR"
  log_info "Creating: $OUTPUT_DIR"
  mkdir -p "$OUTPUT_DIR" || exit 1
fi
log_pass "Output directory writable: $OUTPUT_DIR"

# Test write permission
TEST_FILE="${OUTPUT_DIR}/.test_write_$$"
if touch "$TEST_FILE" && rm "$TEST_FILE"; then
  log_pass "Output directory is writable"
else
  log_fail "Output directory is NOT writable"
  exit 1
fi

BARCODE_DIRS=(
  "${PROJECT_ROOT}/Data/01_inputs/barcodes/tissue_barcodes/deepseq_488B"
  "${PROJECT_ROOT}/Data/01_inputs/barcodes/tissue_barcodes/deepseq_489"
)

for dir in "${BARCODE_DIRS[@]}"; do
  if [[ ! -d "$dir" ]]; then
    log_fail "Barcode directory not found: $dir"
    exit 1
  fi
  
  if ! touch "$dir/.test_write_$$" 2>/dev/null; then
    log_fail "Barcode directory not writable: $dir"
    exit 1
  fi
  rm "$dir/.test_write_$$"
  log_pass "Barcode directory writable: $dir"
done

################################################################################
# STEP 5: Verify Existing Barcode Files
################################################################################

log_info "========== STEP 5: VERIFY EXISTING BARCODE FILES =========="

for tissue in 488B 489; do
  barcode_file="${PROJECT_ROOT}/Data/01_inputs/barcodes/tissue_barcodes/deepseq_${tissue}/deepseq_${tissue}.no_edge_effect.barcodes.tsv"
  if [[ -f "$barcode_file" ]]; then
    count=$(wc -l < "$barcode_file")
    log_pass "deepseq_${tissue}: $count BC8 barcodes"
  else
    log_fail "BC8 barcode file not found: $barcode_file"
  fi
done

################################################################################
# STEP 6: Verify Scripts Exist
################################################################################

log_info "========== STEP 6: VERIFY SCRIPTS EXIST =========="

SCRIPTS=(
  "analysis/src/cnv_calling/numbat/merge_bam_for_numbat.qsub.sh"
  "analysis/src/cnv_calling/numbat/extract_barcodes_from_bam.sh"
  "analysis/src/numbat/numbat/prepare_numbat_atac_inputs.sh"
)

for script in "${SCRIPTS[@]}"; do
  script_path="${PROJECT_ROOT}/$script"
  if [[ -f "$script_path" ]]; then
    log_pass "Script exists: $script"
    
    # Test syntax
    if bash -n "$script_path" 2>/dev/null; then
      log_pass "  ✓ Syntax OK"
    else
      log_fail "  ✗ Syntax error in $script"
      bash -n "$script_path" || true
      exit 1
    fi
  else
    log_fail "Script not found: $script"
    exit 1
  fi
done

################################################################################
# STEP 7: Test BAM Merge Script Dry-Run
################################################################################

log_info "========== STEP 7: TEST MERGE LOGIC =========="

# Just verify the first BAM can be read
TEST_BAM="${BAM_SOURCE_DIR}/deepseq_chr1.filter.bam"
if samtools view -H "$TEST_BAM" >/dev/null 2>&1; then
  log_pass "Can read BAM header: deepseq_chr1.filter.bam"
else
  log_fail "Cannot read BAM: deepseq_chr1.filter.bam"
  exit 1
fi

# Try sampling a read
if samtools view "$TEST_BAM" 2>/dev/null | head -1 | grep -q ""; then
  log_pass "Can read BAM records"
else
  log_warn "Cannot read BAM records (might be empty or compression issue)"
fi

################################################################################
# STEP 8: Summary and Next Steps
################################################################################

log_info "========== STEP 8: SUMMARY =========="
log_pass "All pre-flight checks passed! ✓"
log_info ""
log_info "Next steps:"
log_info "1. Submit merge job:"
log_info "   cd ${PROJECT_ROOT}"
log_info "   qsub analysis/src/cnv_calling/numbat/merge_bam_for_numbat.qsub.sh"
log_info ""
log_info "2. Monitor job:"
log_info "   qstat -u preshita | grep merge"
log_info "   tail -f analysis/qsub_logs/merge_bam_numbat.out"
log_info ""
log_info "3. After merge completes, extract barcodes:"
log_info "   bash analysis/src/cnv_calling/numbat/extract_barcodes_from_bam.sh"
log_info ""
log_info "4. Then resubmit deepseq jobs:"
log_info "   qsub analysis/src/cnv_calling/numbat/numbat/run_numbat_deepseq_488B.qsub.sh"
log_info "   qsub analysis/src/cnv_calling/numbat/numbat/run_numbat_deepseq_489.qsub.sh"
log_info ""

log_pass "TEST COMPLETE - Ready to proceed!"
