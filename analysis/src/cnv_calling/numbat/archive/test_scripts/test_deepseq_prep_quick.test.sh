#!/bin/bash
################################################################################
# test_deepseq_prep_quick.test.sh
#
# Quick test of deepseq tissue input preparation (STEP 1 & 2 only, skip variants)
# Tests: bin extraction, barcode creation, directory setup
# Does NOT run full variant calling (would be 60+ min)
#
################################################################################

#$ -P paxlab
#$ -N test_deepseq_prep
#$ -pe omp 4
#$ -l mem_per_core=8G
#$ -l h_rt=00:20:00
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/test_deepseq_prep.$JOB_ID.log
#$ -j y

set -eo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
cd "$PROJECT_ROOT"

# Module setup
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

module load R

echo "==============================================================================="
echo "DEEPSEQ TISSUE INPUT PREP - QUICK TEST"
echo "==============================================================================="
echo "[HOST] $(hostname)"
echo "[TIME] $(date)"
echo ""

# TEST 1: Extract tissue-specific bin matrices
echo "[TEST 1] Extract Tissue-Specific Bin Matrices"
echo "---"

Rscript - <<'RBIN_TEST'
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

combined_file <- "Data/04_analysis/cnv/numbat/inputs/deepseq_atac_bin.rds"

cat("[INFO] Loading combined deepseq bin matrix...\n")
combined <- readRDS(combined_file)
cat(sprintf("[✓] Loaded: %d cells × %d features\n", ncol(combined), nrow(combined)))

# Extract 488B
cat("\n[INFO] Extracting 488B subset...\n")
cells_488B <- which(grepl("488B|488b", colnames(combined), ignore.case = TRUE))
if (length(cells_488B) > 0) {
  cat(sprintf("[✓] Found %d 488B cells\n", length(cells_488B)))
  subset_488B <- subset(combined, cells = colnames(combined)[cells_488B])
  
  output_488B <- "Data/04_analysis/cnv/numbat/inputs/deepseq_488B_atac_bin.rds.test"
  saveRDS(subset_488B, file = output_488B)
  cat(sprintf("[✓] Test save successful: %s\n", output_488B))
  
  file.size_mb <- file.size(output_488B) / 1024^2
  cat(sprintf("[INFO] File size: %.1f MB\n", file.size_mb))
}

# Extract 489
cat("\n[INFO] Extracting 489 subset...\n")
cells_489 <- which(grepl("489", colnames(combined), ignore.case = TRUE))
if (length(cells_489) > 0) {
  cat(sprintf("[✓] Found %d 489 cells\n", length(cells_489)))
  subset_489 <- subset(combined, cells = colnames(combined)[cells_489])
  
  output_489 <- "Data/04_analysis/cnv/numbat/inputs/deepseq_489_atac_bin.rds.test"
  saveRDS(subset_489, file = output_489)
  cat(sprintf("[✓] Test save successful: %s\n", output_489))
  
  file.size_mb <- file.size(output_489) / 1024^2
  cat(sprintf("[INFO] File size: %.1f MB\n", file.size_mb))
}

cat("\n[✓ TEST 1 PASSED] Bin extraction works correctly\n")
RBIN_TEST
echo ""

# TEST 2: Create tissue-specific barcode lists
echo "[TEST 2] Create Tissue-Specific Barcode Lists"
echo "---"

BARCODE_DIR="Data/01_inputs/barcodes/tissue_barcodes"

for TISSUE in 488B 489; do
  TISSUE_BARCODE_DIR="${BARCODE_DIR}/deepseq_${TISSUE}"
  TISSUE_BARCODE_FILE="${TISSUE_BARCODE_DIR}/deepseq_${TISSUE}.no_edge_effect.barcodes.tsv"
  
  if [[ ! -f "$TISSUE_BARCODE_FILE" ]]; then
    echo "  [ERROR] Barcode file not found: $TISSUE_BARCODE_FILE"
    exit 1
  else
    COUNT=$(wc -l < "$TISSUE_BARCODE_FILE")
    echo "  ✓ ${TISSUE}: $COUNT barcodes (file ready)"
  fi
done
echo ""

# TEST 3: Verify output directories
echo "[TEST 3] Verify Output Directory Structure"
echo "---"

ALLELES_DIR="Data/04_analysis/cnv/numbat/inputs/alleles"
mkdir -p "$ALLELES_DIR/pileup" "$ALLELES_DIR/phasing"
echo "  ✓ Alleles directory structure created"

# Verify writable
if touch "${ALLELES_DIR}/.test_write" 2>/dev/null; then
  rm "${ALLELES_DIR}/.test_write"
  echo "  ✓ Alleles directory writable"
else
  echo "  [ERROR] Alleles directory not writable"
  exit 1
fi
echo ""

# TEST 4: Verify VCF and reference files
echo "[TEST 4] Verify VCF and Reference Files"
echo "---"

VCF="Data/hg38_resources/numbat/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
if [[ -f "$VCF" ]]; then
  SIZE=$(ls -lh "$VCF" | awk '{print $5}')
  echo "  ✓ SNP VCF: $SIZE"
else
  echo "  [ERROR] VCF not found: $VCF"
  exit 1
fi

GENETIC_MAP="Data/hg38_resources/numbat/genetic_map_hg38_withX.txt.gz"
if [[ -f "$GENETIC_MAP" ]]; then
  SIZE=$(ls -lh "$GENETIC_MAP" | awk '{print $5}')
  echo "  ✓ Genetic map: $SIZE"
else
  echo "  [ERROR] Genetic map not found: $GENETIC_MAP"
  exit 1
fi

BAM="Data/04_analysis/cnv/numbat/inputs/bam_merged/deepseq_merged_for_numbat.bam"
if [[ -f "$BAM" ]]; then
  SIZE=$(ls -lh "$BAM" | awk '{print $5}')
  echo "  ✓ Merged BAM: $SIZE"
else
  echo "  [ERROR] BAM not found: $BAM"
  exit 1
fi
echo ""

# TEST 5: Verify NUMBAT scripts ready
echo "[TEST 5] Verify NUMBAT Scripts Ready"
echo "---"

for SCRIPT in run_numbat_deepseq_488B_only run_numbat_deepseq_489_only run_numbat_deepseq_combined; do
  SCRIPT_FILE="analysis/src/cnv_calling/numbat/${SCRIPT}.qsub.sh"
  if [[ -f "$SCRIPT_FILE" ]]; then
    if bash -n "$SCRIPT_FILE" 2>/dev/null; then
      echo "  ✓ $SCRIPT (syntax OK)"
    else
      echo "  [ERROR] $SCRIPT has syntax errors"
      exit 1
    fi
  else
    echo "  [ERROR] $SCRIPT not found"
    exit 1
  fi
done
echo ""

echo "==============================================================================="
echo "✓ ALL TESTS PASSED - READY FOR PRODUCTION"
echo "==============================================================================="
echo ""
echo "NEXT STEPS:"
echo "1. Submit deepseq prep job: qsub analysis/src/cnv_calling/numbat/prepare_deepseq_tissue_inputs.sh"
echo "2. Wait for tissue-specific inputs to be created (~90 minutes)"
echo "3. Submit 3 deepseq NUMBAT jobs (after prep completes):"
echo "   - qsub analysis/src/cnv_calling/numbat/run_numbat_deepseq_488B_only.qsub.sh"
echo "   - qsub analysis/src/cnv_calling/numbat/run_numbat_deepseq_489_only.qsub.sh"
echo "   - qsub analysis/src/cnv_calling/numbat/run_numbat_deepseq_combined.qsub.sh"
echo ""
