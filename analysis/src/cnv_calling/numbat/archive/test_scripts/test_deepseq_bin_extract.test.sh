#!/bin/bash
#$ -P paxlab
#$ -N test_deepseq_bin_extract
#$ -pe omp 4
#$ -l mem_per_core=8G
#$ -l h_rt=00:30:00
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/test_deepseq_bin_extract.$JOB_ID.log
#$ -j n

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

echo "[$(date +'%F %T')] TEST: Deepseq Bin Extraction"
echo "[NODE] $(hostname)"
echo ""

# Test 1: Load combined bin matrix
echo "[TEST 1] Loading combined deepseq bin matrix..."
Rscript - <<'RTEST1' "$PROJECT_ROOT"
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
args <- commandArgs(trailingOnly = TRUE)
project_root <- args[1]

combined_file <- file.path(project_root, "Data/04_analysis/cnv/numbat/inputs/deepseq_atac_bin.rds")
cat(sprintf("[INFO] File: %s\n", combined_file))
cat(sprintf("[INFO] Size: %.1f MB\n", file.size(combined_file) / 1024^2))

combined <- readRDS(combined_file)
cat(sprintf("[✓] Loaded: %d cells × %d features\n", ncol(combined), nrow(combined)))
cat(sprintf("[✓] Metadata columns: %s\n", paste(colnames(combined@meta.data), collapse=", ")))

# Check for tissue info
if ("tissue" %in% colnames(combined@meta.data)) {
  cat("[✓] Tissue column exists\n")
  tissues <- table(combined@meta.data$tissue)
  print(tissues)
} else {
  cat("[WARN] No tissue column - will extract by cell name pattern\n")
}

# Test extraction for 488B
cat("\n[TEST 2] Extracting 488B subset...\n")
cells_488B <- which(grepl("488B|488b", colnames(combined), ignore.case=TRUE))
if (length(cells_488B) > 0) {
  cat(sprintf("[✓] Found %d 488B cells\n", length(cells_488B)))
  subset_488B <- subset(combined, cells = colnames(combined)[cells_488B])
  cat(sprintf("[✓] Subset created: %d cells\n", ncol(subset_488B)))
}

# Test extraction for 489
cat("\n[TEST 3] Extracting 489 subset...\n")
cells_489 <- which(grepl("489", colnames(combined), ignore.case=TRUE))
if (length(cells_489) > 0) {
  cat(sprintf("[✓] Found %d 489 cells\n", length(cells_489)))
  subset_489 <- subset(combined, cells = colnames(combined)[cells_489])
  cat(sprintf("[✓] Subset created: %d cells\n", ncol(subset_489)))
}

cat("\n[✓ TESTS PASSED] Ready for tissue extraction\n")
RTEST1

echo ""
echo "[$(date +'%F %T')] TEST COMPLETE"
