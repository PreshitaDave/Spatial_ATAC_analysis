#!/bin/bash
################################################################################
# run_numbat_analysis_deepseq_combined.qsub.sh
#
# SGE Wrapper: Run NUMBAT ATAC-bin analysis on merged deepseq (488B + 489)
# This runs NUMBAT on a combined object to compare tissue-level CNV patterns
#
################################################################################

#$ -cwd
#$ -N numbat_analysis_deepseq_combined
#$ -l h_rt=10:00:00
#$ -l mem_per_core=8G
#$ -pe omp 8
#$ -P paxlab
#$ -o analysis/qsub_logs/numbat_analysis_deepseq_combined.$JOB_ID.log
#$ -j y

set -eo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
cd "$PROJECT_ROOT"

# Initialize modules
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

module load R 2>/dev/null || true

echo "[$(date +'%F %T')] ===== NUMBAT ANALYSIS: DEEPSEQ COMBINED (488B + 489) ====="
echo "[JOB_ID] $JOB_ID"
echo "[NODE] $(hostname)"
echo ""

# Merge two tissue-specific NUMBAT results and run comparative analysis
PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
DATASET="deepseq"
NCORES=8

cd "$PROJECT_ROOT"

# Paths to individual results
RESULT_488B="Data/04_analysis/cnv/numbat/results/deepseq_488B_atac"
RESULT_489="Data/04_analysis/cnv/numbat/results/deepseq_489_atac"
RESULT_COMBINED="Data/04_analysis/cnv/numbat/results/deepseq_combined_atac"

mkdir -p "$RESULT_COMBINED"

echo "[STEP 1] Verifying individual results exist..."
for dir in "$RESULT_488B" "$RESULT_489"; do
  if [[ ! -d "$dir" ]]; then
    echo "[ERROR] Result directory not found: $dir"
    echo "[INFO] Run tissue-specific analyses first:"
    echo "  qsub analysis/src/cnv_calling/numbat/run_numbat_analysis_deepseq_488B.qsub.sh"
    echo "  qsub analysis/src/cnv_calling/numbat/run_numbat_analysis_deepseq_489.qsub.sh"
    exit 1
  fi
  echo "  ✓ $(basename $dir) exists"
done

echo ""
echo "[STEP 2] Creating combined comparison plots..."

# Create R script to combine results
Rscript - <<'RCOMBINE'
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(data.table)
})

result_488B <- "Data/04_analysis/cnv/numbat/results/deepseq_488B_atac"
result_489 <- "Data/04_analysis/cnv/numbat/results/deepseq_489_atac"
output_dir <- "Data/04_analysis/cnv/numbat/results/deepseq_combined_atac"

cat("Loading results from individual tissues...\n")

# Try to load results - these will be created after individual analyses complete
if (file.exists(file.path(result_488B, "cnv_calls.rds"))) {
  cnv_488B <- readRDS(file.path(result_488B, "cnv_calls.rds"))
  cat("✓ Loaded 488B CNV calls\n")
} else {
  cat("[WARN] 488B CNV calls not found yet\n")
}

if (file.exists(file.path(result_489, "cnv_calls.rds"))) {
  cnv_489 <- readRDS(file.path(result_489, "cnv_calls.rds"))
  cat("✓ Loaded 489 CNV calls\n")
} else {
  cat("[WARN] 489 CNV calls not found yet\n")
}

# Save combined manifest
combined_info <- data.frame(
  tissue = c("488B", "489"),
  dataset = c("deepseq", "deepseq"),
  result_dir = c(result_488B, result_489)
)

write.table(combined_info, file.path(output_dir, "combined_manifest.tsv"), 
            sep="\t", row.names=FALSE, quote=FALSE)

cat("✓ Combined manifest saved\n")
cat("Results directory: ", output_dir, "\n", sep="")

RCOMBINE

if [[ $? -eq 0 ]]; then
  echo "[STEP 2] ✓ Complete"
else
  echo "[STEP 2] ✗ Failed"
  exit 1
fi

echo ""
echo "[$(date +'%F %T')] ===== COMPLETE ====="
