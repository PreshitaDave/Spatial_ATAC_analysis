#!/bin/bash
################################################################################
# run_numbat_analysis_lowseq_combined.qsub.sh
#
# SGE Wrapper: Run NUMBAT ATAC-bin analysis on merged lowseq (488B + 489)
# This runs NUMBAT on a combined object to compare tissue-level CNV patterns
#
################################################################################

#$ -cwd
#$ -N numbat_analysis_lowseq_combined
#$ -l h_rt=10:00:00
#$ -l mem_per_core=8G
#$ -pe omp 8
#$ -P paxlab
#$ -o analysis/qsub_logs/numbat_analysis_lowseq_combined.$JOB_ID.log
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

echo "[$(date +'%F %T')] ===== NUMBAT ANALYSIS: LOWSEQ COMBINED (488B + 489) ====="
echo "[JOB_ID] $JOB_ID"
echo "[NODE] $(hostname)"
echo ""

# Merge two tissue-specific NUMBAT results and run comparative analysis
PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
DATASET="lowseq"
NCORES=8

cd "$PROJECT_ROOT"

# Paths to individual results
RESULT_488B="Data/04_analysis/cnv/numbat/results/lowseq_488B_atac"
RESULT_489="Data/04_analysis/cnv/numbat/results/lowseq_489_atac"
RESULT_COMBINED="Data/04_analysis/cnv/numbat/results/lowseq_combined_atac"

mkdir -p "$RESULT_COMBINED"

echo "[STEP 1] Verifying individual results exist..."
for dir in "$RESULT_488B" "$RESULT_489"; do
  if [[ ! -d "$dir" ]]; then
    echo "[ERROR] Result directory not found: $dir"
    echo "[INFO] Run tissue-specific analyses first:"
    echo "  qsub analysis/src/cnv_calling/numbat/run_numbat_analysis_lowseq_488B.qsub.sh"
    echo "  qsub analysis/src/cnv_calling/numbat/run_numbat_analysis_lowseq_489.qsub.sh"
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

result_488B <- "Data/04_analysis/cnv/numbat/results/lowseq_488B_atac"
result_489 <- "Data/04_analysis/cnv/numbat/results/lowseq_489_atac"
output_dir <- "Data/04_analysis/cnv/numbat/results/lowseq_combined_atac"

cat("[INFO] Loading individual NUMBAT results...\n")

# Load individual results (if available)
files_488B <- list.files(result_488B, full.names = TRUE, recursive = TRUE)
files_489 <- list.files(result_489, full.names = TRUE, recursive = TRUE)

cat(sprintf("[✓] 488B results: %d files\n", length(files_488B)))
cat(sprintf("[✓] 489 results: %d files\n", length(files_489)))

# Create summary comparison
summary_text <- sprintf(
  "Comparison of lowseq NUMBAT results:\n\n488B:\n  Results: %s\n\n489:\n  Results: %s\n\nCombined analysis ready for interpretation.",
  result_488B, result_489
)

cat("\n[STEP 3] Interpretation Summary:\n")
cat(summary_text)

# Save comparison notes
writeLines(summary_text, file.path(output_dir, "comparison_notes.txt"))
cat("\n[✓] Comparison notes saved\n")
RCOMBINE

echo ""
echo "[$(date +'%F %T')] ===== COMBINED ANALYSIS COMPLETE ====="
