#!/bin/bash
# ============================================================================
# submit_archr_variant_pilot.sh
#
# Submit parallel qsub jobs to build ArchR projects and generate UMAP/gene-score
# plots for all variants of a pilot tissue (lowseq_489).
#
# Then submit an aggregation job that depends on all 4 pilot jobs,
# which generates comparison plots across the variants.
#
# Usage:
#   bash submit_archr_variant_pilot.sh
#
# ============================================================================

set -e

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
SCRIPT_DIR="$PROJECT_ROOT/analysis/src/build_tissue/binsize_test_archr"
LOG_DIR="$PROJECT_ROOT/analysis/qsub_logs"

mkdir -p "$LOG_DIR"

log_msg() {
  echo "[$(date '+%F %T')] $1"
}

log_msg "===== Submitting ArchR variant pilot jobs (lowseq_489) ====="

# Pilot tissue(s) - change this array to expand to other tissues
TISSUES=("lowseq_489")
TILESIZES=(500 5000)
BINARIZE_OPTIONS=("FALSE" "TRUE")

# Array to hold job IDs
JOB_IDS=()

# Submit per-variant jobs
log_msg "Submitting variant builder jobs..."

for tissue in "${TISSUES[@]}"; do
  for tilesize in "${TILESIZES[@]}"; do
    for binarize in "${BINARIZE_OPTIONS[@]}"; do

      JOB_NAME="archr_variant_${tissue}_${tilesize}bp_binarize${binarize}"
      LOG_FILE="$LOG_DIR/archr_variant_${tissue}_${tilesize}bp_binarize${binarize}.log"
      ERR_FILE="$LOG_DIR/archr_variant_${tissue}_${tilesize}bp_binarize${binarize}.err"

      log_msg "Submitting: $JOB_NAME"

      JOB_ID=$(qsub -N "$JOB_NAME" \
        -pe omp 8 \
        -l mem_per_core=8G \
        -l h_rt=04:00:00 \
        -P paxlab \
        -o "$LOG_FILE" \
        -e "$ERR_FILE" \
        -cwd \
        -b y \
        "Rscript $SCRIPT_DIR/3_build_archr_variant_project.R $tissue $tilesize $binarize" 2>&1 | awk '{print $3}')

      JOB_IDS+=("$JOB_ID")
      log_msg "  Submitted with job ID: $JOB_ID"

    done
  done
done

log_msg "Submitted ${#JOB_IDS[@]} variant builder jobs"

# Build hold string for dependency
HOLD_STR=$(IFS=, ; echo "${JOB_IDS[*]}")

log_msg "Building hold string for aggregator job: $HOLD_STR"

# Submit aggregation job that depends on all variant jobs
AGG_JOB_NAME="archr_variant_aggregator_pilot"
AGG_LOG_FILE="$LOG_DIR/archr_variant_aggregator_pilot.log"
AGG_ERR_FILE="$LOG_DIR/archr_variant_aggregator_pilot.err"

log_msg "Submitting aggregation job: $AGG_JOB_NAME"

AGG_JOB_ID=$(qsub -N "$AGG_JOB_NAME" \
  -pe omp 4 \
  -l mem_per_core=4G \
  -l h_rt=00:30:00 \
  -P paxlab \
  -o "$AGG_LOG_FILE" \
  -e "$AGG_ERR_FILE" \
  -cwd \
  -hold_jid "$HOLD_STR" \
  -b y \
  "Rscript $SCRIPT_DIR/4_compare_archr_variant_umap.R lowseq_489" 2>&1 | awk '{print $3}')

log_msg "Submitted aggregation job with ID: $AGG_JOB_ID"

log_msg "===== Job submission complete ====="
log_msg "Variant builder jobs: ${JOB_IDS[@]}"
log_msg "Aggregator job: $AGG_JOB_ID (depends on all variant jobs)"
log_msg ""
log_msg "Monitor progress with:"
log_msg "  qstat"
log_msg "  tail -f $LOG_DIR/archr_variant_*.log"
log_msg ""
log_msg "Final outputs will be generated in:"
log_msg "  $PROJECT_ROOT/analysis/binsize_comparison/"
