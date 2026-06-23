#!/bin/bash
# ============================================================================
# submit_distribution_modeling.sh
#
# Submit qsub jobs to characterize per-tile count distributions for spatial ATAC-seq
# across different tile sizes.
#
# Jobs are lighter-weight than the UMAP pipeline (no LSI/clustering/UMAP),
# so they can run in parallel with faster turnaround.
#
# Usage:
#   bash submit_distribution_modeling.sh
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

log_msg "===== Submitting distribution modeling jobs ====="

# Pilot tissue (same as UMAP validation, for direct comparison)
TISSUE="lowseq_489"
TILESIZES=(500 5000)

log_msg "Submitting distribution modeling jobs for: $TISSUE"

for tilesize in "${TILESIZES[@]}"; do
  JOB_NAME="distribution_modeling_${TISSUE}_${tilesize}bp"
  LOG_FILE="$LOG_DIR/distribution_modeling_${TISSUE}_${tilesize}bp.log"
  ERR_FILE="$LOG_DIR/distribution_modeling_${TISSUE}_${tilesize}bp.err"

  log_msg "Submitting: $JOB_NAME"

  qsub -N "$JOB_NAME" \
    -pe omp 4 \
    -l mem_per_core=8G \
    -l h_rt=02:00:00 \
    -P paxlab \
    -o "$LOG_FILE" \
    -e "$ERR_FILE" \
    -cwd \
    -b y \
    "Rscript $SCRIPT_DIR/model_spot_count_distributions_simplified.R $TISSUE $tilesize" 2>&1 | tail -1

  log_msg "  Job submitted for $TISSUE $tilesize bp"
done

log_msg "===== Job submission complete ====="
log_msg ""
log_msg "Monitor progress with:"
log_msg "  qstat -u preshita"
log_msg "  tail -f $LOG_DIR/distribution_modeling_*.log"
log_msg ""
log_msg "Output directory:"
log_msg "  $PROJECT_ROOT/analysis/binsize_comparison/distribution_modeling/"
