#!/bin/bash
# ============================================================================
# submit_all_tissues.sh
#
# Submit distribution modeling + UMAP jobs for all tissues
# (lowseq_489 already done, submitting 488B/deepseq variants)
#
# ============================================================================

set -u

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
SCRIPT_DIR="$PROJECT_ROOT/analysis/src/build_tissue/binsize_test_archr"
LOG_DIR="$PROJECT_ROOT/analysis/qsub_logs"

mkdir -p "$LOG_DIR"

log_msg() {
  echo "[$(date '+%F %T')] $1"
}

log_msg "===== Submitting expanded tissue analyses ====="

# Tissues to process (lowseq_489 already done, so skip it)
TISSUES=("lowseq_488B" "deepseq_488B" "deepseq_489")
TILESIZES=(500 5000)
BINARIZE_OPTIONS=("FALSE" "TRUE")

dist_job_count=0
umap_job_count=0

for tissue in "${TISSUES[@]}"; do
  log_msg ""
  log_msg "===== Processing: $tissue ====="

  # Submit distribution modeling jobs
  log_msg "Submitting distribution modeling jobs..."
  for tilesize in "${TILESIZES[@]}"; do
    JOB_NAME="distribution_modeling_${tissue}_${tilesize}bp"
    LOG_FILE="$LOG_DIR/${JOB_NAME}.log"
    ERR_FILE="$LOG_DIR/${JOB_NAME}.err"

    qsub -N "$JOB_NAME" \
      -pe omp 4 \
      -l mem_per_core=8G \
      -l h_rt=02:00:00 \
      -P paxlab \
      -o "$LOG_FILE" \
      -e "$ERR_FILE" \
      -cwd \
      -b y \
      "Rscript $SCRIPT_DIR/8_model_spot_count_distributions.R $tissue $tilesize" 2>&1 | tail -1

    ((dist_job_count++))
    log_msg "  Submitted: $JOB_NAME"
  done

  # Submit UMAP/ArchR variant jobs
  log_msg "Submitting ArchR variant (UMAP) jobs..."
  for tilesize in "${TILESIZES[@]}"; do
    for binarize in "${BINARIZE_OPTIONS[@]}"; do
      JOB_NAME="build_archr_${tissue}_${tilesize}bp_bin${binarize}"
      LOG_FILE="$LOG_DIR/${JOB_NAME}.log"
      ERR_FILE="$LOG_DIR/${JOB_NAME}.err"

      qsub -N "$JOB_NAME" \
        -pe omp 8 \
        -l mem_per_core=10G \
        -l h_rt=03:00:00 \
        -P paxlab \
        -o "$LOG_FILE" \
        -e "$ERR_FILE" \
        -cwd \
        -b y \
        "Rscript $SCRIPT_DIR/3_build_archr_variant_project.R $tissue $tilesize $binarize" 2>&1 | tail -1

      ((umap_job_count++))
      log_msg "  Submitted: $JOB_NAME"
    done
  done
done

log_msg ""
log_msg "===== Job submission complete ====="
log_msg "Submitted: $dist_job_count distribution modeling jobs + $umap_job_count UMAP jobs"
log_msg ""
log_msg "Monitor progress with:"
log_msg "  qstat -u preshita"
log_msg "  tail -f $LOG_DIR/distribution_modeling_*.log"
log_msg "  tail -f $LOG_DIR/build_archr_*.log"
log_msg ""
log_msg "Output directories:"
log_msg "  Distribution: $PROJECT_ROOT/analysis/binsize_comparison/distribution_modeling/"
log_msg "  UMAP/ArchR:   $PROJECT_ROOT/analysis/binsize_comparison/archr_projects/"
