#!/bin/bash
# ============================================================================
# submit_spatial_smoothing.sh
#
# Submit spatial smoothing jobs for all 4 tissues via qsub
# ============================================================================

set -u

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
SCRIPT_DIR="${PROJECT_ROOT}/analysis/src/build_tissue/binsize_test_archr"
SCRIPT="${SCRIPT_DIR}/add_spatial_smoothing_lsi.R"

# Array of tissues
TISSUES=("lowseq_489" "lowseq_488B" "deepseq_488B" "deepseq_489")

LOG_DIR="${PROJECT_ROOT}/analysis/qsub_logs"
mkdir -p "${LOG_DIR}"

echo "[$(date '+%F %T')] Submitting spatial smoothing jobs..."

for tissue in "${TISSUES[@]}"; do
  echo "[$(date '+%F %T')] Submitting: ${tissue}"

  JOB_NAME="spatial_smooth_${tissue}"
  LOG_FILE="${LOG_DIR}/${JOB_NAME}.log"
  ERR_FILE="${LOG_DIR}/${JOB_NAME}.err"

  qsub -N "${JOB_NAME}" \
    -pe omp 8 \
    -l mem_per_core=8G \
    -l h_rt=03:00:00 \
    -P paxlab \
    -o "${LOG_FILE}" \
    -e "${ERR_FILE}" \
    -cwd \
    -b y \
    "Rscript ${SCRIPT} ${tissue}" 2>&1 | tail -1

  sleep 1  # Small delay between submissions
done

echo "[$(date '+%F %T')] All jobs submitted."
echo "[$(date '+%F %T')] Monitor with: qstat -u preshita"
echo "[$(date '+%F %T')] Outputs will be in: ${PROJECT_ROOT}/analysis/binsize_comparison/spatial_smoothing/"
