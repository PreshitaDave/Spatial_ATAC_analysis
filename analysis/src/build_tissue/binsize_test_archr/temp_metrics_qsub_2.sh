#!/bin/bash
#$ -N recompute_metrics
#$ -l h_rt=01:00:00
#$ -pe omp 4
#$ -P paxlab
#$ -l mem_per_core=8G
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/recompute_metrics_deepseq_488B_500_TRUE.log
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/recompute_metrics_deepseq_488B_500_TRUE.err

set -euo pipefail

# Initialize module system
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

module load R

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
SCRIPT_DIR="${PROJECT_ROOT}/analysis/src/build_tissue/binsize_test_archr"

echo "[2026-06-22 11:48:06] [start] Recomputing metrics: deepseq_488B 500bp binarize=TRUE"

cd "${PROJECT_ROOT}"

Rscript "${SCRIPT_DIR}/create_arrow_variants.R" "deepseq_488B" "500" "TRUE"

exit_code=$?
if [[ ${exit_code} -eq 0 ]]; then
  echo "[2026-06-22 11:48:06] [done] Successfully recomputed metrics for deepseq_488B"
else
  echo "[2026-06-22 11:48:06] [error] Failed with exit code: ${exit_code}"
  exit ${exit_code}
fi
