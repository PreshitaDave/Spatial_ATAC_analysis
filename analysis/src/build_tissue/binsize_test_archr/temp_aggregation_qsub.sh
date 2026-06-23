#!/bin/bash
#$ -N aggregate_metrics
#$ -l h_rt=00:30:00
#$ -pe omp 1
#$ -P paxlab
#$ -l mem_per_core=4G
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/aggregate_metrics.log
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/aggregate_metrics.err
#$ -hold_jid 6143563,6143564,6143565,6143566,6143567,6143568,6143569,6143570,6143571,6143572,6143573,6143574,6143575,6143576,6143577,6143578

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

echo "[2026-06-22 11:48:11] [start] Running sparsity comparison aggregation"

cd "${PROJECT_ROOT}"

Rscript "${SCRIPT_DIR}/compare_arrow_sparsity.R"

exit_code=$?
if [[ ${exit_code} -eq 0 ]]; then
  echo "[2026-06-22 11:48:11] [done] Sparsity comparison completed"
else
  echo "[2026-06-22 11:48:11] [error] Comparison failed with exit code: ${exit_code}"
  exit ${exit_code}
fi
