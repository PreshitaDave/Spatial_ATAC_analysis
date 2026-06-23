#!/bin/bash
#$ -N create_arrow
#$ -l h_rt=24:00:00
#$ -pe omp 8
#$ -P paxlab
#$ -l mem_per_core=8G
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/create_arrow_deepseq_489_500_FALSE.log
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/create_arrow_deepseq_489_500_FALSE.err

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

echo "[2026-06-21 23:19:59] [start] Creating arrow: deepseq_489 500bp binarize=FALSE"

cd "${PROJECT_ROOT}"

Rscript "${SCRIPT_DIR}/create_arrow_variants.R" "deepseq_489" "500" "FALSE"

exit_code=$?
if [[ ${exit_code} -eq 0 ]]; then
  echo "[2026-06-21 23:19:59] [done] Successfully created arrow for deepseq_489"
else
  echo "[2026-06-21 23:19:59] [error] Failed with exit code: ${exit_code}"
  exit ${exit_code}
fi
