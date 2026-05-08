#!/bin/bash -l
#$ -P paxlab
#$ -N allelo_deepseq
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/allelo_deepseq.$JOB_ID.log
#$ -pe omp 8
#$ -l mem_per_core=8G
#$ -l h_rt=48:00:00

set -euo pipefail

# Ensure environment modules are initialized on all SCC nodes/shell modes.
if ! command -v module >/dev/null 2>&1; then
  # shellcheck source=/dev/null
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi

echo "[$(date)] Job ${JOB_ID} starting on $(hostname)"
echo "  NSLOTS=${NSLOTS}"

module load R

cd /projectnb/paxlab/presh/projects/spatial_atac

Rscript --no-save --no-restore \
  analysis/src/alleloscope/run_alleloscope_deepseq.R

echo "[$(date)] Job ${JOB_ID} complete"
