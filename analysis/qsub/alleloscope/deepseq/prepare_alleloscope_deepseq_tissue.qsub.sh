#!/bin/bash -l
#$ -P paxlab
#$ -N allelo_deep_tissue_prep
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/allelo_deep_tissue_prep.$JOB_ID.log
#$ -pe omp 8
#$ -l mem_per_core=8G
#$ -l h_rt=24:00:00

set -euo pipefail

if ! command -v module >/dev/null 2>&1; then
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi

echo "[$(date)] Job ${JOB_ID} starting on $(hostname)"
echo "  NSLOTS=${NSLOTS}"

module load R

cd /projectnb/paxlab/presh/projects/spatial_atac

Rscript --no-save --no-restore \
  analysis/src/alleloscope/alleloscope/deepseq/prepare_alleloscope_by_tissue.R deepseq

echo "[$(date)] Job ${JOB_ID} complete"
