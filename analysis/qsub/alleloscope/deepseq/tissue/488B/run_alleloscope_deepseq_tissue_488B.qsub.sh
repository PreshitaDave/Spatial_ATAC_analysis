#!/bin/bash -l
#$ -P paxlab
#$ -N allelo_deep_488B
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l h_rt=48:00:00
#$ -l mem_per_core=6G
#$ -j n
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/allelo_deep_488B.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/allelo_deep_488B.$JOB_ID.err

set -euo pipefail

if ! command -v module >/dev/null 2>&1; then
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi

module load R

echo "[$(date '+%F %T')] Job ${JOB_ID} starting on $(hostname) NSLOTS=${NSLOTS:-8}"

cd /projectnb/paxlab/presh/projects/spatial_atac
Rscript --no-save --no-restore \
  analysis/src/alleloscope/alleloscope/deepseq/tissue/488B/run_alleloscope_deepseq_tissue_488B.R

echo "[$(date '+%F %T')] Job ${JOB_ID} complete"
