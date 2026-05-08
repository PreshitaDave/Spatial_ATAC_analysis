#!/bin/bash -l
#$ -P paxlab
#$ -N lowseq_tfix
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 4
#$ -l h_rt=12:00:00
#$ -l mem_per_core=16G
#$ -j n
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/lowseq_tfix.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/lowseq_tfix.$JOB_ID.err

set -euo pipefail

if ! command -v module >/dev/null 2>&1; then
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi

module load R

echo "[$(date '+%F %T')] [start] Job ${JOB_ID} host=$(hostname) NSLOTS=${NSLOTS:-4}"

cd /projectnb/paxlab/presh/projects/spatial_atac
Rscript analysis/src/data_org/prepare_lowseq_alleloscope_tissue_from_existing.R

echo "[$(date '+%F %T')] [done] Job ${JOB_ID} complete"
