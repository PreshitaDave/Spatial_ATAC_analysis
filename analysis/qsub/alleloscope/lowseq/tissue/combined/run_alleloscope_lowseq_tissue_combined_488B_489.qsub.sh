#!/bin/bash -l
#$ -P paxlab
#$ -N allelo_low_comb
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l h_rt=48:00:00
#$ -l mem_per_core=6G
#$ -j n
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/allelo_low_comb.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/allelo_low_comb.$JOB_ID.err

set -euo pipefail

if ! command -v module >/dev/null 2>&1; then
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi

module load R

echo "[$(date '+%F %T')] Job ${JOB_ID} starting on $(hostname) NSLOTS=${NSLOTS:-8}"

cd /projectnb/paxlab/presh/projects/spatial_atac
Rscript --no-save --no-restore \
  analysis/src/cnv_calling/alleloscope/alleloscope/lowseq/tissue/combined/run_alleloscope_lowseq_tissue_combined_488B_489.R

echo "[$(date '+%F %T')] Job ${JOB_ID} complete"
