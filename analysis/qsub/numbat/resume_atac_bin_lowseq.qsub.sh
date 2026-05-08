#!/bin/bash -l
#$ -P paxlab
#$ -N numbat_atac_low
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=8G
#$ -l h_rt=24:00:00
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_atac_low.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_atac_low.$JOB_ID.err
#$ -j n

set -euo pipefail

module load R

export DATASET=lowseq
export NCORES="${NSLOTS:-8}"
export NUMBAT_REPO="/projectnb/paxlab/presh/projects/spatial_atac/Data/numbat/numbat_repo"

echo "[$(date '+%F %T')] START NUMBAT ATAC-bin analysis (dataset=${DATASET})"

bash /projectnb/paxlab/presh/projects/spatial_atac/analysis/src/numbat/run_numbat_atac_bin.sh

echo "[$(date '+%F %T')] DONE"
