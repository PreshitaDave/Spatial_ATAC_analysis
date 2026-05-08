#!/bin/bash -l
#$ -P paxlab
#$ -N numbat_low
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=8G
#$ -l h_rt=48:00:00
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_low.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_low.$JOB_ID.err
#$ -j n

set -euo pipefail

module load R

cd /projectnb/paxlab/presh/projects/spatial_atac

echo "[$(date '+%F %T')] START NUMBAT lowseq analysis-only (inputs already prepared)"
export DATASET=lowseq
export NCORES="${NSLOTS:-8}"
export NUMBAT_REPO="/projectnb/paxlab/presh/projects/spatial_atac/Data/numbat/numbat_repo"

# Skip prepare_numbat_atac_inputs.sh — inputs already exist
bash analysis/src/numbat/run_numbat_atac_bin.sh

echo "[$(date '+%F %T')] DONE NUMBAT lowseq"
