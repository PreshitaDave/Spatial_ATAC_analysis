#!/bin/bash -l
#$ -P paxlab
#$ -N tissue_total
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l h_rt=24:00:00
#$ -l mem_per_core=8G
#$ -j n
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/tissue_total.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/tissue_total.$JOB_ID.err

set -euo pipefail

echo "[$(date '+%F %T')] [start] Job ${JOB_ID} host=$(hostname) NSLOTS=${NSLOTS:-8}"

cd /projectnb/paxlab/presh/projects/spatial_atac

export PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
export THREADS="${NSLOTS:-8}"
export FORCE="0"
export BUILD_BAMS="1"
export BUILD_FRAGMENTS="1"

bash analysis/src/alleloscope/build_tissue_files.sh

echo "[$(date '+%F %T')] [done] Job ${JOB_ID} complete"
