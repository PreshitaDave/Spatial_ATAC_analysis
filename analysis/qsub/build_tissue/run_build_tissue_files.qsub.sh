#!/bin/bash -l
#$ -P paxlab
#$ -N tissue_files
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l h_rt=24:00:00
#$ -l mem_per_core=8G
#$ -j n
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/tissue_files.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/tissue_files.$JOB_ID.err

set -euo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"

# samtools 1.12 ships with bgzip and tabix in the same bin dir
export PATH="/share/pkg.7/samtools/1.12/install/bin:${PATH}"

echo "[$(date '+%F %T')] Job ${JOB_ID} starting on $(hostname)"
echo "[$(date '+%F %T')] NSLOTS=${NSLOTS:-8}"
echo "[$(date '+%F %T')] PATH=${PATH}"
echo "[$(date '+%F %T')] samtools: $(command -v samtools 2>/dev/null || echo MISSING)"
echo "[$(date '+%F %T')] bgzip:    $(command -v bgzip 2>/dev/null || echo MISSING)"
echo "[$(date '+%F %T')] tabix:    $(command -v tabix 2>/dev/null || echo MISSING)"

export THREADS="${NSLOTS:-8}"
export PROJECT_ROOT
export FORCE="${FORCE:-0}"

cd "${PROJECT_ROOT}"
bash analysis/src/build_tissue/build_tissue_files.sh

echo "[$(date '+%F %T')] Job ${JOB_ID} complete"
