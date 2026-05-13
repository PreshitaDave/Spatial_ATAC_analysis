#!/bin/bash -l
#$ -P paxlab
#$ -N tissue_vars
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l h_rt=24:00:00
#$ -l mem_per_core=8G
#$ -j n
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/tissue_vars.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/tissue_vars.$JOB_ID.err

set -euo pipefail

echo "[$(date '+%F %T')] Job $JOB_ID host=$(hostname) NSLOTS=$NSLOTS"

if ! command -v module >/dev/null 2>&1; then
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi
module load R

export DATASETS="${DATASETS:-deepseq,lowseq}"
export CHR_START="${CHR_START:-1}"
export CHR_END="${CHR_END:-22}"
export N_WORKERS="${N_WORKERS:-$NSLOTS}"
export OUT_DIR="${OUT_DIR:-/projectnb/paxlab/presh/projects/spatial_atac/Data/05_results/variant_calling/tissue_variants/tables}"

Rscript analysis/src/pipeline/variant_calling/11_comparing_tissue_variants.R
