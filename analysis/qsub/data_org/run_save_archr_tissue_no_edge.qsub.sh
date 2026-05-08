#!/bin/bash -l
#$ -P paxlab
#$ -N archr_noedge
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l h_rt=24:00:00
#$ -l mem_per_core=8G
#$ -j n
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/archr_noedge.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/archr_noedge.$JOB_ID.err

set -euo pipefail

if ! command -v module >/dev/null 2>&1; then
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi

module load R

echo "[$(date '+%F %T')] [start] Job ${JOB_ID} host=$(hostname) NSLOTS=${NSLOTS:-8}"

cd /projectnb/paxlab/presh/projects/spatial_atac
Rscript analysis/src/data_org/save_archr_tissue_no_edge.R

echo "[$(date '+%F %T')] [done] Job ${JOB_ID} complete"
