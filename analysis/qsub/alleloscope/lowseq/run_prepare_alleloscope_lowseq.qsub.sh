#!/bin/bash -l
#$ -P paxlab
#$ -N allelo_low
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=8G
#$ -l h_rt=48:00:00
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/allelo_low.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/allelo_low.$JOB_ID.err
#$ -j n

set -euo pipefail

if ! command -v module >/dev/null 2>&1; then
	[[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi

module load R

cd /projectnb/paxlab/presh/projects/spatial_atac
echo "[$(date '+%F %T')] Starting lowseq Alleloscope prep (JOB_ID=${JOB_ID}, NSLOTS=${NSLOTS:-8})"
export OMP_NUM_THREADS="${NSLOTS:-8}"
R --vanilla -q -f analysis/src/alleloscope/lowseq/prepare_alleloscope_lowseq.R
echo "[$(date '+%F %T')] Finished lowseq Alleloscope prep"
