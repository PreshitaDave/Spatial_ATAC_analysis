#!/bin/bash -l
#$ -N rebuild_frag_489
#$ -P paxlab
#$ -l h_rt=04:00:00
#$ -l mem_per_core=16G
#$ -pe omp 4
#$ -cwd
#$ -j y
#$ -o analysis/qsub_logs/rebuild_fragments_489.$JOB_ID.out
#$ -e analysis/qsub_logs/rebuild_fragments_489.$JOB_ID.err

set -euo pipefail

echo "[$(date '+%F %T')] Job $JOB_ID host=$(hostname) NSLOTS=$NSLOTS" \
  >> analysis/qsub_logs/rebuild_fragments_489.$JOB_ID.out

if ! command -v module >/dev/null 2>&1; then
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi
module load R

Rscript analysis/src/data_org/rebuild_fragments_tissue_489.R
EXIT=$?
echo "R exit code: $EXIT"
exit $EXIT
