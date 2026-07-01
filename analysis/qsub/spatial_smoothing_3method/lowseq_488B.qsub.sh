#!/bin/bash -l
#$ -P paxlab
#$ -N smooth_3method_lowseq_488B
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l h_rt=48:00:00
#$ -l mem_per_core=4G
#$ -j n
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/smooth_3method_lowseq_488B.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/smooth_3method_lowseq_488B.$JOB_ID.err

set -euo pipefail

if ! command -v module >/dev/null 2>&1; then
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi

module load R

echo "[$(date '+%F %T')] Job ${JOB_ID} starting: 3-method spatial smoothing comparison (lowseq_488B)"
Rscript analysis/src/build_tissue/binsize_test_archr/compare_spatial_smoothing_methods.R lowseq_488B
echo "[$(date '+%F %T')] Job ${JOB_ID} completed"
