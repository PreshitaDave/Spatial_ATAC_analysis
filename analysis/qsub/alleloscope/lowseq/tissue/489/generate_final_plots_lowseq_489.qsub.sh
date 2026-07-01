#!/bin/bash -l
#$ -P paxlab
#$ -N allelo_489_plots
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l h_rt=04:00:00
#$ -l mem_per_core=8G
#$ -j n
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/allelo_489_plots.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/allelo_489_plots.$JOB_ID.err

set -euo pipefail

if ! command -v module >/dev/null 2>&1; then
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi

module load R

echo "[$(date '+%F %T')] Job ${JOB_ID} starting on $(hostname) NSLOTS=${NSLOTS:-8}"

cd /projectnb/paxlab/presh/projects/spatial_atac

echo "[$(date '+%F %T')] ========== Generating Final Plots for Lowseq_489 =========="
echo "[$(date '+%F %T')] Input object: Data/alleloscope/lowseq_tissue_from_existing/489/output/rds/Obj_final.rds"
echo "[$(date '+%F %T')] Output: Data/alleloscope/lowseq_tissue_from_existing/489/output/plots/"

Rscript --no-save --no-restore \
  analysis/src/cnv_calling/alleloscope/alleloscope/lowseq/tissue/489/generate_final_plots_lowseq_489.R

EXITCODE=$?

echo "[$(date '+%F %T')] ========== Script completed with exit code: $EXITCODE =========="

if [[ $EXITCODE -eq 0 ]]; then
  echo "[$(date '+%F %T')] ✓ SUCCESS: Heatmap and plots generated"
  ls -lh /projectnb/paxlab/presh/projects/spatial_atac/Data/alleloscope/lowseq_tissue_from_existing/489/output/plots/step6*
else
  echo "[$(date '+%F %T')] ✗ FAILED: Check logs for details"
fi

echo "[$(date '+%F %T')] Job ${JOB_ID} complete"

exit $EXITCODE
