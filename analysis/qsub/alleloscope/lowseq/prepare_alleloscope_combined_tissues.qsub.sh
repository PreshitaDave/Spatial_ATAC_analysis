#!/bin/bash
#$ -P paxlab
#$ -N allelo_low_prep_combined
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l h_rt=24:00:00
#$ -l mem_per_core=6G
#$ -j n
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/allelo_low_prep_combined.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/allelo_low_prep_combined.$JOB_ID.err

set -euo pipefail

# CRITICAL: Initialize module system before using 'module' command
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

module load R

echo "[$(date '+%F %T')] Job ${JOB_ID} starting on $(hostname) NSLOTS=${NSLOTS:-4}"

cd /projectnb/paxlab/presh/projects/spatial_atac

echo "[$(date '+%F %T')] Running combined tissue prep..."
Rscript --no-save --no-restore \
  analysis/src/alleloscope/alleloscope/lowseq/prepare_alleloscope_combined_tissues.R

echo "[$(date '+%F %T')] Combined prep complete, verifying outputs..."
if [[ -f Data/04_analysis/cnv/alleloscope/lowseq_tissue_from_existing/combined_488B_489/alt_all.mtx ]] && \
   [[ -f Data/04_analysis/cnv/alleloscope/lowseq_tissue_from_existing/combined_488B_489/ref_all.mtx ]] && \
   [[ -f Data/04_analysis/cnv/alleloscope/lowseq_tissue_from_existing/combined_488B_489/barcodes.tsv ]]; then
  echo "[$(date '+%F %T')] ✓ All combined inputs verified"
else
  echo "[$(date '+%F %T')] ✗ ERROR: Missing combined input files"
  exit 1
fi

echo "[$(date '+%F %T')] Job ${JOB_ID} complete"
