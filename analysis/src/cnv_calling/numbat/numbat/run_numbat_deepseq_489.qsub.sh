#!/bin/bash -l
#$ -P paxlab
#$ -N numbat_deep_489
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=8G
#$ -l h_rt=48:00:00
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_deep_489.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_deep_489.$JOB_ID.err
#$ -j n

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
module load samtools

# cellsnp-lite lives in the calicost conda env; eagle binary is in external/
export PATH="/projectnb/paxlab/presh/env/calicost_env/bin:/projectnb/paxlab/presh/software/external/Eagle_v2.4.1:${PATH}"

cd /projectnb/paxlab/presh/projects/spatial_atac

echo "[$(date '+%F %T')] START NUMBAT deepseq_489 prep + ATAC-bin run"
export DATASET=deepseq
export TISSUE=489
export NCORES="${NSLOTS:-8}"
export NUMBAT_REPO="/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/cnv/numbat/numbat_repo"
export OUT_DIR="/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/cnv/numbat/results/deepseq_489/atac_only_run2"
export CHROMS="${CHROMS:-chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22}"

# Reference directory with chr{N}.genotypes.bcf naming that Eagle/pileup_and_phase.R expects
export PHASE_PANEL="/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/cnv/numbat/reference/phased_panel_bcf_links"

echo "[$(date '+%F %T')] Parameters:"
echo "  DATASET=${DATASET}"
echo "  TISSUE=${TISSUE}"
echo "  NCORES=${NCORES}"
echo "  OUT_DIR=${OUT_DIR}"
echo "  NUMBAT_REPO=${NUMBAT_REPO}"

# Run prepare and numbat analysis
echo "[$(date '+%F %T')] Running prepare_numbat_atac_inputs.sh..."
bash analysis/src/numbat/numbat/prepare_numbat_atac_inputs.sh

echo "[$(date '+%F %T')] Running run_numbat_atac_bin.sh..."
bash analysis/src/numbat/numbat/run_numbat_atac_bin.sh

echo "[$(date '+%F %T')] DONE NUMBAT deepseq_489 (atac_run2: 220kb bins)"
