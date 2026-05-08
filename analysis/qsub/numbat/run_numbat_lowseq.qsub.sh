#!/bin/bash -l
#$ -P paxlab
#$ -N numbat_low
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=8G
#$ -l h_rt=48:00:00
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_low.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_low.$JOB_ID.err
#$ -j n

set -euo pipefail

module load R
module load samtools

# cellsnp-lite lives in the calicost conda env; eagle binary is in external/
export PATH="/projectnb/paxlab/presh/env/calicost_env/bin:/projectnb/paxlab/presh/software/external/Eagle_v2.4.1:${PATH}"

cd /projectnb/paxlab/presh/projects/spatial_atac

echo "[$(date '+%F %T')] START NUMBAT lowseq prep + ATAC-bin run"
export DATASET=lowseq
export NCORES="${NSLOTS:-8}"
export NUMBAT_REPO="/projectnb/paxlab/presh/projects/spatial_atac/Data/numbat/numbat_repo"
export CHROMS="${CHROMS:-chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22}"
# Panel symlink dir with chr{N}.genotypes.bcf naming that Eagle/pileup_and_phase.R expects
export PHASE_PANEL="/projectnb/paxlab/presh/projects/spatial_atac/Data/numbat/reference/phased_panel_bcf_links"

bash analysis/src/numbat/prepare_numbat_atac_inputs.sh
bash analysis/src/numbat/run_numbat_atac_bin.sh

echo "[$(date '+%F %T')] DONE NUMBAT lowseq"
