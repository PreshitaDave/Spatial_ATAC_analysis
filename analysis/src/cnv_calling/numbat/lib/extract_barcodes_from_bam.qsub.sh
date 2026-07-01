#!/bin/bash
#$ -P paxlab
#$ -N extract_bc16_deepseq
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 4
#$ -l mem_per_core=8G
#$ -l h_rt=8:00:00
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/extract_bc16_deepseq.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/extract_bc16_deepseq.$JOB_ID.err
#$ -j n

set -eo pipefail

# Initialize module system
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] ===== EXTRACT BC16 BARCODES FROM DEEPSEQ MERGED BAM ====="
echo "[COMPUTE NODE] $(hostname)"
echo ""

# Run extraction script
export PROJECT_ROOT
bash "${PROJECT_ROOT}/analysis/src/cnv_calling/numbat/extract_barcodes_from_bam.sh"

echo ""
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ===== EXTRACTION COMPLETE ====="
