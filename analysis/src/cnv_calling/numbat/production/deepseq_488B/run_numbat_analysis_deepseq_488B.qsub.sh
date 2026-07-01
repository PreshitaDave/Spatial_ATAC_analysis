#!/bin/bash
################################################################################
# run_numbat_analysis_deepseq_488B.qsub.sh
#
# SGE Wrapper: Run NUMBAT ATAC-bin analysis for deepseq 488B
#
################################################################################

#$ -cwd
#$ -N numbat_analysis_deepseq_488B
#$ -l h_rt=08:00:00
#$ -l mem_per_core=8G
#$ -pe omp 8
#$ -P paxlab
#$ -o analysis/qsub_logs/numbat_analysis_deepseq_488B.$JOB_ID.log
#$ -j y

set -eo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
cd "$PROJECT_ROOT"

# Initialize modules
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

echo "[$(date +'%F %T')] ===== NUMBAT ANALYSIS: DEEPSEQ 488B ====="
echo "[JOB_ID] $JOB_ID"
echo "[NODE] $(hostname)"
echo ""

bash analysis/src/cnv_calling/numbat/run_numbat_analysis_atac.sh deepseq 488B

echo ""
echo "[$(date +'%F %T')] ===== COMPLETE ====="
