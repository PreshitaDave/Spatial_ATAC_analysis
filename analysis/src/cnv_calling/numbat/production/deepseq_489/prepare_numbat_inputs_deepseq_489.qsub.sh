#!/bin/bash
#$ -cwd
#$ -N numbat_prep_deepseq_489
#$ -l h_rt=08:00:00
#$ -l mem_per_core=8G
#$ -pe omp 8
#$ -P paxlab
#$ -o analysis/qsub_logs/numbat_prep_deepseq_489.$JOB_ID.out
#$ -e analysis/qsub_logs/numbat_prep_deepseq_489.$JOB_ID.err

################################################################################
# SGE Wrapper: Prepare NUMBAT inputs for deepseq tissue 489
################################################################################

set -eo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
cd "$PROJECT_ROOT"

echo "NUMBAT Input Preparation: deepseq 489 (Job $JOB_ID @ $(hostname))"

# Initialize module system
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  [[ -f "$profile_file" ]] && { . "$profile_file" 2>/dev/null || true; break; }
done
set -u

module load R 2>/dev/null || true

# Add calicost_env to PATH for cellsnp-lite and vcftools
export PATH="/projectnb/paxlab/presh/env/calicost_env/bin:$PATH"
which cellsnp-lite >/dev/null 2>&1 && which Rscript && bash analysis/src/cnv_calling/numbat/prepare_numbat_inputs.sh deepseq 489 || { echo "ERROR: R or cellsnp-lite not found"; exit 1; }
