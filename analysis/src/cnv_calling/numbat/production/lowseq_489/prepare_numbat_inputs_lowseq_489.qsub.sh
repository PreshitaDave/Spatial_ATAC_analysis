#!/bin/bash
#$ -cwd
#$ -N numbat_prep_lowseq_489
#$ -l h_rt=08:00:00
#$ -l mem_per_core=8G
#$ -pe omp 8
#$ -P paxlab
#$ -o analysis/qsub_logs/numbat_prep_lowseq_489.$JOB_ID.out
#$ -e analysis/qsub_logs/numbat_prep_lowseq_489.$JOB_ID.err

################################################################################
# SGE Wrapper: Prepare NUMBAT inputs for lowseq tissue 489
# 
# Generates:
#   1. Pileup/Phase: allele counts
#   2. ATAC Bins: 220kb bin × cell matrix
#   3. Reference: aggregated from all cells
#
# Runtime: ~3 hours (smaller tissue)
# Resources: 8 cores, 64GB memory
################################################################################

set -eo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
cd "$PROJECT_ROOT"

echo "==============================================="
echo "NUMBAT Input Preparation: lowseq 489"
echo "==============================================="
echo "Job ID: $JOB_ID"
echo "Host: $(hostname)"
echo "Start Time: $(date)"
echo ""

# Initialize module system
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

# Verify required tools
echo "[STEP 1] Verifying environment..."
module load R 2>/dev/null || true
which Rscript && echo "  ✓ R loaded" || { echo "  ✗ R not found"; exit 1; }

# Add calicost_env to PATH for cellsnp-lite and vcftools
export PATH="/projectnb/paxlab/presh/env/calicost_env/bin:$PATH"
which cellsnp-lite && echo "  ✓ cellsnp-lite loaded" || { echo "  ✗ cellsnp-lite not found"; exit 1; }

# Run main preparation script
echo "[STEP 2] Running input preparation (lowseq 489)..."
bash analysis/src/cnv_calling/numbat/numbat/prepare_numbat_inputs.sh lowseq 489

echo ""
echo "==============================================="
echo "Completed: $(date)"
echo "==============================================="
