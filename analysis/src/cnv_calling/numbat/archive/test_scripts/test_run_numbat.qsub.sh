#!/bin/bash
#$ -N test_run_numbat
#$ -l h_rt=01:00:00
#$ -pe omp 4
#$ -P paxlab
#$ -l mem_per_core=8G
#$ -o analysis/qsub_logs/test_run_numbat_$JOB_ID.log

set +eo pipefail  # Capture all output even if there are errors

# Module initialization
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

module load R
cd /projectnb/paxlab/presh/projects/spatial_atac

echo "[WRAPPER] Starting run_numbat() test"
echo "[WRAPPER] Job: $JOB_ID"
echo "[WRAPPER] Time: $(date)"
echo ""

# Use timeout to kill if it hangs longer than 50 minutes
timeout 3000s Rscript analysis/src/cnv_calling/numbat/test_run_numbat_call.R \
  "lowseq_489" \
  "Data/04_analysis/cnv/numbat/inputs/lowseq_489/atac_bin/lowseq_489_atac_bin.rds" \
  "Data/04_analysis/cnv/numbat/inputs/lowseq_489/alleles/lowseq_489_atac_allele_counts.tsv.gz" \
  "4" \
  "Data/04_analysis/cnv/numbat/results/lowseq_489/test_run_numbat/"

EXIT_CODE=$?

echo ""
echo "[WRAPPER] Script exited with code: $EXIT_CODE"
if [[ $EXIT_CODE -eq 124 ]]; then
  echo "[WRAPPER] ✗ TIMEOUT: Job ran longer than 50 minutes"
elif [[ $EXIT_CODE -eq 0 ]]; then
  echo "[WRAPPER] ✓ COMPLETED: Job finished successfully"
else
  echo "[WRAPPER] ✗ FAILED: Job exited with error code $EXIT_CODE"
fi
echo "[WRAPPER] Time: $(date)"
