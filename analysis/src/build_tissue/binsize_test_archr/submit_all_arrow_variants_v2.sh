#!/bin/bash
# ============================================================================
# submit_all_arrow_variants_v2.sh
#
# Submit 16 parallel qsub jobs for all tissue x tilesize x binarize combinations
# Each job will create proper .log and .err files
#
# ============================================================================

set -euo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
SCRIPT_DIR="${PROJECT_ROOT}/analysis/src/build_tissue/binsize_test_archr"
LOG_DIR="${PROJECT_ROOT}/analysis/qsub_logs"

mkdir -p "${LOG_DIR}"

# Define combinations
TISSUES=("deepseq_488B" "deepseq_489" "lowseq_488B" "lowseq_489")
TILESIZES=(500 5000)
BINARIZE_OPTIONS=("FALSE" "TRUE")

echo "$(date '+%F %T') [start] Submitting 16 arrow variant creation jobs"
echo "Project: ${PROJECT_ROOT}"
echo "Logs: ${LOG_DIR}"
echo ""

JOB_COUNT=0

# Submit all combinations
for tissue in "${TISSUES[@]}"; do
  for tilesize in "${TILESIZES[@]}"; do
    for binarize in "${BINARIZE_OPTIONS[@]}"; do
      JOB_COUNT=$((JOB_COUNT + 1))
      
      # Create a PERMANENT qsub script for this job (don't delete after submit)
      LOGNAME="create_arrow_${tissue}_${tilesize}_${binarize}"
      QSUB_SCRIPT="${SCRIPT_DIR}/temp_qsub_${JOB_COUNT}.sh"
      
      cat > "${QSUB_SCRIPT}" <<QSUB_HEREDOC
#!/bin/bash
#\$ -N create_arrow
#\$ -l h_rt=24:00:00
#\$ -pe omp 8
#\$ -P paxlab
#\$ -l mem_per_core=8G
#\$ -o ${LOG_DIR}/${LOGNAME}.log
#\$ -e ${LOG_DIR}/${LOGNAME}.err

set -euo pipefail

# Initialize module system
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "\$profile_file" ]]; then
    . "\$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

module load R

PROJECT_ROOT="${PROJECT_ROOT}"
SCRIPT_DIR="\${PROJECT_ROOT}/analysis/src/build_tissue/binsize_test_archr"

echo "[$(date '+%F %T')] [start] Creating arrow: ${tissue} ${tilesize}bp binarize=${binarize}"

cd "\${PROJECT_ROOT}"

Rscript "\${SCRIPT_DIR}/create_arrow_variants.R" "${tissue}" "${tilesize}" "${binarize}"

exit_code=\$?
if [[ \${exit_code} -eq 0 ]]; then
  echo "[$(date '+%F %T')] [done] Successfully created arrow for ${tissue}"
else
  echo "[$(date '+%F %T')] [error] Failed with exit code: \${exit_code}"
  exit \${exit_code}
fi
QSUB_HEREDOC
      
      # Submit the job
      JOB_ID=$(qsub "${QSUB_SCRIPT}" 2>&1 | awk '{print $3}')
      echo "[${JOB_COUNT}] Job ${JOB_ID}: ${tissue} ${tilesize}bp binarize=${binarize} -> ${LOGNAME}.{log,err}"
      
      # Don't delete the qsub script - keep for reference
    done
  done
done

echo ""
echo "$(date '+%F %T') [done] Submitted ${JOB_COUNT} jobs"
echo ""
echo "Monitor with:"
echo "  qstat -u preshita"
echo "  tail -f ${LOG_DIR}/create_arrow_*.log"
echo ""
echo "Check log files when jobs complete:"
echo "  ls -lh ${LOG_DIR}/create_arrow_*.{log,err}"
