#!/bin/bash
# ============================================================================
# submit_metrics_recompute.sh
#
# Recompute metrics for all tissue x tilesize x binarize combinations.
# Arrow files already exist and will NOT be recreated.
# After all 16 metrics jobs complete, runs the comparison aggregation job.
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

echo "$(date '+%F %T') [start] Submitting 16 metrics recomputation jobs"
echo "Project: ${PROJECT_ROOT}"
echo "Logs: ${LOG_DIR}"
echo ""

JOB_COUNT=0
JOB_IDS=()

# Submit all combinations
for tissue in "${TISSUES[@]}"; do
  for tilesize in "${TILESIZES[@]}"; do
    for binarize in "${BINARIZE_OPTIONS[@]}"; do
      JOB_COUNT=$((JOB_COUNT + 1))

      # Create a qsub script for this job
      LOGNAME="recompute_metrics_${tissue}_${tilesize}_${binarize}"
      QSUB_SCRIPT="${SCRIPT_DIR}/temp_metrics_qsub_${JOB_COUNT}.sh"

      cat > "${QSUB_SCRIPT}" <<QSUB_HEREDOC
#!/bin/bash
#\$ -N recompute_metrics
#\$ -l h_rt=01:00:00
#\$ -pe omp 4
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

echo "[$(date '+%F %T')] [start] Recomputing metrics: ${tissue} ${tilesize}bp binarize=${binarize}"

cd "\${PROJECT_ROOT}"

Rscript "\${SCRIPT_DIR}/create_arrow_variants.R" "${tissue}" "${tilesize}" "${binarize}"

exit_code=\$?
if [[ \${exit_code} -eq 0 ]]; then
  echo "[$(date '+%F %T')] [done] Successfully recomputed metrics for ${tissue}"
else
  echo "[$(date '+%F %T')] [error] Failed with exit code: \${exit_code}"
  exit \${exit_code}
fi
QSUB_HEREDOC

      # Submit the job and capture ID
      JOB_ID=$(qsub "${QSUB_SCRIPT}" 2>&1 | awk '{print $3}')
      JOB_IDS+=("${JOB_ID}")
      echo "[${JOB_COUNT}] Job ${JOB_ID}: ${tissue} ${tilesize}bp binarize=${binarize} -> ${LOGNAME}.{log,err}"
    done
  done
done

echo ""
echo "$(date '+%F %T') [done] Submitted ${JOB_COUNT} metrics recomputation jobs"
echo ""

# Create the aggregation job that depends on all metrics jobs
AGGREGATION_LOGNAME="aggregate_metrics"
AGGREGATION_SCRIPT="${SCRIPT_DIR}/temp_aggregation_qsub.sh"
AGGREGATION_HOLD_JID=$(printf '%s,' "${JOB_IDS[@]}" | sed 's/,$//')

cat > "${AGGREGATION_SCRIPT}" <<AGGREGATION_HEREDOC
#!/bin/bash
#\$ -N aggregate_metrics
#\$ -l h_rt=00:30:00
#\$ -pe omp 1
#\$ -P paxlab
#\$ -l mem_per_core=4G
#\$ -o ${LOG_DIR}/${AGGREGATION_LOGNAME}.log
#\$ -e ${LOG_DIR}/${AGGREGATION_LOGNAME}.err
#\$ -hold_jid ${AGGREGATION_HOLD_JID}

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

echo "[$(date '+%F %T')] [start] Running sparsity comparison aggregation"

cd "\${PROJECT_ROOT}"

Rscript "\${SCRIPT_DIR}/compare_arrow_sparsity.R"

exit_code=\$?
if [[ \${exit_code} -eq 0 ]]; then
  echo "[$(date '+%F %T')] [done] Sparsity comparison completed"
else
  echo "[$(date '+%F %T')] [error] Comparison failed with exit code: \${exit_code}"
  exit \${exit_code}
fi
AGGREGATION_HEREDOC

AGGREGATION_JOB_ID=$(qsub "${AGGREGATION_SCRIPT}" 2>&1 | awk '{print $3}')
echo "[aggregation] Job ${AGGREGATION_JOB_ID}: compare_arrow_sparsity (depends on all metrics jobs)"

echo ""
echo "Monitor with:"
echo "  qstat -u preshita"
echo "  tail -f ${LOG_DIR}/recompute_metrics_*.log"
echo "  tail -f ${LOG_DIR}/aggregate_metrics.log"
echo ""
echo "When done, check results:"
echo "  cat ${LOG_DIR}/aggregate_metrics.log"
echo "  cat analysis/binsize_comparison/sparsity_recommendations.txt"
echo "  # View plots: analysis/binsize_comparison/sparsity_comparison_plots.pdf"
