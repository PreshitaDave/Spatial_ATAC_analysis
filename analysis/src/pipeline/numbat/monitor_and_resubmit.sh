#!/bin/bash
# Auto-monitor ArchR job and resubmit NUMBAT when ready

cd /projectnb/paxlab/presh/projects/spatial_atac

ARCHIR_JOB=5649235
NUMBAT_SCRIPT="analysis/qsub/pipeline/run_numbat_lowseq_489.qsub.sh"
LOG_FILE="analysis/qsub_logs/build_tissue/numbat_autosubmit_5649235_monitor.log"

echo "[$(date +'%F %T')] Starting ArchR→NUMBAT auto-submit monitor" | tee "$LOG_FILE"
echo "[$(date +'%F %T')] Monitoring ArchR job: $ARCHIR_JOB" | tee -a "$LOG_FILE"

# Poll for ArchR job completion (check every 30 seconds, max 24 hours)
ELAPSED=0
MAX_WAIT=86400  # 24 hours in seconds
POLL_INTERVAL=30

while [ $ELAPSED -lt $MAX_WAIT ]; do
  # Check if ArchR job still exists in queue
  if ! qstat -j $ARCHIR_JOB &>/dev/null; then
    echo "[$(date +'%F %T')] ArchR job $ARCHIR_JOB completed (no longer in queue)" | tee -a "$LOG_FILE"
    break
  fi
  
  echo "[$(date +'%F %T')] ArchR job still running... (elapsed: ${ELAPSED}s)" | tee -a "$LOG_FILE"
  sleep $POLL_INTERVAL
  ELAPSED=$((ELAPSED + POLL_INTERVAL))
done

# Verify ArchR completed successfully
echo "[$(date +'%F %T')] Verifying ArchR projects were created..." | tee -a "$LOG_FILE"
if find Data/01_outputs/archR_objects -name "*_archR_project" -type d | grep -q .; then
  echo "[$(date +'%F %T')] ✓ ArchR projects found" | tee -a "$LOG_FILE"
  
  # Delete old NUMBAT job if still queued
  if qstat -j 5649110 &>/dev/null; then
    echo "[$(date +'%F %T')] Deleting old NUMBAT job 5649110..." | tee -a "$LOG_FILE"
    qdel 5649110
  fi
  
  # Submit NUMBAT job
  echo "[$(date +'%F %T')] Submitting NUMBAT job..." | tee -a "$LOG_FILE"
  JOB_OUTPUT=$(qsub "$NUMBAT_SCRIPT" 2>&1)
  
  if echo "$JOB_OUTPUT" | grep -q "Your job"; then
    NUMBAT_JOB=$(echo "$JOB_OUTPUT" | grep -oP '\d+' | head -1)
    echo "[$(date +'%F %T')] ✓ NUMBAT job submitted: $NUMBAT_JOB" | tee -a "$LOG_FILE"
    echo "[$(date +'%F %T')] Output: $JOB_OUTPUT" | tee -a "$LOG_FILE"
  else
    echo "[$(date +'%F %T')] ✗ Failed to submit NUMBAT job" | tee -a "$LOG_FILE"
    echo "[$(date +'%F %T')] Error: $JOB_OUTPUT" | tee -a "$LOG_FILE"
  fi
else
  echo "[$(date +'%F %T')] ✗ ArchR projects NOT found - ArchR job may have failed" | tee -a "$LOG_FILE"
  echo "[$(date +'%F %T')] Check log: tail -f analysis/qsub_logs/build_tissue/archR_qc_cluster_${ARCHIR_JOB}.log" | tee -a "$LOG_FILE"
fi

echo "[$(date +'%F %T')] Monitor script complete" | tee -a "$LOG_FILE"
