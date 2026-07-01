#!/bin/bash
################################################################################
# submit_deepseq_prepare_sequential.sh
#
# Submit deepseq prepare jobs sequentially to prevent data contamination
# This script ensures deepseq_489 only starts AFTER deepseq_488B completes
#
################################################################################

set -eo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
cd "$PROJECT_ROOT"

echo "=========================================="
echo "SUBMITTING DEEPSEQ PREPARE JOBS (SEQUENTIAL)"
echo "=========================================="
echo ""

echo "[STEP 1] Submit deepseq_488B prepare job (no dependency)..."
JOB_488B=$(qsub analysis/src/cnv_calling/numbat/prepare_numbat_inputs_deepseq_488B.qsub.sh)
JOB_ID_488B=$(echo "$JOB_488B" | grep -oE "[0-9]+")
echo "  Job ID: $JOB_ID_488B"
echo "  Status: Queued"
echo ""

echo "[STEP 2] Submit deepseq_489 prepare job (depends on 488B completing)..."
DEEPSEQ_489_SCRIPT="analysis/src/cnv_calling/numbat/prepare_numbat_inputs_deepseq_489.qsub.sh"

# Add dependency to the qsub call
qsub -hold_jid "$JOB_ID_488B" "$DEEPSEQ_489_SCRIPT" > /tmp/qsub_489.out
JOB_ID_489=$(cat /tmp/qsub_489.out | grep -oE "[0-9]+")
echo "  Job ID: $JOB_ID_489"
echo "  Status: Held (waiting for job $JOB_ID_488B to complete)"
echo ""

echo "[MONITORING]"
echo "  deepseq_488B will run first (Job $JOB_ID_488B)"
echo "  deepseq_489 will start automatically when 488B finishes (Job $JOB_ID_489)"
echo ""
echo "Watch status with:"
echo "  qstat -u preshita"
echo ""
