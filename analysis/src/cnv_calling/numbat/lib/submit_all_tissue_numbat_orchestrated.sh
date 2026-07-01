#!/bin/bash
################################################################################
# NUMBAT ORCHESTRATED SUBMISSION SCRIPT
# Master controller for deepseq and lowseq NUMBAT pipeline
# Handles all job dependencies and serial execution requirements
#
# Job Naming Convention: numbat_prep_{dataset}_{tissue}
#   Examples:
#     - numbat_prep_lowseq_488B
#     - numbat_prep_deepseq_489
#     - numbat_analysis_lowseq_combined
################################################################################

set -eo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
cd "$PROJECT_ROOT"

echo "================================================================================"
echo "NUMBAT ORCHESTRATED SUBMISSION - $(date +'%F %T')"
echo "================================================================================"
echo ""

# ============================================================================
# STEP 1: Kill buggy lowseq jobs and clean contaminated output
# ============================================================================
echo "[STEP 1] Cleaning up buggy lowseq jobs..."
echo ""

# Kill jobs if they exist
for JOB in 5692040 5692041; do
  if qstat -j $JOB >/dev/null 2>&1; then
    echo "  Killing job $JOB..."
    qdel $JOB 2>/dev/null || true
    sleep 2
  fi
done

# Clean contaminated output
CONTAMINATED_DIR="$PROJECT_ROOT/Data/04_analysis/cnv/numbat/inputs/alleles/pileup"
if [[ -d "$CONTAMINATED_DIR" ]]; then
  echo "  Removing contaminated pileup directories..."
  rm -rf "$CONTAMINATED_DIR"/lowseq_*_atac 2>/dev/null || true
  echo "  ✓ Cleaned"
fi
echo ""

# ============================================================================
# STEP 2: Check current deepseq status
# ============================================================================
echo "[STEP 2] Checking deepseq_488B status..."
if qstat -j 5692082 >/dev/null 2>&1; then
  echo "  ✓ Job 5692082 (numbat_prep_deepseq_488B) still RUNNING"
  DEEPSEQ_488B_RUNNING=true
else
  echo "  ✗ Job 5692082 not found or completed"
  DEEPSEQ_488B_RUNNING=false
fi
echo ""

# ============================================================================
# STEP 3: Submit lowseq jobs (serial: 488B then 489)
# ============================================================================
echo "[STEP 3] Submitting lowseq prepare jobs (SERIAL: 488B → 489)..."
echo ""

echo "  Submitting: numbat_prep_lowseq_488B"
LOWSEQ_488B_JOB=$(qsub "$PROJECT_ROOT/analysis/src/cnv_calling/numbat/prepare_numbat_inputs_lowseq_488B.qsub.sh" | awk '{print $1}')
echo "    Job ID: $LOWSEQ_488B_JOB"
sleep 1

echo ""
echo "  Submitting: numbat_prep_lowseq_489 (depends on $LOWSEQ_488B_JOB)"
LOWSEQ_489_JOB=$(qsub -hold_jid $LOWSEQ_488B_JOB "$PROJECT_ROOT/analysis/src/cnv_calling/numbat/prepare_numbat_inputs_lowseq_489.qsub.sh" | awk '{print $1}')
echo "    Job ID: $LOWSEQ_489_JOB (will start after $LOWSEQ_488B_JOB)"
echo ""

# ============================================================================
# STEP 4: Submit deepseq_489 prepare (depends on deepseq_488B if still running)
# ============================================================================
echo "[STEP 4] Submitting deepseq prepare jobs..."
echo ""

if [[ "$DEEPSEQ_488B_RUNNING" == true ]]; then
  echo "  deepseq_488B (5692082) still running"
  echo "  Submitting: numbat_prep_deepseq_489 (depends on 5692082)"
  DEEPSEQ_489_JOB=$(qsub -hold_jid 5692082 "$PROJECT_ROOT/analysis/src/cnv_calling/numbat/prepare_numbat_inputs_deepseq_489.qsub.sh" | awk '{print $1}')
  echo "    Job ID: $DEEPSEQ_489_JOB (will start after 5692082)"
else
  echo "  ✗ deepseq_488B not running - cannot submit deepseq_489 with dependency"
  echo "  Will need manual submission when 488B completes"
  DEEPSEQ_489_JOB=""
fi
echo ""

# ============================================================================
# STEP 5: Summary
# ============================================================================
echo "================================================================================"
echo "SUBMISSION SUMMARY"
echo "================================================================================"
echo ""
echo "LOWSEQ JOBS:"
echo "  numbat_prep_lowseq_488B       Job $LOWSEQ_488B_JOB"
echo "  numbat_prep_lowseq_489        Job $LOWSEQ_489_JOB (→ depends on $LOWSEQ_488B_JOB)"
echo ""
echo "DEEPSEQ JOBS:"
echo "  numbat_prep_deepseq_488B      Job 5692082 (already running)"
if [[ -n "$DEEPSEQ_489_JOB" ]]; then
  echo "  numbat_prep_deepseq_489       Job $DEEPSEQ_489_JOB (→ depends on 5692082)"
else
  echo "  numbat_prep_deepseq_489       NEEDS MANUAL SUBMISSION"
fi
echo ""
echo "EXECUTION ORDER:"
echo "  1. numbat_prep_deepseq_488B (5692082) running ~80%"
echo "     └─ ETA complete: 11:00-11:30 AM"
echo ""
echo "  2. numbat_prep_lowseq_488B ($LOWSEQ_488B_JOB) submitted"
echo "     └─ ETA complete: 12:00-1:00 PM"
echo ""
echo "  3. numbat_prep_deepseq_489 ($DEEPSEQ_489_JOB) depends on 5692082"
echo "     └─ ETA complete: 2:00-2:30 PM"
echo ""
echo "  4. numbat_prep_lowseq_489 ($LOWSEQ_489_JOB) depends on $LOWSEQ_488B_JOB"
echo "     └─ ETA complete: 2:00-3:00 PM"
echo ""
echo "  5. Analysis jobs (488B) start when prep completes"
echo "  6. Analysis jobs (489) start when prep completes"
echo "  7. Combined analysis jobs run after both analyses complete"
echo ""
echo "================================================================================"
echo "MONITOR JOBS:"
echo "  qstat -u preshita"
echo "  OR"
echo "  watch -n 5 'qstat -u preshita | grep numbat'"
echo "================================================================================"
echo ""
echo "✓ Master submission complete - $(date +'%T')"
