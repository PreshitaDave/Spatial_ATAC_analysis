#!/bin/bash
# ============================================================================
# submit_all_v2.sh
#
# Submit the QC-corrected (v2) arrow creation + ArchR project build for every
# remaining tissue/tilesize/binarize combination in binsize_comparison.
# deepseq_488B_5000bp_binarizeFALSE is already done and is skipped here.
#
# Each combo: create_arrow (16 cores) -hold_jid-> build_arch (8 cores)
# All combos submitted up front; SGE's scheduler queues/runs as capacity allows.
#
# Usage: bash submit_all_v2.sh
# ============================================================================

set -eu

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
QSUB_DIR="$PROJECT_ROOT/analysis/qsub/build_tissue"
LOG_DIR="$PROJECT_ROOT/analysis/qsub_logs/build_tissue"

mkdir -p "$LOG_DIR"

log_msg() {
  echo "[$(date '+%F %T')] $1"
}

# tissue:tilesize:binarize combos already done (skip these)
DONE_COMBOS=("deepseq_488B:5000:FALSE" "deepseq_488B:500:FALSE")

is_done() {
  local combo="$1:$2:$3"
  for d in "${DONE_COMBOS[@]}"; do
    [[ "$d" == "$combo" ]] && return 0
  done
  return 1
}

TISSUES=("deepseq_488B" "deepseq_489" "lowseq_488B" "lowseq_489")
TILESIZES=(500 5000)
BINARIZE_OPTIONS=("FALSE" "TRUE")

submitted=0

for tissue in "${TISSUES[@]}"; do
  for tilesize in "${TILESIZES[@]}"; do
    for binarize in "${BINARIZE_OPTIONS[@]}"; do
      if is_done "$tissue" "$tilesize" "$binarize"; then
        log_msg "SKIP (already done): $tissue $tilesize $binarize"
        continue
      fi

      ARROW_LOG="$LOG_DIR/create_arrow_v2_${tissue}_${tilesize}_${binarize}_\$JOB_ID.log"
      BUILD_LOG="$LOG_DIR/build_archr_v2_${tissue}_${tilesize}_${binarize}_\$JOB_ID.log"

      JOB1=$(qsub -N "carr_${tissue}_${tilesize}${binarize}" \
        -o "$ARROW_LOG" \
        "$QSUB_DIR/run_create_arrow_variants_v2.qsub.sh" "$tissue" "$tilesize" "$binarize" \
        | awk '{print $3}')

      JOB2=$(qsub -N "barc_${tissue}_${tilesize}${binarize}" \
        -hold_jid "$JOB1" \
        -o "$BUILD_LOG" \
        "$QSUB_DIR/run_build_archr_variant_project_v2.qsub.sh" "$tissue" "$tilesize" "$binarize" \
        | awk '{print $3}')

      log_msg "Submitted: $tissue $tilesize $binarize  -> arrow=$JOB1 build=$JOB2 (held on $JOB1)"
      submitted=$((submitted + 1))
    done
  done
done

log_msg ""
log_msg "===== Submission complete: $submitted combos (x2 jobs each) ====="
log_msg "Monitor with: qstat -u $(whoami)"
