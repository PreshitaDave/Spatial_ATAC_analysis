#!/bin/bash
set -e

PROJECT_ROOT="${PROJECT_ROOT:=/projectnb/paxlab/presh/projects/spatial_atac}"
LOG_DIR="${PROJECT_ROOT}/analysis/qsub_logs"

cd "${PROJECT_ROOT}"

{
  echo "[$(date '+%F %T')] [start] job_id=$JOB_ID host=$(hostname)"
  echo "[$(date '+%F %T')] [start] cwd=$PWD"
  
  DEEPSEQ_FRAG="${PROJECT_ROOT}/Data/01_inputs/fragments/deepseq.fragments.sort.filtered.bed.gz"
  
  echo "[$(date '+%F %T')] [step] === Fragment file checks ==="
  ls -lh "${DEEPSEQ_FRAG}" 2>&1
  file "${DEEPSEQ_FRAG}" 2>&1
  
  echo "[$(date '+%F %T')] [step] === Fragment header (first 3 lines) ==="
  gzip -dc "${DEEPSEQ_FRAG}" 2>&1 | head -3
  
  echo "[$(date '+%F %T')] [step] === Barcode length distribution (sampling 100k lines) ==="
  gzip -dc "${DEEPSEQ_FRAG}" 2>&1 | head -100000 \
    | awk 'NF>=4{bc=$4; sub(/-1$/,"",bc); print length(bc)}' \
    | sort | uniq -c
  
  echo "[$(date '+%F %T')] [step] === Total unique barcodes (full file) ==="
  gzip -dc "${DEEPSEQ_FRAG}" 2>&1 \
    | awk 'NF>=4{bc=$4; sub(/-1$/,"",bc); n[bc]++} END{print "Unique barcodes: " length(n); print "Total rows: " NR}' 
  
  echo "[$(date '+%F %T')] [step] === BAM file checks ==="
  DEEPSEQ_BAM_ORIG="${PROJECT_ROOT}/Data/03_intermediate/variant_calling/monopogen/variant_calling/deepseq/Bam/chr1.filter.targeted.bam"
  DEEPSEQ_BAM_CB16="${PROJECT_ROOT}/Data/03_intermediate/variant_calling/monopogen/variant_calling/deepseq/Bam_cb16/deepseq_chr1.filter.bam"
  
  ls -lh "${DEEPSEQ_BAM_ORIG}" 2>&1 | tail -1
  ls -lh "${DEEPSEQ_BAM_CB16}" 2>&1 | tail -1
  
  echo "[$(date '+%F %T')] [done] diagnostic complete"
  
} 2>&1 | tee "${LOG_DIR}/diagnose_deepseq_barcodes.$(date '+%s').out"
