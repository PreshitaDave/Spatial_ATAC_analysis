#!/bin/bash -l
#$ -P paxlab
#$ -N frag_489
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=8G
#$ -l h_rt=12:00:00
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/frag_489.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/frag_489.$JOB_ID.err
#$ -j n

set -euo pipefail
module load samtools

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
BARCODE_DIR="${PROJECT_ROOT}/Data/alleloscope/barcodes"
DATA_DIR="${PROJECT_ROOT}/Data"
TISSUE_FRAG_DIR="${DATA_DIR}/variant_calling/lowseq/tissue"

log() { echo "[$(date '+%F %T')] $*"; }

# NOTE: lowseq.fragments.sort.filtered.bed only contains 488B reads.
# Tissue 489 fragments live in variant_calling/lowseq/tissue/lowseq_489.fragments.tsv.gz
# We use that as the source and reformat to BED for Alleloscope.

# Build -1-suffixed barcode lookup files for awk
log "Building 489 barcode lookup files"
awk '{print $1 "-1"}' "${BARCODE_DIR}/deepseq_489.barcodes.tsv" | sort > /tmp/deepseq_489_bc1.txt
log "Deepseq 489 barcodes: $(wc -l < /tmp/deepseq_489_bc1.txt)"

# ---- lowseq 489: source is tissue-specific fragments TSV ----
SRC_LOWSEQ="${TISSUE_FRAG_DIR}/lowseq_489.fragments.tsv.gz"
OUT_LOWSEQ="${DATA_DIR}/lowseq_489.fragments.sort.filtered.bed"
log "Source lowseq 489: ${SRC_LOWSEQ} ($(ls -lh ${SRC_LOWSEQ} | awk '{print $5}'))"
log "Writing lowseq 489 fragments BED -> ${OUT_LOWSEQ}.gz"
# fragments.tsv.gz is already chr/start/end/barcode/count format — just recompress
zcat "${SRC_LOWSEQ}" \
    | sort -k1,1 -k2,2n -k3,3n --parallel="${NSLOTS:-8}" -S 8G \
    | bgzip -@ "${NSLOTS:-8}" > "${OUT_LOWSEQ}.gz"
tabix -p bed "${OUT_LOWSEQ}.gz"
log "Lowseq 489 fragment lines: $(zcat ${OUT_LOWSEQ}.gz | wc -l)"
log "Wrote: ${OUT_LOWSEQ}.gz + .tbi"

# ---- deepseq 489 ----
OUT_DEEPSEQ="${DATA_DIR}/deepseq_489.fragments.sort.filtered.bed"
log "Filtering deepseq fragments to tissue 489 -> ${OUT_DEEPSEQ}.gz"
zcat "${DATA_DIR}/deepseq.fragments.sort.filtered.bed.gz" \
    | awk 'NR==FNR { bc[$1]=1; next } $4 in bc' \
    /tmp/deepseq_489_bc1.txt - \
    > "${OUT_DEEPSEQ}"
log "Deepseq 489 fragments: $(wc -l < "${OUT_DEEPSEQ}")"

log "Sorting and bgzip-compressing deepseq 489 fragments"
sort -k1,1 -k2,2n -k3,3n --parallel="${NSLOTS:-8}" -S 8G "${OUT_DEEPSEQ}" \
    | bgzip -@ "${NSLOTS:-8}" > "${OUT_DEEPSEQ}.gz"
tabix -p bed "${OUT_DEEPSEQ}.gz"
rm -f "${OUT_DEEPSEQ}"
log "Wrote: ${OUT_DEEPSEQ}.gz + .tbi"

rm -f /tmp/lowseq_489_bc1.txt /tmp/deepseq_489_bc1.txt
log "DONE tissue 489 fragments"
