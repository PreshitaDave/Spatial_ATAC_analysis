#!/bin/bash -l
#$ -N validate_489
#$ -P paxlab
#$ -l h_rt=02:00:00
#$ -l mem_per_core=4G
#$ -pe omp 2
#$ -j n
#$ -o analysis/qsub_logs/validate_489.$JOB_ID.out
#$ -e analysis/qsub_logs/validate_489.$JOB_ID.err
#$ -cwd

set -euo pipefail

safe_head_gz() {
  local gz_file="$1"
  gzip -dc "$gz_file" | sed -n '1,3p' || true
}

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
cd "$PROJECT_ROOT"

DEEP_BAM="/projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D1942/D1942_deepseq/processed_data/D1942_deepseq.bam"
LOW_BAM="/projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D1942/D01942_NG05827_lowseq/lowseq.sorted.bam"

DEEP_FRAG="${PROJECT_ROOT}/Data/01_inputs/fragments/deepseq_489.fragments.sort.filtered.bed.gz"
LOW_FRAG="${PROJECT_ROOT}/Data/01_inputs/fragments/lowseq_489.fragments.sort.filtered.bed.gz"

DEEP_BC="${PROJECT_ROOT}/Data/01_inputs/barcodes/tissue_barcodes/deepseq_489.barcodes.tsv"
LOW_BC="${PROJECT_ROOT}/Data/01_inputs/barcodes/tissue_barcodes/lowseq_489.barcodes.tsv"
SPATIAL="${PROJECT_ROOT}/Data/01_inputs/spatial/tissue_positions_list.csv"

TMP_DIR="${PROJECT_ROOT}/Data/07_tmp/validate_489_${JOB_ID}"
mkdir -p "$TMP_DIR"

echo "[$(date '+%F %T')] start validate_489 job=$JOB_ID host=$(hostname)" 
echo "[$(date '+%F %T')] project_root=$PROJECT_ROOT"

echo "=== File type checks ==="
file "$DEEP_BAM" "$LOW_BAM" "$DEEP_FRAG" "$LOW_FRAG" | cat

echo "=== Symlink checks ==="
ls -lh "$DEEP_FRAG" "$LOW_FRAG" | cat

echo "=== Fragment smoke test (first 3 lines) ==="
safe_head_gz "$DEEP_FRAG"
echo "---"
safe_head_gz "$LOW_FRAG"

echo "=== Build fragment-derived unique barcodes ==="
gzip -dc "$DEEP_FRAG" | awk 'NF>=4{bc=$4; sub(/-1$/, "", bc); print bc}' | sort -u > "$TMP_DIR/deep_frag_barcodes.tsv"
gzip -dc "$LOW_FRAG"  | awk 'NF>=4{bc=$4; sub(/-1$/, "", bc); print bc}' | sort -u > "$TMP_DIR/low_frag_barcodes.tsv"

sort "$DEEP_BC" > "$TMP_DIR/deep_489_barcodes.sorted.tsv"
sort "$LOW_BC"  > "$TMP_DIR/low_489_barcodes.sorted.tsv"

echo "=== Barcode counts ==="
echo "deep_489_barcode_file=$(wc -l < \"$TMP_DIR/deep_489_barcodes.sorted.tsv\")"
echo "low_489_barcode_file=$(wc -l < \"$TMP_DIR/low_489_barcodes.sorted.tsv\")"
echo "deep_frag_unique=$(wc -l < \"$TMP_DIR/deep_frag_barcodes.tsv\")"
echo "low_frag_unique=$(wc -l < \"$TMP_DIR/low_frag_barcodes.tsv\")"

echo "=== Overlap: barcode lists vs fragment-derived ==="
comm -12 "$TMP_DIR/deep_489_barcodes.sorted.tsv" "$TMP_DIR/deep_frag_barcodes.tsv" > "$TMP_DIR/deep_overlap.tsv"
comm -12 "$TMP_DIR/low_489_barcodes.sorted.tsv"  "$TMP_DIR/low_frag_barcodes.tsv"  > "$TMP_DIR/low_overlap.tsv"
echo "deep_overlap=$(wc -l < \"$TMP_DIR/deep_overlap.tsv\")"
echo "low_overlap=$(wc -l < \"$TMP_DIR/low_overlap.tsv\")"

echo "=== Deep vs Low 489 barcode concordance ==="
comm -12 "$TMP_DIR/deep_489_barcodes.sorted.tsv" "$TMP_DIR/low_489_barcodes.sorted.tsv" > "$TMP_DIR/deep_low_intersection.tsv"
echo "deep_low_intersection=$(wc -l < \"$TMP_DIR/deep_low_intersection.tsv\")"
echo "deep_only=$(comm -23 \"$TMP_DIR/deep_489_barcodes.sorted.tsv\" \"$TMP_DIR/low_489_barcodes.sorted.tsv\" | wc -l)"
echo "low_only=$(comm -13 \"$TMP_DIR/deep_489_barcodes.sorted.tsv\" \"$TMP_DIR/low_489_barcodes.sorted.tsv\" | wc -l)"

echo "=== Spatial overlap checks (in_tissue==1, inferred 489 by y<=4000) ==="
awk -F',' 'NF>=6 && $2==1 && $6+0<=4000 {print $1}' "$SPATIAL" | sort -u > "$TMP_DIR/spatial_489_barcodes.tsv"
echo "spatial_489=$(wc -l < \"$TMP_DIR/spatial_489_barcodes.tsv\")"
echo "deep_vs_spatial_489=$(comm -12 \"$TMP_DIR/deep_489_barcodes.sorted.tsv\" \"$TMP_DIR/spatial_489_barcodes.tsv\" | wc -l)"
echo "low_vs_spatial_489=$(comm -12 \"$TMP_DIR/low_489_barcodes.sorted.tsv\" \"$TMP_DIR/spatial_489_barcodes.tsv\" | wc -l)"

echo "=== Optional BAM header check (if samtools available) ==="
if command -v samtools >/dev/null 2>&1; then
  echo "samtools_found=yes"
  samtools view -H "$DEEP_BAM" | head -20 | cat
  echo "---"
  samtools view -H "$LOW_BAM" | head -20 | cat
else
  echo "samtools_found=no"
fi

echo "=== Artifacts ==="
echo "tmp_dir=$TMP_DIR"
echo "[$(date '+%F %T')] done"
