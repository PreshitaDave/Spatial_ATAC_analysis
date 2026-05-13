#!/bin/bash -l
#$ -N validate_tt
#$ -P paxlab
#$ -l h_rt=02:00:00
#$ -l mem_per_core=4G
#$ -pe omp 2
#$ -j n
#$ -o analysis/qsub_logs/validate_tissue_tech.$JOB_ID.out
#$ -e analysis/qsub_logs/validate_tissue_tech.$JOB_ID.err
#$ -cwd

set -euo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
cd "$PROJECT_ROOT"

BARCODE_DIR="$PROJECT_ROOT/Data/01_inputs/barcodes/tissue_barcodes"
FRAG_DIR="$PROJECT_ROOT/Data/01_inputs/fragments"
BAM_DIR="$PROJECT_ROOT/Data/01_inputs/bam"
TMP_DIR="$PROJECT_ROOT/Data/07_tmp/validate_tissue_tech_${JOB_ID}"
OUT_TSV="$BARCODE_DIR/input_matrix_tissue_tech.tsv"

mkdir -p "$TMP_DIR"
mkdir -p "$BARCODE_DIR"

echo "[$(date '+%F %T')] start validate_tissue_tech job=$JOB_ID host=$(hostname)"

# Ensure canonical 488B aliases exist (no data copy, only links)
if [[ ! -e "$FRAG_DIR/deepseq_488B.fragments.sort.filtered.bed.gz" && -e "$FRAG_DIR/deepseq.fragments.sort.filtered.bed.gz" ]]; then
  ln -s "$FRAG_DIR/deepseq.fragments.sort.filtered.bed.gz" "$FRAG_DIR/deepseq_488B.fragments.sort.filtered.bed.gz"
  [[ -e "$FRAG_DIR/deepseq.fragments.sort.filtered.bed.gz.tbi" ]] && \
    ln -sf "$FRAG_DIR/deepseq.fragments.sort.filtered.bed.gz.tbi" "$FRAG_DIR/deepseq_488B.fragments.sort.filtered.bed.gz.tbi"
  echo "[$(date '+%F %T')] linked deepseq_488B.fragments -> deepseq.fragments"
fi

if [[ ! -e "$FRAG_DIR/lowseq_488B.fragments.sort.filtered.bed" && -e "$FRAG_DIR/lowseq.fragments.sort.filtered.bed" ]]; then
  ln -s "$FRAG_DIR/lowseq.fragments.sort.filtered.bed" "$FRAG_DIR/lowseq_488B.fragments.sort.filtered.bed"
  echo "[$(date '+%F %T')] linked lowseq_488B.fragments -> lowseq.fragments"
fi

printf "tech\ttissue\tbarcode_file\tbarcode_n\tfragment_file\tfragment_exists\tbam_file\tbam_exists\tfragment_unique_barcodes\toverlap_n\toverlap_pct\taction\n" > "$OUT_TSV"

process_combo() {
  local tech="$1"
  local tissue="$2"

  local bc_file="$BARCODE_DIR/${tech}_${tissue}.barcodes.tsv"
  local frag_file=""
  local bam_file=""
  local action="ok"

  if [[ "$tech" == "deepseq" && "$tissue" == "488B" ]]; then
    frag_file="$FRAG_DIR/deepseq_488B.fragments.sort.filtered.bed.gz"
    bam_file="$BAM_DIR/deepseq_488B.bam"
  elif [[ "$tech" == "deepseq" && "$tissue" == "489" ]]; then
    frag_file="$FRAG_DIR/deepseq_489.fragments.sort.filtered.bed.gz"
    bam_file="$BAM_DIR/deepseq_489.bam"
  elif [[ "$tech" == "lowseq" && "$tissue" == "488B" ]]; then
    frag_file="$FRAG_DIR/lowseq_488B.fragments.sort.filtered.bed"
    bam_file="$BAM_DIR/lowseq_488B.bam"
  else
    frag_file="$FRAG_DIR/lowseq_489.fragments.sort.filtered.bed.gz"
    bam_file="$BAM_DIR/lowseq_489.bam"
  fi

  local bc_n=0
  local frag_exists=0
  local bam_exists=0
  local frag_uniq_n=0
  local overlap_n=0
  local overlap_pct="0.00"

  [[ -e "$bc_file" ]] && bc_n=$(wc -l < "$bc_file") || action="missing_barcode"
  [[ -e "$frag_file" ]] && frag_exists=1 || action="missing_fragment"
  [[ -e "$bam_file" ]] && bam_exists=1

  if [[ "$frag_exists" -eq 0 && "$bam_exists" -eq 1 ]]; then
    # Fallback attempt: only if tooling is present. Keep this small and explicit.
    if command -v samtools >/dev/null 2>&1 && command -v bedtools >/dev/null 2>&1; then
      local out_frag_gz="$frag_file"
      mkdir -p "$(dirname "$out_frag_gz")"
      bedtools bamtobed -i "$bam_file" | awk 'NF>=4{print $1"\t"$2"\t"$3"\t"$4"\t1}' | gzip -c > "$out_frag_gz"
      if [[ -s "$out_frag_gz" ]]; then
        frag_exists=1
        action="derived_from_bam"
      else
        action="derive_failed"
      fi
    else
      action="missing_fragment_no_tools"
    fi
  fi

  if [[ "$bc_n" -gt 0 && "$frag_exists" -eq 1 ]]; then
    local bc_sorted="$TMP_DIR/${tech}_${tissue}.bc.sorted.tsv"
    local frag_sorted="$TMP_DIR/${tech}_${tissue}.frag.sorted.tsv"

    sort "$bc_file" > "$bc_sorted"
    if [[ "$frag_file" == *.gz ]]; then
      gzip -dc "$frag_file" | awk 'NF>=4{bc=$4; sub(/-1$/, "", bc); print bc}' | sort -u > "$frag_sorted"
    else
      cat "$frag_file" | awk 'NF>=4{bc=$4; sub(/-1$/, "", bc); print bc}' | sort -u > "$frag_sorted"
    fi

    frag_uniq_n=$(wc -l < "$frag_sorted")
    overlap_n=$(comm -12 "$bc_sorted" "$frag_sorted" | wc -l)
    overlap_pct=$(awk -v a="$overlap_n" -v b="$bc_n" 'BEGIN{if (b>0) printf "%.2f", (100*a)/b; else printf "0.00"}')
  fi

  printf "%s\t%s\t%s\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%s\t%s\n" \
    "$tech" "$tissue" "$bc_file" "$bc_n" "$frag_file" "$frag_exists" "$bam_file" "$bam_exists" \
    "$frag_uniq_n" "$overlap_n" "$overlap_pct" "$action" >> "$OUT_TSV"
}

process_combo deepseq 488B
process_combo deepseq 489
process_combo lowseq 488B
process_combo lowseq 489

echo "[$(date '+%F %T')] wrote $OUT_TSV"
cat "$OUT_TSV"
echo "[$(date '+%F %T')] done"
