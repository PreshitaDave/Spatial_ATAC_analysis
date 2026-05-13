#!/bin/bash -l
#$ -P paxlab
#$ -N rebuild_fragments_bam_tissue
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 4
#$ -l h_rt=12:00:00
#$ -l mem_per_core=8G
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/rebuild_fragments_from_bam_by_tissue.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/rebuild_fragments_from_bam_by_tissue.$JOB_ID.err

set -euo pipefail
PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
OUT_FRAG_DIR="$PROJECT_ROOT/Data/01_inputs/fragments"
OUT_BARCODE_DIR="$PROJECT_ROOT/Data/01_inputs/barcodes/tissue_barcodes"
mkdir -p "$OUT_FRAG_DIR" "$OUT_BARCODE_DIR"

# Load samtools
if [[ -f /etc/profile.d/modules.sh ]]; then
  source /etc/profile.d/modules.sh
fi
module load samtools || true

log(){ echo "[$(date '+%F %T')] $*"; }
TS=$(date +%s)

# Datasets and tissues
datasets=(deepseq lowseq)
tissues=(488B 489)

for dataset in "${datasets[@]}"; do
  log "Processing dataset: $dataset"

  if [[ "$dataset" == "deepseq" ]]; then
    # Prefer rewritten 16bp BAMs
    DEEP_BAM_DIR="$PROJECT_ROOT/Data/03_intermediate/variant_calling/monopogen/variant_calling/deepseq/Bam_cb16"
    if [[ -d "$DEEP_BAM_DIR" ]]; then
      log "Found deepseq Bam_cb16: $DEEP_BAM_DIR"
      BAM_FILES=("$DEEP_BAM_DIR"/*.bam)
    else
      BAM_FILES=("$PROJECT_ROOT/Data/03_intermediate/variant_calling/monopogen/variant_calling/deepseq/Bam"/*.bam)
    fi
  else
    # lowseq: use per-tissue BAM symlinks in Data/01_inputs/bam
    BAM_FILES=()
    for t in "${tissues[@]}"; do
      p="$PROJECT_ROOT/Data/01_inputs/bam/${dataset}_${t}.bam"
      if [[ -f "$p" ]]; then
        BAM_FILES+=("$p")
      fi
    done
  fi

  for tissue in "${tissues[@]}"; do
    BARCODE_LIST="$OUT_BARCODE_DIR/${dataset}_${tissue}.barcodes.tsv"
    if [[ ! -f "$BARCODE_LIST" ]]; then
      log "Barcode list missing for ${dataset}_${tissue}: $BARCODE_LIST — skipping tissue"
      continue
    fi

    OUT_FRAG="$OUT_FRAG_DIR/${dataset}_${tissue}.fragments.from_bam.${TS}.bed.gz"
    TMP_FRAG="$OUT_FRAG_DIR/${dataset}_${tissue}.fragments.from_bam.${TS}.tmp"
    : > "$TMP_FRAG"

    log "Building fragments for ${dataset}_${tissue} into $OUT_FRAG (using ${#BAM_FILES[@]} BAMs)"

    for f in "${BAM_FILES[@]}"; do
      if [[ ! -f "$f" ]]; then
        log "BAM not found: $f — skipping"
        continue
      fi
      log "Processing BAM: $f"
      # Extract reads, capture CB tag, produce bed-like fragment lines
      samtools view -F 0x4 "$f" \
        | awk -v BARC="$BARCODE_LIST" 'BEGIN{while((getline line < BARC)>0) bar[line]=1} {bc=""; for(i=12;i<=NF;i++) if($i~/^CB:Z:/){bc=substr($i,6); break} if(bc=="") next; if(!(bc in bar)) next; chr=$3; start=$4-1; len=length($10); end=start+len; print chr"\t"start"\t"end"\t"bc""-1"\t"1}' >> "$TMP_FRAG"
    done

    # Sort and compress to final path (do not overwrite existing files)
    log "Sorting and compressing fragments to $OUT_FRAG"
    sort -k1,1 -k2,2n "$TMP_FRAG" | gzip -c > "$OUT_FRAG"
    rm -f "$TMP_FRAG"

    # Barcode length summary
    gzip -dc "$OUT_FRAG" | awk 'NF>=4{bc=$4; sub(/-1$/,"",bc); print length(bc)}' | sort | uniq -c > "$OUT_BARCODE_DIR/${dataset}_${tissue}_fragment_barcode_lengths.${TS}.txt"

    # Write nFrags cache
    gzip -dc "$OUT_FRAG" | awk 'NF>=4{bc=$4; sub(/-1$/,"",bc); n[bc]++} END{for(b in n) print b "\t" n[b]}' > "$OUT_BARCODE_DIR/${dataset}_${tissue}_nFrags_from_bam.${TS}.tsv"
    gzip -f "$OUT_BARCODE_DIR/${dataset}_${tissue}_nFrags_from_bam.${TS}.tsv"
    log "Wrote nFrags cache and barcode length summary for ${dataset}_${tissue}"
  done

done

# Quick lowseq existence check for tissue 489
if [[ -f "$PROJECT_ROOT/Data/01_inputs/fragments/lowseq_489.fragments.sort.filtered.bed.gz" ]]; then
  log "lowseq_489 fragment file exists: yes"
else
  log "lowseq_489 fragment file exists: no"
fi

log "done"
