#!/bin/bash -l
#$ -P paxlab
#$ -N rebuild_fragments_bam
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 4
#$ -l h_rt=08:00:00
#$ -l mem_per_core=8G
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/rebuild_fragments_from_bam.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/rebuild_fragments_from_bam.$JOB_ID.err

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

# Deepseq: use 16bp rewritten BAMs if present (Bam_cb16)
DEEP_BAM_DIR="$PROJECT_ROOT/Data/03_intermediate/variant_calling/monopogen/variant_calling/deepseq/Bam_cb16"
if [[ -d "$DEEP_BAM_DIR" ]]; then
  log "Found deepseq Bam_cb16: $DEEP_BAM_DIR"
  BAM_FILES=("$DEEP_BAM_DIR"/*.bam)
else
  # fallback to original deepseq Bam
  BAM_FILES=("$PROJECT_ROOT/Data/03_intermediate/variant_calling/monopogen/variant_calling/deepseq/Bam"/*.bam)
fi

# Create fragments from BAMs (chromosome-split BAMs) into a single fragments file
TMP_FRAG="$OUT_FRAG_DIR/deepseq.fragments.from_bam.tmp.tsv"
: > "$TMP_FRAG"
for f in "${BAM_FILES[@]}"; do
  log "Processing $f"
  if command -v samtools >/dev/null 2>&1; then
    samtools view -F 0x4 "$f" \
      | awk '{OFS="\t"; chr=$3; start=$4-1; len=length($10); end=start+len; bc=""; for(i=12;i<=NF;i++){if($i~/^CB:Z:/){bc=substr($i,6); break}}; if(bc=="") next; print chr, start, end, bc "-1", 1}' \
      >> "$TMP_FRAG"
  else
    log "samtools not found on this node; aborting"
    exit 1
  fi
done

# sort and uniq counts approximate (keep as bed-like)
sort -k1,1 -k2,2n "$TMP_FRAG" | gzip -c > "$OUT_FRAG_DIR/deepseq.fragments.sort.filtered.bed.gz"
rm -f "$TMP_FRAG"
log "Wrote deepseq fragments to $OUT_FRAG_DIR/deepseq.fragments.sort.filtered.bed.gz"

# Validate barcode lengths and write nFrags cache
gzip -dc "$OUT_FRAG_DIR/deepseq.fragments.sort.filtered.bed.gz" | awk 'NF>=4{bc=$4; sub(/-1$/,"",bc); print length(bc)}' | sort | uniq -c > "$OUT_BARCODE_DIR/deepseq_fragment_barcode_lengths.txt"

awk 'NF>=4{bc=$4; sub(/-1$/,"",bc); n[bc]++} END{for(b in n) print b "\t" n[b]}' <(gzip -dc "$OUT_FRAG_DIR/deepseq.fragments.sort.filtered.bed.gz") > "$OUT_BARCODE_DIR/deepseq_nFrags_from_fragments.tsv"
gzip -f "$OUT_BARCODE_DIR/deepseq_nFrags_from_fragments.tsv"
log "Wrote deepseq nFrags cache to $OUT_BARCODE_DIR/deepseq_nFrags_from_fragments.tsv.gz"

# Lowseq: similar check but do not rewrite unless requested
LOW_BAM="$PROJECT_ROOT/Data/01_inputs/bam/lowseq_489.bam"
if [[ -f "$LOW_BAM" ]]; then
  log "Checking lowseq BAM $LOW_BAM"
  # quick barcode length check from BAM (sample first 10000 reads)
  if command -v samtools >/dev/null 2>&1; then
    samtools view "$LOW_BAM" | head -10000 | awk '{for(i=12;i<=NF;i++) if($i~/^CB:Z:/) {bc=substr($i,6); print length(bc); break}}' | sort | uniq -c > "$OUT_BARCODE_DIR/lowseq_sample_barcode_lengths.txt"
  fi
fi

log "done"
