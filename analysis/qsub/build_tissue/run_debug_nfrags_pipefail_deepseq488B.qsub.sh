#!/bin/bash -l
#$ -P paxlab
#$ -N dbg_pipefail
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 1
#$ -l h_rt=00:10:00
#$ -l mem_per_core=4G
#$ -j n

set -euo pipefail
FILE="/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/fragments/deepseq_488B/deepseq_488B.fragments.sort.filtered.bed.gz"

echo "[debug2] host=$(hostname)"
echo "[debug2] start"

echo "[debug2] test A: gzip|head"
set +e
gzip -dc "$FILE" | head -n 3
stA=$?
set -e
echo "[debug2] status A=$stA"

echo "[debug2] test B: gzip|awk|head"
set +e
gzip -dc "$FILE" | awk 'NF>=4{bc=$4; sub(/-1$/, "", bc); n[bc]++} END{for(b in n) print b "\t" n[b]}' | head -n 3
stB=$?
set -e
echo "[debug2] status B=$stB"

echo "[debug2] test C: gzip|awk > file"
TMP_OUT="/projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/debug2_nfrags_${JOB_ID}.tsv"
set +e
gzip -dc "$FILE" | awk 'NF>=4{bc=$4; sub(/-1$/, "", bc); n[bc]++} END{for(b in n) print b "\t" n[b]}' > "$TMP_OUT"
stC=$?
set -e
echo "[debug2] status C=$stC"

if [[ -f "$TMP_OUT" ]]; then
  wc -l "$TMP_OUT"
  head -n 3 "$TMP_OUT"
fi

echo "[debug2] done"
