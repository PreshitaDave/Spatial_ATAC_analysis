#!/bin/bash -l
#$ -P paxlab
#$ -N dbg_nfrags_awk
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 1
#$ -l h_rt=00:10:00
#$ -l mem_per_core=4G
#$ -j n

set -euo pipefail

FILE="/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/fragments/deepseq_488B/deepseq_488B.fragments.sort.filtered.bed.gz"
TMP_OUT="/projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/debug_nfrags_${JOB_ID}.tsv"

echo "[debug] host=$(hostname)"
echo "[debug] file=$FILE"
ls -lh "$FILE"

echo "[debug] first fragment lines"
gzip -dc "$FILE" | head -n 3

echo "[debug] awk sample output"
gzip -dc "$FILE" | awk 'NF>=4{bc=$4; sub(/-1$/, "", bc); n[bc]++} END{for(b in n) print b "\t" n[b]}' | head -n 5

echo "[debug] writing full awk output to $TMP_OUT"
gzip -dc "$FILE" | awk 'NF>=4{bc=$4; sub(/-1$/, "", bc); n[bc]++} END{for(b in n) print b "\t" n[b]}' > "$TMP_OUT"

echo "[debug] awk output line count"
wc -l "$TMP_OUT"

echo "[debug] awk output first lines"
head -n 3 "$TMP_OUT"
