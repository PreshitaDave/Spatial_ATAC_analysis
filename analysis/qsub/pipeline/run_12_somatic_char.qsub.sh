#!/bin/bash -l
#$ -N somatic_char
#$ -P paxlab
#$ -l h_rt=12:00:00
#$ -l mem_per_core=8G
#$ -pe omp 8
#$ -j n
#$ -o analysis/qsub_logs/somatic_char.$JOB_ID.out
#$ -e analysis/qsub_logs/somatic_char.$JOB_ID.err
#$ -cwd

set -euo pipefail

echo "[$(date '+%F %T')] Job $JOB_ID host=$(hostname) NSLOTS=$NSLOTS" \
  >> analysis/qsub_logs/somatic_char.$JOB_ID.out

if ! command -v module >/dev/null 2>&1; then
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi
module load R

# ── Script-11 must have run first ──────────────────────────────────────────
# Outputs consumed: analysis/comparison/tissue_variants/{dataset}_{set}.tsv
#
# Run script 11 first if outputs are missing:
#   qsub analysis/qsub/pipeline/run_11_comparing_tissue_variants.qsub.sh
# ──────────────────────────────────────────────────────────────────────────

export COMP_DIR="${COMP_DIR:-analysis/comparison/tissue_variants}"
export OUT_DIR="${OUT_DIR:-analysis/comparison/somatic_char}"
export SKIP_TRINUC="${SKIP_TRINUC:-1}"   # skip BSgenome trinucleotide (too slow on 3M+ variants)
export DATASETS="${DATASETS:-deepseq,lowseq}"
export CHR_START="${CHR_START:-1}"
export CHR_END="${CHR_END:-22}"
export N_WORKERS="${N_WORKERS:-$NSLOTS}"
export UPSTREAM_BP="${UPSTREAM_BP:-5000}"

Rscript analysis/src/12_somatic_snv_characterization.R
