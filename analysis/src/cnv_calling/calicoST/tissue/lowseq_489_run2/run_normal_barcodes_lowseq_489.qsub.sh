#!/bin/bash -l
#$ -P paxlab
#$ -N normal_barcodes_lowseq_489
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 4
#$ -l mem_per_core=16G
#$ -l h_rt=02:00:00
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/normal_barcodes_lowseq_489_$JOB_ID.log
#$ -m bea
#$ -M preshita@bu.edu

set -euo pipefail

module load R

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
RSCRIPT="Rscript"
SCRIPT="${PROJECT_ROOT}/analysis/src/build_tissue/binsize_test_archr/9_extract_normal_barcodes_calicost.R"
OUT="${PROJECT_ROOT}/analysis/binsize_comparison/normal_barcodes/lowseq_489_normal_barcodes.csv"

echo "=== Normal barcode extraction: lowseq_489 ==="
echo "Job ID: ${JOB_ID}"
echo "Start: $(date)"
echo "Host: $(hostname)"

${RSCRIPT} "${SCRIPT}" lowseq_489

echo ""
echo "Output: ${OUT}"
wc -l "${OUT}" 2>/dev/null && echo "Barcodes written successfully." || echo "WARNING: output file not found"
echo "End: $(date)"
