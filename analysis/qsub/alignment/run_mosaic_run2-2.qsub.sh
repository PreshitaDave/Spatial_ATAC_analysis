#!/bin/bash -l
#$ -P paxlab
#$ -N mosaic_run2
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=16G
#$ -l h_rt=06:00:00
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/run_mosaic_run2-2_$JOB_ID.log
#$ -m bea
#$ -M preshita@bu.edu

set -euo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
ALIGN_DIR="${PROJECT_ROOT}/analysis/src/alignment"
ENV_PY="/projectnb/paxlab/presh/env/conda_env/mosaicfield_env/bin/python3"
NBCONVERT="/projectnb/paxlab/presh/env/conda_env/mosaicfield_env/bin/jupyter-nbconvert"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_NB="mosaic_run2-2_executed_${TIMESTAMP}.ipynb"

echo "[$(date '+%F %T')] START mosaic_run2-2.ipynb full rerun (corrected calibration)"
echo "Job ID: ${JOB_ID}, Host: $(hostname), Cores: ${NSLOTS:-8}"

# Notebook uses paths relative to its own directory (./mosaicfield_outputs, ./MOSAICField/...)
cd "${ALIGN_DIR}"

"${NBCONVERT}" \
    --to notebook \
    --execute \
    --output "${OUTPUT_NB}" \
    --ExecutePreprocessor.timeout=21600 \
    --ExecutePreprocessor.kernel_name=python3 \
    mosaic_run2-2.ipynb

echo "[$(date '+%F %T')] DONE — executed notebook saved to ${ALIGN_DIR}/${OUTPUT_NB}"
