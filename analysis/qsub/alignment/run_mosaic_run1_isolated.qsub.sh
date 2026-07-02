#!/bin/bash -l
#$ -P paxlab
#$ -N mosaic_run1
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=16G
#$ -l h_rt=06:00:00
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/run_mosaic_run1_isolated_$JOB_ID.log
#$ -m bea
#$ -M preshita@bu.edu

set -euo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
RUN1_DIR="${PROJECT_ROOT}/analysis/src/alignment/run1_workspace"
NBCONVERT="/projectnb/paxlab/presh/env/conda_env/mosaicfield_env/bin/jupyter-nbconvert"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_NB="mosaic_run1_isolated_executed_${TIMESTAMP}.ipynb"

echo "[$(date '+%F %T')] START mosaic_run1_isolated.ipynb (original free-affine pipeline, isolated)"
echo "Job ID: ${JOB_ID}, Host: $(hostname), Cores: ${NSLOTS:-8}"

# This is a COPY of the pristine mosaic_run1.ipynb, run from an isolated working
# directory (run1_workspace/) so its outputs (./mosaicfield_outputs, ./MOSAICField
# symlink) never collide with the shared analysis/src/alignment/mosaicfield_outputs/
# used by other (rigid-transform) work. atac_affine_aligned.h5ad/xenium_affine_aligned.h5ad
# were already regenerated here via step1_regenerate_run1_inputs.py (median cells/spot=5,
# confirming this is genuinely the original, uncorrected-scale pipeline).
cd "${RUN1_DIR}"

"${NBCONVERT}" \
    --to notebook \
    --execute \
    --output "${OUTPUT_NB}" \
    --ExecutePreprocessor.timeout=21600 \
    --ExecutePreprocessor.kernel_name=python3 \
    mosaic_run1_isolated.ipynb

echo "[$(date '+%F %T')] DONE — executed notebook saved to ${RUN1_DIR}/${OUTPUT_NB}"
