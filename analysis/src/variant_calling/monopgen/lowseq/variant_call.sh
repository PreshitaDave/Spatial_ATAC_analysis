#!/bin/bash -l
#$ -P paxlab
#$ -N monopogen_germline_deepseq
#$ -pe omp 4
#$ -l h_rt=24:00:00
#$ -l mem_per_core=8G
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/germline_deepseq
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/germline_deepseq_error

set -euo pipefail

# --- config ---
MONOPOGEN_PATH="/projectnb/paxlab/presh/software/Monopogen"

OUT_DIR="/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/lowseq"
REF_FA="/projectnb/paxlab/presh/projects/spatial_atac/Data/hg38_resources/Homo_sapiens_assembly38.fasta"

REGION_LST="${MONOPOGEN_PATH}/resource/GRCh38.region.lst"
# Use your canonical region list if preferred
# REGION_LST="/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/GRCh38.chr1_22_X_Y.region.lst"
WORK_DIR="/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/lowseq"
REF_FA="/projectnb/paxlab/presh/projects/spatial_atac/Data/hg38_resources/Homo_sapiens_assembly38.fasta"
IMPUTATION_PANEL="/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/phased_hg38_1000G_ref/"
STEP="all"
THREADS="${NSLOTS:-16}"

# --- env ---
module load python3
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${MONOPOGEN_PATH}/apps"

mkdir -p "${OUT_DIR}"
cd "${WORK_DIR}"

# --- run ---
python "${MONOPOGEN_PATH}/src/Monopogen.py" germline \
  -r "${REGION_LST}" \
  -s "${STEP}" \
  -o "${OUT_DIR}" \
  -g "${REF_FA}" \
  -p "${IMPUTATION_PANEL}" \
  -a "${MONOPOGEN_PATH}/apps" \
  -t "${THREADS}"  --norun TRUE

echo "Done: ${OUT_DIR}"


