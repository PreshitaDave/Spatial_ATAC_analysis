#!/bin/bash
set -euo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
LOG_DIR="${PROJECT_ROOT}/analysis/qsub_logs"
SCRIPT="${PROJECT_ROOT}/analysis/qsub/pipeline/alleloscope/gatk/run_gatk_hc_chr.qsub.sh"
DATASETS="${DATASETS:-deepseq lowseq}"
CHR_START="${CHR_START:-1}"
CHR_END="${CHR_END:-22}"
PE_SLOTS="${PE_SLOTS:-4}"
MEM_PER_CORE="${MEM_PER_CORE:-8G}"
H_RT="${H_RT:-24:00:00}"

mkdir -p "${LOG_DIR}"

for DATASET in ${DATASETS}; do
  for CHR_NUM in $(seq "${CHR_START}" "${CHR_END}"); do
    CHR="chr${CHR_NUM}"
    if [[ "${DATASET}" == "deepseq" ]]; then
      DS_TAG="d"
    else
      DS_TAG="l"
    fi
    JOB_NAME="gatk_${DS_TAG}_c${CHR_NUM}"
    echo "Submitting ${JOB_NAME}"
    qsub \
      -P paxlab \
      -pe omp "${PE_SLOTS}" \
      -l mem_per_core="${MEM_PER_CORE}" \
      -l h_rt="${H_RT}" \
      -N "${JOB_NAME}" \
      -o "${LOG_DIR}" \
      -e "${LOG_DIR}" \
      -v "DATASET=${DATASET},CHR=${CHR}" \
      "${SCRIPT}"
  done
done
