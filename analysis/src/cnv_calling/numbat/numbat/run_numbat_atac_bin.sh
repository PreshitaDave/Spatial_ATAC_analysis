#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${SCRIPT_DIR}/numbat_common.sh"

DATASET="${DATASET:-lowseq}"
NCORES="${NCORES:-4}"
NUMBAT_BIN=$(resolve_numbat_bin_dir)
RUN_SCRIPT="${NUMBAT_BIN}/run_numbat_multiome.R"
POSTPROCESS_SCRIPT="${SCRIPT_DIR}/postprocess_numbat_results.R"

COUNTMAT="${COUNTMAT:-${NUMBAT_INPUT_DIR}/${DATASET}_atac_bin.rds}"
ALLELEDF="${ALLELEDF:-${NUMBAT_INPUT_DIR}/alleles/${DATASET}_atac_allele_counts.tsv.gz}"
REF="${REF:-${NUMBAT_REF_DIR}/lambdas_ATAC_bincnt.rds}"
GTF_BIN="${GTF_BIN:-${NUMBAT_REF_DIR}/var220kb.rds}"
OUT_DIR="${OUT_DIR:-${NUMBAT_RESULT_DIR}/${DATASET}/atac_only_run2}"
PLOT_DIR="${PLOT_DIR:-${OUT_DIR}/plots}"
ITERATION="${ITERATION:-2}"
PARAMS="${PARAMS:-${NUMBAT_RESULT_DIR}/${DATASET}/atac_only_run2/par_numbat.rds}"

require_file "${RUN_SCRIPT}"
require_file "${POSTPROCESS_SCRIPT}"
require_file "${COUNTMAT}"
require_file "${ALLELEDF}"
require_file "${REF}"
require_file "${GTF_BIN}"

mkdir -p "${OUT_DIR}"

log "Running NUMBAT multiome in ATAC-bin mode for dataset=${DATASET}"
run_cmd "Rscript '${RUN_SCRIPT}' --countmat '${COUNTMAT}' --alleledf '${ALLELEDF}' --out_dir '${OUT_DIR}' --ref '${REF}' --gtf '${GTF_BIN}' --parL '${PARAMS}'"

log "Generating NUMBAT summary plots for dataset=${DATASET}"
run_cmd "Rscript '${POSTPROCESS_SCRIPT}' --out_dir '${OUT_DIR}' --plot_dir '${PLOT_DIR}' --iteration '${ITERATION}' --gtf '${GTF_BIN}' --dataset '${DATASET}'"

log "NUMBAT ATAC-bin run done. Output: ${OUT_DIR}"
