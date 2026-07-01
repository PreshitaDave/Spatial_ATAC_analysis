#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
# Use real path to avoid symlink issues with mkdir
NUMBAT_DATA_DIR="${PROJECT_ROOT}/Data/04_analysis/cnv/numbat"
NUMBAT_INPUT_DIR="${NUMBAT_DATA_DIR}/inputs"
NUMBAT_REF_DIR="${NUMBAT_DATA_DIR}/reference"
NUMBAT_RESULT_DIR="${NUMBAT_DATA_DIR}/results"
NUMBAT_LOG_DIR="${NUMBAT_DATA_DIR}/logs"

DRY_RUN="${DRY_RUN:-0}"
DEFAULT_CHROMS="${DEFAULT_CHROMS:-chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22}"

mkdir -p "${NUMBAT_INPUT_DIR}" "${NUMBAT_REF_DIR}" "${NUMBAT_RESULT_DIR}" "${NUMBAT_LOG_DIR}"

log() {
  echo "[$(date '+%F %T')] $*"
}

run_cmd() {
  if [[ "${DRY_RUN}" == "1" ]]; then
    log "DRY_RUN: $*"
  else
    log "RUN: $*"
    eval "$*"
  fi
}

require_file() {
  local f="$1"
  if [[ ! -f "${f}" ]]; then
    if [[ "${DRY_RUN}" == "1" ]]; then
      log "DRY_RUN WARNING: required file not found: ${f}"
      return
    fi
    log "ERROR: required file not found: ${f}"
    exit 1
  fi
}

resolve_numbat_bin_dir() {
  if [[ -n "${NUMBAT_BIN_DIR:-}" && -d "${NUMBAT_BIN_DIR}" ]]; then
    echo "${NUMBAT_BIN_DIR}"
    return
  fi

  if [[ -n "${NUMBAT_REPO:-}" && -d "${NUMBAT_REPO}/inst/bin" ]]; then
    echo "${NUMBAT_REPO}/inst/bin"
    return
  fi

  local from_pkg
  from_pkg=$(R --vanilla --slave -e "cat(system.file('bin', package='numbat'))")
  if [[ -n "${from_pkg}" && -d "${from_pkg}" ]]; then
    echo "${from_pkg}"
    return
  fi

  if [[ "${DRY_RUN}" == "1" ]]; then
    log "DRY_RUN WARNING: could not resolve numbat bin scripts. Set NUMBAT_REPO or NUMBAT_BIN_DIR."
    echo "__NUMBAT_BIN_MISSING__"
    return
  fi
  log "ERROR: could not resolve numbat bin scripts. Set NUMBAT_REPO or NUMBAT_BIN_DIR."
  exit 1
}

split_csv() {
  local csv="$1"
  local old_ifs="${IFS}"
  IFS=',' read -r -a _items <<< "${csv}"
  IFS="${old_ifs}"
  printf '%s\n' "${_items[@]}"
}

join_by() {
  local delimiter="$1"
  shift
  local first="${1:-}"
  shift || true
  printf '%s' "${first}"
  local item
  for item in "$@"; do
    printf '%s%s' "${delimiter}" "${item}"
  done
}

collect_existing_files() {
  local pattern="$1"
  shift
  local -a found=()
  local chrom
  while IFS= read -r chrom; do
    [[ -z "${chrom}" ]] && continue
    local candidate
    printf -v candidate "${pattern}" "${chrom}"
    if [[ -f "${candidate}" ]]; then
      found+=("${candidate}")
    fi
  done < <(split_csv "${CHROMS:-${DEFAULT_CHROMS}}")

  printf '%s\n' "${found[@]}"
}

require_nonempty_list() {
  local label="$1"
  shift
  if [[ "$#" -eq 0 ]]; then
    log "ERROR: no files resolved for ${label}"
    exit 1
  fi
}
