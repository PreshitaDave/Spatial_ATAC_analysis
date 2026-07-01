#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${SCRIPT_DIR}/numbat_common.sh"

NUMBAT_REPO_DIR="${NUMBAT_REPO_DIR:-${NUMBAT_DATA_DIR}/numbat_repo}"
INSTALL_R_PKG="${INSTALL_R_PKG:-1}"

if [[ ! -d "${NUMBAT_REPO_DIR}/.git" ]]; then
  run_cmd "git clone https://github.com/kharchenkolab/numbat.git '${NUMBAT_REPO_DIR}'"
else
  run_cmd "git -C '${NUMBAT_REPO_DIR}' pull --ff-only"
fi

if [[ "${INSTALL_R_PKG}" == "1" ]]; then
  run_cmd "R --vanilla --slave -e \"if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='https://cloud.r-project.org'); remotes::install_github('kharchenkolab/numbat')\""
fi

log "NUMBAT repo ready at ${NUMBAT_REPO_DIR}"
log "Set NUMBAT_REPO=${NUMBAT_REPO_DIR} to force script usage from repo copy"
