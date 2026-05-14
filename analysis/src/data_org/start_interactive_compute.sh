#!/usr/bin/env bash
set -euo pipefail

# Start a persistent terminal (tmux/screen) and request an interactive compute node.
#
# Usage:
#   bash analysis/src/data_org/start_interactive_compute.sh --demo
#   bash analysis/src/data_org/start_interactive_compute.sh --self-test
#   bash analysis/src/data_org/start_interactive_compute.sh --run
#
# Notes:
# - --demo prints exactly what would run (no allocation requested).
# - --self-test checks command availability.
# - --run performs the full workflow.

SESSION_NAME="spatial_atac_agent"
QRSH_CMD="qrsh -l h_rt=16:00:00 -pe omp 1 -P paxlab -l mem_per_core=8G"
MODE="demo"
INSIDE="0"

for arg in "$@"; do
  case "$arg" in
    --demo)
      MODE="demo"
      ;;
    --self-test)
      MODE="self-test"
      ;;
    --run)
      MODE="run"
      ;;
    --inside)
      INSIDE="1"
      ;;
    *)
      echo "Unknown option: $arg" >&2
      exit 2
      ;;
  esac
done

log() {
  echo "[$(date +'%F %T')] $*"
}

in_persistent_shell() {
  [[ -n "${TMUX:-}" || -n "${STY:-}" ]]
}

is_login_host() {
  local h
  h="$(hostname 2>/dev/null || true)"
  [[ "$h" == scc1* || "$h" == *login* ]]
}

print_plan() {
  cat <<EOF
Plan:
1. Start persistent shell:
   tmux new -As ${SESSION_NAME}
   (fallback) screen -S ${SESSION_NAME}
2. Request interactive compute node:
   ${QRSH_CMD}
3. Validate node:
   hostname
EOF
}

self_test() {
  local ok=0
  log "Self-test: checking required commands"

  if command -v tmux >/dev/null 2>&1; then
    echo "PASS: tmux found"
  else
    echo "WARN: tmux not found"
  fi

  if command -v screen >/dev/null 2>&1; then
    echo "PASS: screen found"
  else
    echo "WARN: screen not found"
  fi

  if command -v qrsh >/dev/null 2>&1; then
    echo "PASS: qrsh found"
  else
    echo "FAIL: qrsh not found"
    ok=1
  fi

  echo "Host: $(hostname 2>/dev/null || echo unknown)"
  if is_login_host; then
    echo "Info: currently on login node"
  else
    echo "Info: currently on compute/non-login node"
  fi

  return "$ok"
}

start_persistent_then_run() {
  local script_path
  script_path="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"

  if command -v tmux >/dev/null 2>&1; then
    log "Starting/attaching tmux session: ${SESSION_NAME}"
    exec tmux new -As "${SESSION_NAME}" "bash ${script_path} --run --inside"
  fi

  if command -v screen >/dev/null 2>&1; then
    log "Starting screen session: ${SESSION_NAME}"
    exec screen -S "${SESSION_NAME}" bash -lc "bash ${script_path} --run --inside"
  fi

  echo "ERROR: neither tmux nor screen is available." >&2
  exit 1
}

run_workflow() {
  if [[ "$INSIDE" != "1" ]] && ! in_persistent_shell; then
    start_persistent_then_run
  fi

  log "Persistent session confirmed"
  log "Current host: $(hostname 2>/dev/null || echo unknown)"

  if is_login_host; then
    log "Requesting interactive compute node via qrsh"
    exec ${QRSH_CMD}
  else
    log "Already on compute/non-login host; no qrsh request needed"
    hostname
  fi
}

case "$MODE" in
  demo)
    print_plan
    ;;
  self-test)
    self_test
    ;;
  run)
    run_workflow
    ;;
  *)
    echo "Invalid mode: $MODE" >&2
    exit 2
    ;;
esac
