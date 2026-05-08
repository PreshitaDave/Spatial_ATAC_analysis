#!/usr/bin/env bash
set -euo pipefail

# ----------------------------
# Defaults
# ----------------------------
CORES=1
MEM_PER_CORE=8G
RUNTIME=08:00:00
PORT=8898
PROJECT=paxlab
GPU=0
GPU_RESOURCE_KEY=gpus
NOTEBOOK_DIR=/projectnb/paxlab/presh
LOGFILE="/tmp/jupyter_scc_${USER}.log"

# If no conda env is provided, this binary is used directly
JUPYTER=/projectnb/paxlab/presh/env/scenicplus/bin/jupyter

# Optional conda controls
CONDA_ENV=
CONDA_BASE=
MODULE_TO_LOAD=

usage() {
  cat << 'USAGE'
Usage:
  bash start_jupyter_compute.sh [options]

Options:
  -c, --cores N                 CPU cores (default: 8)
  -m, --mem-per-core SIZE       Memory per core, e.g. 8G (default: 8G)
  -r, --runtime HH:MM:SS        Walltime (default: 08:00:00)
  -p, --port PORT               Jupyter port (default: 8898)
  -P, --project NAME            SCC project/account (default: paxlab)
  -g, --gpu N                   Number of GPUs (default: 0)
      --gpu-resource-key KEY    GPU resource key for qrsh (default: gpus)
  -d, --notebook-dir PATH       Notebook root dir (default: /projectnb/paxlab/presh)
  -j, --jupyter PATH            Jupyter executable if not using conda activation
  -l, --logfile PATH            Log file path

      --conda-env NAME_OR_PATH  Conda env name or full env path
      --conda-base PATH         Conda base path, e.g. /path/to/miniconda3
      --module NAME             Optional module to load on compute node (example: miniconda)

  -h, --help                    Show help

Examples:
  bash start_jupyter_compute.sh --runtime 12:00:00 --cores 16 --project paxlab
  bash start_jupyter_compute.sh --gpu 1 --project paxlab
  bash start_jupyter_compute.sh --conda-env scenicplus --module miniconda
  bash start_jupyter_compute.sh --conda-env /projectnb/paxlab/presh/env/scenicplus
  bash start_jupyter_compute.sh --conda-base /path/to/miniconda3 --conda-env scenicplus
USAGE
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

# ----------------------------
# Parse args
# ----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -c|--cores) CORES="${2:-}"; shift 2 ;;
    -m|--mem-per-core) MEM_PER_CORE="${2:-}"; shift 2 ;;
    -r|--runtime) RUNTIME="${2:-}"; shift 2 ;;
    -p|--port) PORT="${2:-}"; shift 2 ;;
    -P|--project) PROJECT="${2:-}"; shift 2 ;;
    -g|--gpu) GPU="${2:-}"; shift 2 ;;
    --gpu-resource-key) GPU_RESOURCE_KEY="${2:-}"; shift 2 ;;
    -d|--notebook-dir) NOTEBOOK_DIR="${2:-}"; shift 2 ;;
    -j|--jupyter) JUPYTER="${2:-}"; shift 2 ;;
    -l|--logfile) LOGFILE="${2:-}"; shift 2 ;;
    --conda-env) CONDA_ENV="${2:-}"; shift 2 ;;
    --conda-base) CONDA_BASE="${2:-}"; shift 2 ;;
    --module) MODULE_TO_LOAD="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

# ----------------------------
# Validate inputs
# ----------------------------
[[ "$CORES" =~ ^[0-9]+$ ]] || die "--cores must be an integer"
[[ "$PORT" =~ ^[0-9]+$ ]] || die "--port must be an integer"
[[ "$GPU" =~ ^[0-9]+$ ]] || die "--gpu must be an integer"
[[ "$RUNTIME" =~ ^[0-9]{1,2}:[0-9]{2}:[0-9]{2}$ ]] || die "--runtime must be HH:MM:SS"

if [[ -z "$CONDA_ENV" ]]; then
  [[ -x "$JUPYTER" ]] || die "Jupyter not executable at: $JUPYTER"
fi

# ----------------------------
# Print request summary
# ----------------------------
echo "Requesting interactive compute node..."
echo "  Project      : $PROJECT"
echo "  Cores        : $CORES"
echo "  Mem/Core     : $MEM_PER_CORE"
echo "  Runtime      : $RUNTIME"
echo "  GPU count    : $GPU"
echo "  GPU key      : $GPU_RESOURCE_KEY"
echo "  Port         : $PORT"
echo "  Notebook dir : $NOTEBOOK_DIR"
echo "  Conda env    : ${CONDA_ENV:-<none>}"
echo "  Conda base   : ${CONDA_BASE:-<auto>}"
echo "  Module load  : ${MODULE_TO_LOAD:-<none>}"
echo ""

GPU_QRSH_ARGS=()
if [[ "$GPU" -gt 0 ]]; then
  GPU_QRSH_ARGS=(-l "${GPU_RESOURCE_KEY}=${GPU}")
fi

qrsh \
  -P "$PROJECT" \
  -l "h_rt=${RUNTIME}" \
  -pe omp "$CORES" \
  -l "mem_per_core=${MEM_PER_CORE}" \
  "${GPU_QRSH_ARGS[@]}" \
  -now no << EOF
set -euo pipefail

if [[ -n "${MODULE_TO_LOAD}" ]]; then
  if command -v module >/dev/null 2>&1; then
    module load "${MODULE_TO_LOAD}"
  else
    echo "WARN: module command not available; skipping module load ${MODULE_TO_LOAD}"
  fi
fi

COMPUTE_NODE=\$(hostname)
echo ">>> Running on compute node: \$COMPUTE_NODE"

if command -v python3 >/dev/null 2>&1; then
  TOKEN=\$(python3 -c "import secrets; print(secrets.token_hex(16))")
else
  TOKEN=\$(openssl rand -hex 16)
fi

JUPYTER_CMD="${JUPYTER}"

if [[ -n "${CONDA_ENV}" ]]; then
  if [[ -d "${CONDA_ENV}" && -x "${CONDA_ENV}/bin/jupyter" ]]; then
    JUPYTER_CMD="${CONDA_ENV}/bin/jupyter"
    echo "Using conda env path directly: ${CONDA_ENV}"
  else
    if [[ -n "${CONDA_BASE}" ]]; then
      source "${CONDA_BASE}/etc/profile.d/conda.sh"
    elif command -v conda >/dev/null 2>&1; then
      eval "\$(conda shell.bash hook)"
    elif [[ -f "\$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
      source "\$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "\$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
      source "\$HOME/anaconda3/etc/profile.d/conda.sh"
    else
      echo "ERROR: Could not initialize conda. Use --conda-base or --module miniconda."
      exit 1
    fi

    conda activate "${CONDA_ENV}"
    JUPYTER_CMD="\$(command -v jupyter)"
  fi
fi

echo "Python: \$(command -v python || true)"
echo "Jupyter: \$JUPYTER_CMD"
echo ""

echo "=========================================================="
echo " Starting Jupyter on \$COMPUTE_NODE:${PORT}"
echo "=========================================================="
echo "VS Code connection URL:"
echo "http://\$COMPUTE_NODE:${PORT}/?token=\$TOKEN"
echo ""
echo "SSH tunnel from your local machine:"
echo "ssh -NL ${PORT}:\$COMPUTE_NODE:${PORT} ${USER}@scc1.bu.edu"
echo "Then connect VS Code to:"
echo "http://localhost:${PORT}/?token=\$TOKEN"
echo "=========================================================="
echo ""

# Choose the best available Jupyter entrypoint in this environment.
if command -v jupyter-notebook >/dev/null 2>&1; then
  BACKEND="notebook"
  jupyter-notebook \
    --no-browser \
    --port="${PORT}" \
    --ip=0.0.0.0 \
    --NotebookApp.token="\$TOKEN" \
    --NotebookApp.notebook_dir="${NOTEBOOK_DIR}" \
    2>&1 | tee "${LOGFILE}"
elif "\$JUPYTER_CMD" lab --help >/dev/null 2>&1; then
  BACKEND="lab"
  "\$JUPYTER_CMD" lab \
    --no-browser \
    --port="${PORT}" \
    --ip=0.0.0.0 \
    --ServerApp.token="\$TOKEN" \
    --ServerApp.root_dir="${NOTEBOOK_DIR}" \
    2>&1 | tee "${LOGFILE}"
elif "\$JUPYTER_CMD" server --help >/dev/null 2>&1; then
  BACKEND="server"
  "\$JUPYTER_CMD" server \
    --no-browser \
    --port="${PORT}" \
    --ip=0.0.0.0 \
    --ServerApp.token="\$TOKEN" \
    --ServerApp.root_dir="${NOTEBOOK_DIR}" \
    2>&1 | tee "${LOGFILE}"
else
  echo "ERROR: No usable Jupyter backend found (notebook/lab/server)." >&2
  echo "Install one of: notebook, jupyterlab, jupyter_server" >&2
  exit 1
fi
EOF