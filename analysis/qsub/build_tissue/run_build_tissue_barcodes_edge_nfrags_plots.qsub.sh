#!/bin/bash -l
#$ -P paxlab
#$ -N edge_nfrags
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 4
#$ -l h_rt=08:00:00
#$ -l mem_per_core=8G
#$ -j n
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/edge_nfrags.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/edge_nfrags.$JOB_ID.err

set -euo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
SCRIPT_PATH="analysis/src/build_tissue/build_tissue_barcodes_edge_nfrags_plots.R"

# Ensure environment modules are available in batch context.
if [[ -f /etc/profile.d/modules.sh ]]; then
	# shellcheck source=/dev/null
	source /etc/profile.d/modules.sh
fi
if ! command -v module >/dev/null 2>&1; then
	echo "[$(date '+%F %T')] [error] module command not found" >&2
	exit 127
fi

module load R
if ! command -v Rscript >/dev/null 2>&1; then
	echo "[$(date '+%F %T')] [error] Rscript still missing after module load R" >&2
	exit 127
fi

echo "[$(date '+%F %T')] [start] job_id=${JOB_ID:-NA} host=$(hostname)"
echo "[$(date '+%F %T')] [start] cwd=${PWD}"
echo "[$(date '+%F %T')] [start] script=${SCRIPT_PATH}"
echo "[$(date '+%F %T')] [start] NSLOTS=${NSLOTS:-4}"
echo "[$(date '+%F %T')] [start] R=$(command -v Rscript 2>/dev/null || echo MISSING)"

cd "${PROJECT_ROOT}"

echo "[$(date '+%F %T')] [step] running edge nFrags script"
Rscript "${SCRIPT_PATH}"

echo "[$(date '+%F %T')] [done] job_id=${JOB_ID:-NA} finished on $(hostname)"