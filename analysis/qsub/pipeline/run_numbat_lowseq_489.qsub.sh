#!/bin/bash
#$ -P paxlab
#$ -N numbat_lowseq_489
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 6
#$ -l h_rt=24:00:00
#$ -l mem_per_core=8G
#$ -j y
#$ -o analysis/qsub_logs/build_tissue/numbat_lowseq_489_$JOB_ID.log

set -eo pipefail

echo "[$(date +'%F %T')] [SETUP] Job ID: $JOB_ID"
echo "[$(date +'%F %T')] [SETUP] Hostname: $(hostname)"
echo "[$(date +'%F %T')] [SETUP] Working directory: $PWD"
echo "[$(date +'%F %T')] [SETUP] Available cores: $NSLOTS"

# Initialize module environment
echo "[$(date +'%F %T')] [MODULE] Initializing module environment..."
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

# Load required modules
echo "[$(date +'%F %T')] [MODULE] Loading R module..."
module load R
echo "[$(date +'%F %T')] [MODULE] R version:"
Rscript --version

# Verify input files
echo "[$(date +'%F %T')] [CHECK] Verifying input files..."
PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"

ARCHR_PATH="$PROJECT_ROOT/Data/01_outputs/archR_objects/lowseq_489"
if [[ ! -d "$ARCHR_PATH" ]]; then
  echo "[ERROR] ArchR directory not found: $ARCHR_PATH"
  exit 1
fi
echo "[OK] ArchR project found"

BARCODE_FILE="$PROJECT_ROOT/Data/01_inputs/barcodes/tissue_barcodes/lowseq_489/lowseq_489.no_edge_effect.barcodes.tsv"
if [[ ! -f "$BARCODE_FILE" ]]; then
  echo "[ERROR] Barcode file not found: $BARCODE_FILE"
  exit 1
fi
BARCODE_COUNT=$(wc -l < "$BARCODE_FILE")
echo "[OK] Barcode file found ($BARCODE_COUNT barcodes)"

# Create output directories
echo "[$(date +'%F %T')] [SETUP] Creating output directories..."
mkdir -p "$PROJECT_ROOT/Data/04_analysis/cnv/numbat/results/lowseq/run_489"
mkdir -p "$PROJECT_ROOT/analysis/qsub_logs/build_tissue"

# Run NUMBAT pipeline
echo "[$(date +'%F %T')] [RUN] Starting NUMBAT pipeline for lowseq_489..."
export PROJECT_ROOT=$PROJECT_ROOT
export NSLOTS=$NSLOTS

Rscript "$PROJECT_ROOT/analysis/src/pipeline/numbat/run_numbat_lowseq_489.R"

EXIT_CODE=$?
echo "[$(date +'%F %T')] [DONE] Script exited with code: $EXIT_CODE"

if [[ $EXIT_CODE -eq 0 ]]; then
  echo "[$(date +'%F %T')] [SUCCESS] NUMBAT pipeline completed successfully!"
  echo "[$(date +'%F %T')] [OUTPUT] Results: $PROJECT_ROOT/Data/04_analysis/cnv/numbat/results/lowseq/run_489/"
else
  echo "[$(date +'%F %T')] [FAILED] NUMBAT pipeline failed with exit code $EXIT_CODE"
fi

exit $EXIT_CODE
