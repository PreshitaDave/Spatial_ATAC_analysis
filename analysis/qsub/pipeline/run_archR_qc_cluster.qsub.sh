#!/bin/bash
#$ -P paxlab
#$ -N archR_qc_cluster
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l h_rt=24:00:00
#$ -l mem_per_core=8G
#$ -j y
#$ -o analysis/qsub_logs/build_tissue/archR_qc_cluster_$JOB_ID.log

set -euo pipefail

echo "[$(date +'%F %T')] [SETUP] Job ID: $JOB_ID"
echo "[$(date +'%F %T')] [SETUP] Hostname: $(hostname)"
echo "[$(date +'%F %T')] [SETUP] Working directory: $PWD"
echo "[$(date +'%F %T')] [SETUP] Available cores: $NSLOTS"

# Load required modules
echo "[$(date +'%F %T')] [MODULE] Loading R module..."
module load R
echo "[$(date +'%F %T')] [MODULE] R version:"
Rscript --version

# Verify input files exist
echo "[$(date +'%F %T')] [CHECK] Verifying input files..."
declare -a OBJECTS=("deepseq_488B" "deepseq_489" "lowseq_488B" "lowseq_489")
for obj in "${OBJECTS[@]}"; do
  fragments="/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/fragments/${obj}/${obj}.fragments.sort.filtered.bed.gz"
  barcodes="/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/${obj}/${obj}.no_edge_effect.barcodes.tsv"
  
  if [[ ! -f "$fragments" ]]; then
    echo "[ERROR] Missing fragments: $fragments"
    exit 1
  fi
  if [[ ! -f "$barcodes" ]]; then
    echo "[ERROR] Missing barcodes: $barcodes"
    exit 1
  fi
  echo "[OK] $obj: fragments ($(du -h "$fragments" | cut -f1)) + barcodes ($(wc -l < "$barcodes") bcs)"
done

# Create output directories
echo "[$(date +'%F %T')] [SETUP] Creating output directories..."
mkdir -p /projectnb/paxlab/presh/projects/spatial_atac/Data/01_outputs/archR_objects
mkdir -p /projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/cnv_analysis
mkdir -p /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/build_tissue

# Run ArchR pipeline
echo "[$(date +'%F %T')] [RUN] Starting ArchR QC pipeline..."
export PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
export NSLOTS=$NSLOTS

Rscript /projectnb/paxlab/presh/projects/spatial_atac/analysis/src/pipeline/archr/0_create_archr_qc_cluster.R

EXIT_CODE=$?
echo "[$(date +'%F %T')] [DONE] Script exited with code: $EXIT_CODE"

if [[ $EXIT_CODE -eq 0 ]]; then
  echo "[$(date +'%F %T')] [SUCCESS] ArchR pipeline completed successfully!"
  echo "[$(date +'%F %T')] [OUTPUT] ArchR objects: /projectnb/paxlab/presh/projects/spatial_atac/Data/01_outputs/archR_objects/"
  echo "[$(date +'%F %T')] [OUTPUT] PDF reports: /projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/cnv_analysis/"
else
  echo "[$(date +'%F %T')] [FAILED] ArchR pipeline failed with exit code $EXIT_CODE"
fi

exit $EXIT_CODE
