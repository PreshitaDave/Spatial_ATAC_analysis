#!/bin/bash
#$ -N test_build_tissue_edge
#$ -l h_rt=04:00:00
#$ -pe omp 1
#$ -P paxlab
#$ -l mem_per_core=8G
#$ -o analysis/qsub_logs/test_build_tissue_edge.log
#$ -j y

echo "[$(date '+%Y-%m-%d %H:%M:%S')] ========================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] JOB START: build_tissue_barcodes_edge_nfrags_plots"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Hostname: $(hostname)"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] PWD: $PWD"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] ========================================="

# Set working directory
cd /projectnb/paxlab/presh/projects/spatial_atac

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 1: Checking R availability..."
which Rscript
echo "[$(date '+%Y-%m-%d %H:%M:%S')] R version:"
Rscript --version

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 2: Checking input files..."
echo "Fragment files:"
ls -lh Data/01_inputs/fragments/*/
echo "Spatial file:"
ls -lh Data/tissue_positions_list.csv*

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 3: Running build_tissue_barcodes_edge_nfrags_plots.R with EDGE_AXIS=row..."
EDGE_AXIS=row Rscript analysis/src/build_tissue/build_tissue_barcodes_edge_nfrags_plots.R 2>&1

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 4: Checking outputs..."
echo "Barcode outputs:"
find Data/01_inputs/barcodes/tissue_barcodes -name "*.barcodes.tsv" -type f -exec ls -lh {} \;

echo "Plot outputs:"
find analysis/plots/variant_qc/edge_effect_nfrags -name "*.png" -type f -exec ls -lh {} \;

echo "[$(date '+%Y-%m-%d %H:%M:%S')] ========================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] JOB COMPLETED"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] ========================================="

