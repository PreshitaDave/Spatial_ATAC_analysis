#!/bin/bash
#$ -N archr_create_no_doublet
#$ -l h_rt=8:00:00
#$ -l mem_per_core=8G
#$ -pe omp 8
#$ -P paxlab
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/

set -euo pipefail

# Initialize module system
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

# Load R module
module load R 2>/dev/null || module load gcc/9.2.0 R/4.2.0 2>/dev/null || true

# Verify R is available
which Rscript || { echo "[ERROR] Rscript not found"; exit 1; }
Rscript --version

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting ArchR project creation (NO DOUBLET REMOVAL)..."

cd /projectnb/paxlab/presh/projects/spatial_atac

# Run the ArchR project creation script
# This will:
# 1. Create ArchR projects from Arrow files for all tissues
# 2. Filter to edge-effect-filtered barcodes
# 3. Apply basic QC (TSS >= 3, nFrags >= 1000) - NO doublet removal
# 4. Add LSI, clustering, and UMAP
# 5. Generate 4-panel QC plots (TSS, nFrags, scatter, UMAP clustering)
Rscript analysis/src/build_tissue/0_create_archr_no_doublet.R

echo "[$(date '+%Y-%m-%d %H:%M:%S')] ArchR project creation complete!"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output files:"
echo "  Projects created in: Data/01_outputs/archR_objects/{tissue}/{tissue}_archR_project_final/"
echo "  QC Plots (4 panels: TSS, nFrags, scatter, UMAP clustering):"
ls -lh analysis/plots/archr_obj/archR_qc_*.pdf 2>/dev/null | tail -10 || echo "  (PDF plots generated)"
echo "  Final barcodes saved in: Data/01_outputs/archR_objects/{tissue}/{tissue}_final_cell_barcodes.txt"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] ArchR projects recreated WITHOUT doublet removal!"
