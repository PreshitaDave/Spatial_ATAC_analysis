#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=02:00:00
#$ -l mem_per_core=8G
#$ -pe omp 1
#$ -P paxlab
#$ -o analysis/qsub_logs/alleloscope/generate_final_plots_lowseq_489.$JOB_ID.log

set -eo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Initialize module system
# ─────────────────────────────────────────────────────────────────────────────
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

module load R

# Verify R is available
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ===== Generate Step 10 Final Plots (Heatmap) =====" >&2
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Verifying R setup..." >&2
which Rscript
Rscript --version

# ─────────────────────────────────────────────────────────────────────────────
# Run the final plots script
# ─────────────────────────────────────────────────────────────────────────────
cd /projectnb/paxlab/presh/projects/spatial_atac

echo "[$(date +'%Y-%m-%d %H:%M:%S')] Starting final plots generation..." >&2

Rscript analysis/src/cnv_calling/alleloscope/alleloscope/lowseq/tissue/489/generate_final_plots_lowseq_489.R

EXIT_CODE=$?

if [[ $EXIT_CODE -eq 0 ]]; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ Final plots generation completed successfully" >&2
  
  # Verify output plots were created
  output_dir="Data/04_analysis/cnv/alleloscope/lowseq_tissue_from_existing/489/output/plots"
  if [[ -f "$output_dir/step6_CNV_coverage_heatmap.png" ]]; then
    size=$(ls -lh "$output_dir/step6_CNV_coverage_heatmap.png" | awk '{print $5}')
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ Heatmap created: $output_dir/step6_CNV_coverage_heatmap.png (size: $size)" >&2
  fi
  
  if [[ -f "$output_dir/rds/cov_obj.rds" ]]; then
    size=$(ls -lh "$output_dir/rds/cov_obj.rds" | awk '{print $5}')
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ Coverage object saved: $output_dir/rds/cov_obj.rds (size: $size)" >&2
  fi
else
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ Final plots generation FAILED (exit code: $EXIT_CODE)" >&2
  exit $EXIT_CODE
fi

echo "[$(date +'%Y-%m-%d %H:%M:%S')] ===== Complete =====" >&2
