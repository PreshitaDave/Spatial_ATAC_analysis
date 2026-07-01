#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=04:00:00
#$ -l mem_per_core=8G
#$ -pe omp 1
#$ -P paxlab
#$ -o analysis/qsub_logs/alleloscope/regenerate_fragments_lowseq_489.$JOB_ID.log

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
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ===== Regenerate Alleloscope chr1000k_fragments.tsv =====" >&2
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Verifying R setup..." >&2
which Rscript
Rscript --version

# ─────────────────────────────────────────────────────────────────────────────
# Run the regeneration script
# ─────────────────────────────────────────────────────────────────────────────
cd /projectnb/paxlab/presh/projects/spatial_atac

echo "[$(date +'%Y-%m-%d %H:%M:%S')] Starting fragment regeneration..." >&2

Rscript regenerate_alleloscope_fragments_489.R

EXIT_CODE=$?

if [[ $EXIT_CODE -eq 0 ]]; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ Regeneration completed successfully" >&2
  
  # Verify output file was created
  output_file="Data/04_analysis/cnv/alleloscope/lowseq_tissue_from_existing/489/chr1000k_fragments.tsv"
  if [[ -f "$output_file" ]]; then
    file_size=$(du -h "$output_file" | cut -f1)
    non_zero=$(awk '{for(i=2;i<=NF;i++) sum+=$i} END {print sum}' "$output_file")
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] Output file: $output_file (size: $file_size)" >&2
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] Total fragment counts: $non_zero" >&2
  fi
else
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ Regeneration FAILED (exit code: $EXIT_CODE)" >&2
  exit $EXIT_CODE
fi

echo "[$(date +'%Y-%m-%d %H:%M:%S')] ===== Complete =====" >&2
