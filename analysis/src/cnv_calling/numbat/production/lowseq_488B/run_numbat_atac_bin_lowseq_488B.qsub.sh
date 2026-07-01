#!/bin/bash
#$ -N numbat_atac_488B
#$ -l h_rt=24:00:00
#$ -pe omp 6
#$ -P paxlab
#$ -l mem_per_core=8G
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_atac_488B_.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_atac_488B_.err

set -eo pipefail

# Initialize module system
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

module load R
which Rscript && Rscript --version

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
cd "$PROJECT_ROOT"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] ════════════════════════════════════════════════════════════"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] NUMBAT ATAC-bin mode: lowseq_488B"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ════════════════════════════════════════════════════════════"

# Use NUMBAT package's run_numbat_multiome.R directly
NUMBAT_REPO="$PROJECT_ROOT/Data/04_analysis/cnv/numbat/numbat_repo"
RUN_SCRIPT="$NUMBAT_REPO/inst/bin/run_numbat_multiome.R"

# Input files
COUNTMAT="$PROJECT_ROOT/Data/04_analysis/cnv/numbat/inputs/lowseq_488B/atac_bin/lowseq_488B_atac_bin.rds"
ALLELEDF="$PROJECT_ROOT/Data/04_analysis/cnv/numbat/inputs/lowseq_488B/alleles/lowseq_488B_atac_allele_counts.tsv.gz"
REF="$PROJECT_ROOT/Data/04_analysis/cnv/numbat/reference/lambdas_ATAC_bincnt.rds"
GTF="$PROJECT_ROOT/Data/04_analysis/cnv/numbat/reference/var220kb.rds"
PARL="$PROJECT_ROOT/Data/04_analysis/cnv/numbat/reference/par_numbatm_custom.rds"
OUT_DIR="$PROJECT_ROOT/Data/04_analysis/cnv/numbat/results/lowseq_488B/atac_only_run2"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] Configuration:"
echo "[$(date +'%Y-%m-%d %H:%M:%S')]   Countmat: $COUNTMAT"
echo "[$(date +'%Y-%m-%d %H:%M:%S')]   Alleles: $ALLELEDF"
echo "[$(date +'%Y-%m-%d %H:%M:%S')]   Output: $OUT_DIR"
echo "[$(date +'%Y-%m-%d %H:%M:%S')]"

# Verify inputs
for f in "$COUNTMAT" "$ALLELEDF" "$REF" "$GTF" "$PARL" "$RUN_SCRIPT"; do
  if [[ ! -f "$f" ]]; then
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ ERROR: File not found: $f"
    exit 1
  fi
done
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ All input files verified"

# Create output directory
mkdir -p "$OUT_DIR"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ Output directory created"
echo "[$(date +'%Y-%m-%d %H:%M:%S')]"

# Run NUMBAT using package script
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ──────────────────────────────────────────────────────"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Running: Rscript $RUN_SCRIPT"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ──────────────────────────────────────────────────────"
echo "[$(date +'%Y-%m-%d %H:%M:%S')]"

Rscript "$RUN_SCRIPT" \
  --countmat "$COUNTMAT" \
  --alleledf "$ALLELEDF" \
  --out_dir "$OUT_DIR" \
  --ref "$REF" \
  --gtf "$GTF" \
  --parL "$PARL"

EXIT_CODE=$?

echo "[$(date +'%Y-%m-%d %H:%M:%S')]"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ──────────────────────────────────────────────────────"
if [[ $EXIT_CODE -eq 0 ]]; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ NUMBAT completed successfully"
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] Output: $OUT_DIR"
  ls -lh "$OUT_DIR"/*.rds "$OUT_DIR"/*.png "$OUT_DIR"/*.tsv.gz 2>/dev/null | head -10
else
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ NUMBAT failed with exit code: $EXIT_CODE"
fi
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ──────────────────────────────────────────────────────"

exit $EXIT_CODE
