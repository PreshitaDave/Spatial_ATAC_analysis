#!/bin/bash
#
# run_numbat_multiome_lowseq_489.qsub.sh
# Purpose: NUMBAT CNV analysis for lowseq_489 using multiome script
# Uses: run_numbat_multiome.R with command-line arguments
#

#$ -N numbat_multiome_489
#$ -l h_rt=08:00:00
#$ -pe omp 8
#$ -P paxlab
#$ -l mem_per_core=8G
#$ -o analysis/qsub_logs/numbat_multiome_489_$JOB_ID.log

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

echo "════════════════════════════════════════════════════════════════════════════"
echo "NUMBAT multiome Analysis: lowseq_489 (Command-line argument mode)"
echo "════════════════════════════════════════════════════════════════════════════"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Job ID: $JOB_ID"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Hostname: $(hostname)"
echo ""

# Configuration
SAMPLE="lowseq_489"
COUNTMAT="Data/04_analysis/cnv/numbat/inputs/${SAMPLE}/atac_bin/${SAMPLE}_atac_bin.rds"
ALLELEDF="Data/04_analysis/cnv/numbat/inputs/${SAMPLE}/alleles/${SAMPLE}_atac_allele_counts.tsv.gz"
OUT_DIR="Data/04_analysis/cnv/numbat/results/${SAMPLE}/multiome_20260519/"
REF="Data/04_analysis/cnv/numbat/reference/lambdas_ATAC_bincnt.rds"
GTF="Data/04_analysis/cnv/numbat/reference/var220kb.rds"
PARL="Data/04_analysis/cnv/numbat/reference/par_numbatm.rds"
R_SCRIPT="analysis/src/cnv_calling/numbat/run_numbat_multiome.R"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] Configuration:"
echo "[$(date +'%Y-%m-%d %H:%M:%S')]   Sample: $SAMPLE"
echo "[$(date +'%Y-%m-%d %H:%M:%S')]   Countmat: $COUNTMAT"
echo "[$(date +'%Y-%m-%d %H:%M:%S')]   Alleles: $ALLELEDF"
echo "[$(date +'%Y-%m-%d %H:%M:%S')]   Output: $OUT_DIR"
echo "[$(date +'%Y-%m-%d %H:%M:%S')]   Reference: $REF"
echo "[$(date +'%Y-%m-%d %H:%M:%S')]"

# Verify inputs exist
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Pre-flight checks..."
for f in "$COUNTMAT" "$ALLELEDF" "$REF" "$GTF" "$PARL" "$R_SCRIPT"; do
  if [[ ! -f "$f" ]]; then
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ ERROR: File not found: $f"
    exit 1
  fi
done
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ All input files verified"
echo ""

# Create output directory
mkdir -p "$OUT_DIR"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ Output directory created: $OUT_DIR"
echo ""

# Run NUMBAT multiome analysis
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ──────────────────────────────────────────────────────"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Calling: Rscript run_numbat_multiome.R"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ──────────────────────────────────────────────────────"
echo ""

Rscript "$R_SCRIPT" \
  --countmat "$COUNTMAT" \
  --alleledf "$ALLELEDF" \
  --out_dir "$OUT_DIR" \
  --ref "$REF" \
  --gtf "$GTF" \
  --parL "$PARL"

EXIT_CODE=$?

echo ""
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ──────────────────────────────────────────────────────"
if [[ $EXIT_CODE -eq 0 ]]; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✓ Job completed successfully (exit code: 0)"
else
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] ✗ Job failed with exit code: $EXIT_CODE"
fi
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ──────────────────────────────────────────────────────"
echo ""

exit $EXIT_CODE
