#!/bin/bash
#$ -N numbat_prep_tissue_spec
#$ -l h_rt=24:00:00
#$ -l mem_per_core=8G
#$ -pe omp 8
#$ -P paxlab
#$ -e analysis/qsub_logs/numbat_prep_tissue_spec.$JOB_ID.err
#$ -o analysis/qsub_logs/numbat_prep_tissue_spec.$JOB_ID.out

set -euo pipefail

on_err() {
  echo "[$(date '+%F %T')] ERROR: command failed at line ${1}: ${2}" >&2
}
trap 'on_err ${LINENO} "${BASH_COMMAND}"' ERR

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
SCRIPT_DIR="analysis/src/cnv_calling/numbat/numbat"
PREP_SCRIPT="${PROJECT_ROOT}/${SCRIPT_DIR}/prepare_numbat_atac_inputs.sh"

# Export NUMBAT_REPO so resolve_numbat_bin_dir works in child scripts
export NUMBAT_REPO="${PROJECT_ROOT}/Data/04_analysis/cnv/numbat/numbat_repo"

# Initialize module system
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

# Verify modules and Rscript
echo "[$(date '+%F %T')] [STEP 0] Verifying compute environment..."
hostname
module load R samtools htslib bedtools
which Rscript && Rscript --version
echo "✓ Compute environment ready"

cd "$PROJECT_ROOT"

# Verify prep script exists
if [[ ! -f "$PREP_SCRIPT" ]]; then
  echo "[ERROR] Prep script not found: $PREP_SCRIPT" >&2
  exit 1
fi

# Pre-set NUMBAT_BIN to avoid R hanging issues during function resolution
export NUMBAT_BIN="/projectnb/paxlab/presh/Rlibs/4.5/numbat/bin"

echo ""
echo "=============================================="
echo "NUMBAT Tissue-Specific Pileup/Phasing Re-run"
echo "=============================================="
echo ""

# Run for each tissue
for TISSUE in 488B 489; do
  echo ""
  echo "[$(date '+%F %T')] [TISSUE=$TISSUE] Starting pileup/phasing preparation..."
  
  export DATASET="lowseq"
  export TISSUE="${TISSUE}"
  export NCORES=8
  
  if DATASET="$DATASET" TISSUE="$TISSUE" NCORES=4 bash "$PREP_SCRIPT" 2>&1 | tee "analysis/qsub_logs/numbat_prep_lowseq_${TISSUE}.log"; then
    echo "[$(date '+%F %T')] ✓ Tissue $TISSUE completed successfully"
  else
    echo "[$(date '+%F %T')] ✗ Tissue $TISSUE failed; continuing to next tissue..."
  fi
done

echo ""
echo "=============================================="
echo "Re-run Summary"
echo "=============================================="
echo ""
echo "[STEP 1] Checking created outputs..."

# Verify tissue-specific outputs were created
for TISSUE in 488B 489; do
  ALLELE_FILE="Data/04_analysis/cnv/numbat/inputs/alleles/lowseq_${TISSUE}_atac_allele_counts.tsv.gz"
  PILEUP_VCF="Data/04_analysis/cnv/numbat/inputs/alleles/pileup/lowseq_${TISSUE}_atac/cellSNP.base.vcf"
  
  if [[ -f "$ALLELE_FILE" ]]; then
    SIZE=$(du -h "$ALLELE_FILE" | cut -f1)
    echo "✓ lowseq_${TISSUE} allele counts: $SIZE"
  else
    echo "✗ lowseq_${TISSUE} allele counts NOT FOUND"
  fi
  
  if [[ -f "$PILEUP_VCF" ]]; then
    VARIANTS=$(grep -vc '^#' "$PILEUP_VCF" || true)
    echo "✓ lowseq_${TISSUE} pileup variants: $VARIANTS"
  else
    echo "✗ lowseq_${TISSUE} pileup VCF NOT FOUND"
  fi
done

echo ""
echo "[STEP 2] Tissue segregation verification..."
TISSUES_WITH_ALLELES=$(ls -1 Data/04_analysis/cnv/numbat/inputs/alleles/*_atac_allele_counts.tsv.gz 2>/dev/null | wc -l)
echo "Total tissue-specific allele files: $TISSUES_WITH_ALLELES (expected: 2)"

PILEUP_DIRS=$(find Data/04_analysis/cnv/numbat/inputs/alleles/pileup -maxdepth 1 -type d -name "*_atac" 2>/dev/null | wc -l)
echo "Total tissue-specific pileup dirs: $PILEUP_DIRS (expected: 2)"

if [[ "$TISSUES_WITH_ALLELES" -eq 2 && "$PILEUP_DIRS" -eq 2 ]]; then
  echo "✓ All tissue-specific files created successfully!"
else
  echo "✗ Tissue segregation incomplete or failed"
fi

echo ""
echo "[$(date '+%F %T')] Preparation complete!"
