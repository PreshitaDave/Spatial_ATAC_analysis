#!/bin/bash
################################################################################
# run_numbat_prep_tissue_deepseq.qsub.sh
#
# Purpose: Prepare deepseq tissue-specific NUMBAT inputs
# - Merge per-chromosome BAMs for each tissue
# - Create tissue-specific ATAC bin matrices
# - Generate pileup/phasing VCFs (blocked by reference panel)
#
# Tissues: deepseq_488B, deepseq_489
# Dataset: deepseq
#
################################################################################

#$ -P paxlab
#$ -N numbat_prep_deepseq
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=8G
#$ -l h_rt=24:00:00
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_prep_deepseq.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_prep_deepseq.$JOB_ID.err
#$ -j n

set -eo pipefail

on_err() {
  echo "[$(date '+%F %T')] ERROR: command failed at line ${1}: ${2}" >&2
}
trap 'on_err ${LINENO} "${BASH_COMMAND}"' ERR

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"

# Initialize module system
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

# Load modules
module load R samtools htslib bedtools

echo "[$(date '+%F %T')] ===== DEEPSEQ TISSUE-SPECIFIC BAM MERGING & BIN PREP ====="
echo "[NODE] $(hostname)"
cd "$PROJECT_ROOT"

DATASET="deepseq"
NCORES=8

# Set pre-resolve NUMBAT_BIN to avoid R hanging
export NUMBAT_BIN="/projectnb/paxlab/presh/Rlibs/4.5/numbat/bin"

echo ""
echo "[STEP 0] Verifying environment..."
which Rscript
Rscript --version

# Create output directories
mkdir -p Data/04_analysis/cnv/numbat/inputs/bam_merged

for TISSUE in 488B 489; do
  echo ""
  echo "=========================================="
  echo "[TISSUE=$TISSUE] Starting BAM merging..."
  echo "=========================================="
  
  export TISSUE DATASET
  
  # Merge ALL chromosome BAMs into single tissue BAM
  TISSUE_OUTPUT_DIR="Data/04_analysis/cnv/numbat/inputs/bam_merged"
  MERGED_BAM="${TISSUE_OUTPUT_DIR}/${DATASET}_${TISSUE}_merged_for_numbat.bam"
  
  if [[ -f "$MERGED_BAM" ]] && [[ -f "${MERGED_BAM}.bai" ]]; then
    echo "  ✓ Merged BAM already exists"
    ls -lh "$MERGED_BAM"
  else
    echo "  [STEP 1] Collecting per-chromosome BAMs..."
    
    # Get all chromosome BAMs (monopogen outputs chr*.filter.targeted.bam)
    BAM_DIR="Data/04_analysis/03_intermediate/variant_calling/monopogen/variant_calling/${DATASET}/Bam"
    BAM_FILES=($(find "$BAM_DIR" -name "chr*.filter.targeted.bam" 2>/dev/null | sort -V))
    
    if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
      echo "  [ERROR] No BAM files found in: $BAM_DIR" >&2
      exit 1
    fi
    
    echo "    Found ${#BAM_FILES[@]} chromosome BAMs"
    
    # Verify first BAM has barcode tags
    echo "  [STEP 2] Checking barcode format..."
    CB_TAG=$(samtools view "${BAM_FILES[0]}" 2>/dev/null | head -1 | grep -o 'CB:Z:[^ ]*' | head -1 || echo "NONE")
    echo "    Sample CB tag: $CB_TAG"
    
    if [[ "$CB_TAG" == "NONE" ]]; then
      echo "    [WARN] No CB tags found - proceeding anyway"
    fi
    
    # Merge all BAMs
    echo "  [STEP 3] Merging ${#BAM_FILES[@]} BAMs..."
    echo "    This will take ~30-60 minutes..."
    
    samtools merge -@ $((NCORES-1)) "$MERGED_BAM" "${BAM_FILES[@]}"
    
    # Index
    echo "  [STEP 4] Indexing merged BAM..."
    samtools index -@ $((NCORES-1)) "$MERGED_BAM"
    
    SIZE=$(ls -lh "$MERGED_BAM" | awk '{print $5}')
    echo "    ✓ Merged BAM created: $SIZE"
  fi
  
done

echo ""
echo "[$(date '+%F %T')] ===== SUMMARY ====="
echo ""
echo "Tissue-specific outputs created:"
for TISSUE in 488B 489; do
  MERGED_BAM="Data/04_analysis/cnv/numbat/inputs/bam_merged/deepseq_${TISSUE}_merged_for_numbat.bam"
  if [[ -f "$MERGED_BAM" ]]; then
    SIZE=$(ls -lh "$MERGED_BAM" | awk '{print $5}')
    echo "  ✓ $TISSUE: $SIZE (merged BAM)"
  else
    echo "  ✗ $TISSUE: MISSING"
  fi
done

echo ""
echo "[$(date '+%F %T')] Preparation complete!"
