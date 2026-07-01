#!/bin/bash
################################################################################
# create_snp_vcf_from_variants.sh
#
# Purpose: Merge monopogen VCFs and create SNP-only VCF for NUMBAT pileup
#          Uses tissue-specific variants from monopogen variant calling
#
# Input:   Data/04_analysis/03_intermediate/variant_calling/monopogen/variant_calling/{dataset}/gatk_hc/chr*.hc.vcf.gz
# Output:  Data/02_references/genome/hg38_resources/numbat/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf
#
# Author:  Spatial ATAC Pipeline
# Date:    May 17, 2026
#
################################################################################

set -eo pipefail

echo "[$(date +'%Y-%m-%d %H:%M:%S')] Creating SNP VCF from monopogen variants..."

PROJECT_ROOT="${PROJECT_ROOT:-.}"
DATASET="${1:-lowseq}"

# Paths
VARIANT_DIR="${PROJECT_ROOT}/Data/04_analysis/03_intermediate/variant_calling/monopogen/variant_calling/${DATASET}/gatk_hc"
REF_DIR="${PROJECT_ROOT}/Data/02_references/genome/hg38_resources/numbat"
OUTPUT_VCF="${REF_DIR}/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf"

# Verify input
if [[ ! -d "$VARIANT_DIR" ]]; then
  echo "[ERROR] Variant directory not found: $VARIANT_DIR" >&2
  exit 1
fi

echo "[INFO] Source dataset: $DATASET"
echo "[INFO] Variant input: $VARIANT_DIR"
echo "[INFO] Output: $OUTPUT_VCF"

# Check for VCF files
VCFCOUNT=$(find "$VARIANT_DIR" -name "chr*.hc.vcf.gz" -type f | wc -l)
if [[ $VCFCOUNT -eq 0 ]]; then
  echo "[ERROR] No chr*.hc.vcf.gz files found in $VARIANT_DIR" >&2
  exit 1
fi

echo "[INFO] Found $VCFCOUNT chromosome VCF files"

# Create output directory
mkdir -p "$REF_DIR"

# Collect list of VCF files in chromosome order
VCFLIST=()
for CHR in {1..22} X; do
  VCF="${VARIANT_DIR}/chr${CHR}.hc.vcf.gz"
  if [[ -f "$VCF" ]]; then
    VCFLIST+=("$VCF")
  else
    echo "[WARN] Missing: chr${CHR}.hc.vcf.gz" >&2
  fi
done

echo "[INFO] Merging ${#VCFLIST[@]} VCF files..."

# Temporary file for merged VCF
TEMP_MERGED=$(mktemp)
trap "rm -f $TEMP_MERGED" EXIT

# Merge all VCFs, keeping only SNPs
echo "[INFO] Extracting SNPs from all chromosomes..."
{
  # Write header from first file
  zcat "${VCFLIST[0]}" | grep "^#" 
  
  # Concatenate and filter all files for SNPs only
  for VCF in "${VCFLIST[@]}"; do
    zcat "$VCF" | grep -v "^#" | awk '$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT](,[ACGT])*$/ {print}'
  done
} | sort -k1,1V -k2,2n > "$TEMP_MERGED"

TOTAL_VARIANTS=$(grep -vc "^#" "$TEMP_MERGED" || true)
echo "[INFO] Total SNPs after filtering: $TOTAL_VARIANTS"

# Create final VCF
if [[ $TOTAL_VARIANTS -gt 0 ]]; then
  mv "$TEMP_MERGED" "$OUTPUT_VCF"
  echo "[✓] SNP VCF created: $OUTPUT_VCF"
  echo "[INFO] VCF size: $(du -h "$OUTPUT_VCF" | cut -f1)"
else
  echo "[ERROR] No SNPs found after filtering!" >&2
  rm -f "$TEMP_MERGED"
  exit 1
fi

echo "[$(date +'%Y-%m-%d %H:%M:%S')] Done!"
