#!/bin/bash
################################################################################
# setup_numbat_references.sh
#
# Purpose: Download and organize NUMBAT reference files needed for pileup/phasing
#
# Reference files needed:
# 1. genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf - SNP VCF from 1000G for pileup
# 2. 1000G_hg38/ - BCF files for phasing with Eagle
# 3. genetic_map_hg38_withX.txt.gz - (already available)
#
# Output: Data/02_references/genome/hg38_resources/numbat/
#
# Author:  Spatial ATAC Pipeline
# Date:    May 17, 2026
#
################################################################################

set -eo pipefail

PROJ_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
REF_DIR="${PROJ_ROOT}/Data/02_references/genome/hg38_resources/numbat"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] Setting up NUMBAT reference files..."
echo "[INFO] Reference directory: ${REF_DIR}"

mkdir -p "${REF_DIR}"

# Create symlink for easier reference from hg38_resources
if [[ ! -L "${PROJ_ROOT}/Data/hg38_resources" ]]; then
  ln -sf Data/02_references/genome/hg38_resources "${PROJ_ROOT}/Data/hg38_resources"
  echo "[INFO] Created symlink: Data/hg38_resources -> Data/02_references/genome/hg38_resources"
fi

# Create symlink to existing genetic map
if [[ ! -f "${REF_DIR}/genetic_map_hg38_withX.txt.gz" ]]; then
  ln -sf /projectnb/paxlab/presh/software/external/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
    "${REF_DIR}/"
  echo "[INFO] Linked genetic map"
fi

# Download 1000G SNP VCF for pileup (smaller version with AF > 5e-2)
echo "[INFO] Checking SNP VCF..."
if [[ ! -f "${REF_DIR}/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf" ]]; then
  echo "[INFO] Downloading 1000G SNP VCF..."
  cd "${REF_DIR}"
  
  # Try official 1000G source
  if wget -q --timeout=30 "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/ALL.wgs.nhomogeneous_calling_region.20170711.sites.vcf.gz" -O temp.vcf.gz 2>/dev/null; then
    gunzip -c temp.vcf.gz > genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf
    rm -f temp.vcf.gz
    echo "[✓] Downloaded SNP VCF successfully"
  else
    echo "[WARN] Could not download SNP VCF from 1000G - manual download may be needed"
    echo "[INFO] Files available at: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase3_reference_panel_GRCh38/"
  fi
fi

# Create 1000G phasing panel directory
echo "[INFO] Setting up 1000G phasing panel..."
if [[ ! -d "${REF_DIR}/1000G_hg38" ]]; then
  mkdir -p "${REF_DIR}/1000G_hg38"
  echo "[INFO] Created 1000G_hg38 directory for BCF files"
  echo "[INFO] NOTE: Individual chr*.bcf files need to be downloaded separately"
  echo "[INFO] Source: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
else
  echo "[✓] 1000G_hg38 directory exists"
fi

# List what exists and what's missing
echo ""
echo "===================== Reference File Status ====================="
ls -lh "${REF_DIR}/" | grep -v "^total" | awk '{print $9, "(" $5 ")"}'
echo ""
if [[ ! -f "${REF_DIR}/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf" ]]; then
  echo "[WARN] ✗ SNP VCF MISSING - download required"
else
  echo "[✓] SNP VCF available"
fi

if [[ -d "${REF_DIR}/1000G_hg38" ]]; then
  BCF_COUNT=$(find "${REF_DIR}/1000G_hg38" -name "*.bcf" 2>/dev/null | wc -l)
  echo "[INFO] BCF files in 1000G_hg38: $BCF_COUNT/22"
  if [[ $BCF_COUNT -lt 22 ]]; then
    echo "[WARN] ✗ Some BCF files missing - download all 22 chromosomes from:"
    echo "      https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
  fi
fi

echo ""
echo "[INFO] Setup complete. Check messages above for any issues."
