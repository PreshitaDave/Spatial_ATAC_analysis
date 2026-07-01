#!/bin/bash
################################################################################
# create_numbat_panel.sh
#
# Convert SNP VCF to BCF panel reference files for NUMBAT Eagle phasing
#
################################################################################

set -eo pipefail

# Load modules
module load samtools htslib 2>/dev/null || true

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
SNP_VCF="$PROJECT_ROOT/Data/hg38_resources/numbat/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
PANEL_DIR="$PROJECT_ROOT/Data/04_analysis/cnv/numbat/inputs/alleles/panel_ref"

echo "[BUILDING NUMBAT REFERENCE PANEL]"
echo "Input VCF: $SNP_VCF"
echo "Output panel: $PANEL_DIR"
echo ""

mkdir -p "$PANEL_DIR"

# Convert VCF to BCF and split by chromosome
for chr in {1..22}; do
  echo "[CHR $chr] Extracting and converting to BCF..."
  bcftools view -r chr$chr "$SNP_VCF" -O b --threads 4 -o "$PANEL_DIR/chr${chr}.genotypes.bcf"
  
  echo "[CHR $chr] Indexing BCF..."
  bcftools index -f "$PANEL_DIR/chr${chr}.genotypes.bcf" --threads 4
  
  SIZE=$(ls -lh "$PANEL_DIR/chr${chr}.genotypes.bcf" | awk '{print $5}')
  echo "[CHR $chr] ✓ Complete ($SIZE)"
done

echo ""
echo "[SUCCESS] Panel reference created with $(ls $PANEL_DIR/*.bcf | wc -l) files!"
ls -lh "$PANEL_DIR" | grep bcf | head -5
