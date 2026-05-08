#!/bin/bash
# merge_all_chr_gl.sh -- merge all chrN.gl.vcf.gz (1-22) into one combined file

set -euo pipefail

module load htslib/1.16
module load bcftools/1.16

# Change to working directory with VCF files
cd /projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/lowseq/germline

# Collect all gl VCF files for chromosomes 1-22
echo "Collecting gl VCF files for chromosomes 1-22..."
VCF_FILES=""
for CHR in {1..22}; do
    VCF_FILE="chr${CHR}.gl.vcf.gz"
    if [[ -f "$VCF_FILE" ]]; then
        echo "  Found: $VCF_FILE"
        VCF_FILES="$VCF_FILES $VCF_FILE"
    else
        echo "  WARNING: $VCF_FILE not found!"
    fi
done

if [[ -z "$VCF_FILES" ]]; then
    echo "ERROR: No gl VCF files found!"
    exit 1
fi

# Merge all VCF files
OUTPUT="combined.gl.vcf.gz"
echo ""
echo "Merging all chromosome gl VCFs into: $OUTPUT"

bcftools concat -Oz $VCF_FILES -o "$OUTPUT"

# Index the merged VCF
echo "Indexing merged VCF..."
tabix -p vcf "$OUTPUT"

echo "Done! Created:"
echo "  - $OUTPUT"
echo "  - ${OUTPUT}.tbi"
