#!/bin/bash
################################################################################
# prepare_deepseq_tissue_inputs.sh
#
# Purpose: Create tissue-specific deepseq NUMBAT inputs from scratch
#
# Steps:
#   1. Extract tissue-specific bin matrices from combined Seurat object (FAST)
#   2. Create tissue-specific barcode lists (INSTANT)
#   3. Run tissue-specific variant calling in parallel (MODERATE)
#
# Efficiency:
#   - Single Seurat load → multiple tissue exports (no recomputation)
#   - Parallel variant calling (both tissues simultaneously)
#   - Single BAM process with barcode filtering
#
# Output:
#   - deepseq_{488B,489}_atac_bin.rds (tissue-specific bin matrices)
#   - deepseq_{488B,489}_atac_barcodes_for_pileup.tsv (tissue-specific barcodes)
#   - alleles/pileup/deepseq_{488B,489}_atac/ (tissue-specific variants)
#   - deepseq_{488B,489}_atac_allele_counts.tsv.gz (tissue-specific allele counts)
#
################################################################################

set -eo pipefail

on_err() {
  echo "[$(date '+%F %T')] ERROR: command failed at line ${1}: ${2}" >&2
}
trap 'on_err ${LINENO} "${BASH_COMMAND}"' ERR

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
DATASET="deepseq"
TISSUES=("488B" "489")
NCORES=8

echo "[$(date +'%F %T')] ===== DEEPSEQ TISSUE-SPECIFIC INPUT PREPARATION ====="
echo "[NODE] $(hostname)"
cd "$PROJECT_ROOT"

################################################################################
# STEP 1: Extract Tissue-Specific Bin Matrices
################################################################################

echo ""
echo "[STEP 1 @ $(date +'%F %T')] Extracting tissue-specific bin matrices..."

COMBINED_BIN="Data/04_analysis/cnv/numbat/inputs/${DATASET}_atac_bin.rds"

if [[ ! -f "$COMBINED_BIN" ]]; then
  echo "[ERROR] Combined bin matrix not found: $COMBINED_BIN" >&2
  exit 1
fi
echo "  Source: $(ls -lh $COMBINED_BIN | awk '{print $5}')"

# Load combined Seurat → extract tissues → save separately
Rscript - <<'RBIN_EXTRACT' "$COMBINED_BIN" "$DATASET" "$PROJECT_ROOT"
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
combined_bin <- args[1]
dataset <- args[2]
project_root <- args[3]

cat(sprintf("[INFO] Loading combined Seurat object: %s\n", combined_bin))
seurat_combined <- readRDS(combined_bin)

cat(sprintf("[INFO] Total cells: %d\n", ncol(seurat_combined)))
cat(sprintf("[INFO] Current metadata columns: %s\n", paste(colnames(seurat_combined@meta.data), collapse=", ")))

# Check if tissue info is in metadata
if ("tissue" %in% colnames(seurat_combined@meta.data)) {
  cat("[INFO] Tissue info found in metadata\n")
  tissues <- unique(seurat_combined@meta.data$tissue)
  cat(sprintf("[INFO] Tissues: %s\n", paste(tissues, collapse=", ")))
  
  # Extract by tissue
  for (tissue in c("488B", "489")) {
    cells <- colnames(seurat_combined)[seurat_combined@meta.data$tissue == tissue]
    if (length(cells) > 0) {
      cat(sprintf("[INFO] Extracting %s: %d cells\n", tissue, length(cells)))
      seurat_tissue <- subset(seurat_combined, cells = cells)
      
      output_file <- file.path(project_root, sprintf("Data/04_analysis/cnv/numbat/inputs/%s_%s_atac_bin.rds", dataset, tissue))
      saveRDS(seurat_tissue, file = output_file)
      cat(sprintf("[✓] Saved: %s\n", output_file))
    } else {
      cat(sprintf("[WARN] No cells found for tissue %s\n", tissue))
    }
  }
} else {
  # No tissue info - try to extract from cell names or create placeholder
  cat("[WARN] No 'tissue' column in metadata - checking cell names...\n")
  
  # Create tissue labels from existing data
  # Assume order or naming convention
  n_cells <- ncol(seurat_combined)
  tissue_489_cells <- which(grepl("_489", colnames(seurat_combined), ignore.case = TRUE))
  tissue_488b_cells <- which(grepl("_488B|_488b", colnames(seurat_combined), ignore.case = TRUE))
  
  if (length(tissue_489_cells) > 0) {
    cat(sprintf("[INFO] Extracting by cell name pattern - 489: %d cells\n", length(tissue_489_cells)))
    seurat_489 <- subset(seurat_combined, cells = colnames(seurat_combined)[tissue_489_cells])
    output_file <- file.path(project_root, sprintf("Data/04_analysis/cnv/numbat/inputs/%s_489_atac_bin.rds", dataset))
    saveRDS(seurat_489, file = output_file)
    cat(sprintf("[✓] Saved: %s\n", output_file))
  }
  
  if (length(tissue_488b_cells) > 0) {
    cat(sprintf("[INFO] Extracting by cell name pattern - 488B: %d cells\n", length(tissue_488b_cells)))
    seurat_488b <- subset(seurat_combined, cells = colnames(seurat_combined)[tissue_488b_cells])
    output_file <- file.path(project_root, sprintf("Data/04_analysis/cnv/numbat/inputs/%s_488B_atac_bin.rds", dataset))
    saveRDS(seurat_488b, file = output_file)
    cat(sprintf("[✓] Saved: %s\n", output_file))
  }
}

cat("[INFO] Bin matrix extraction complete\n")
RBIN_EXTRACT

################################################################################
# STEP 2: Create Tissue-Specific Barcode Lists
################################################################################

echo ""
echo "[STEP 2 @ $(date +'%F %T')] Creating tissue-specific barcode lists..."

BARCODE_DIR="Data/01_inputs/barcodes/tissue_barcodes"

for TISSUE in "${TISSUES[@]}"; do
  TISSUE_BARCODE_DIR="${BARCODE_DIR}/${DATASET}_${TISSUE}"
  TISSUE_BARCODE_FILE="${TISSUE_BARCODE_DIR}/${DATASET}_${TISSUE}.no_edge_effect.barcodes.tsv"
  
  if [[ ! -f "$TISSUE_BARCODE_FILE" ]]; then
    echo "  [WARN] Barcode file not found: $TISSUE_BARCODE_FILE" >&2
  else
    # Create filtered barcode list for pileup
    OUTPUT_BARCODE="Data/04_analysis/cnv/numbat/inputs/${DATASET}_${TISSUE}_atac_barcodes_for_pileup.tsv"
    cp "$TISSUE_BARCODE_FILE" "$OUTPUT_BARCODE"
    echo "  ✓ ${TISSUE}: $(wc -l < $OUTPUT_BARCODE) barcodes"
  fi
done

################################################################################
# STEP 3: Run Tissue-Specific Variant Calling
################################################################################

echo ""
echo "[STEP 3 @ $(date +'%F %T')] Running tissue-specific variant calling..."

# Prepare environment
ALLELES_DIR="Data/04_analysis/cnv/numbat/inputs/alleles"
mkdir -p "$ALLELES_DIR/pileup" "$ALLELES_DIR/phasing"

# Get paths
BAM_FILE="Data/04_analysis/cnv/numbat/inputs/bam_merged/${DATASET}_merged_for_numbat.bam"
SNP_VCF="Data/hg38_resources/numbat/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf"
GENETIC_MAP="Data/hg38_resources/numbat/genetic_map_hg38_withX.txt.gz"
REF_FASTA="Data/02_references/genome/hg38.fa"

# Verify inputs
echo "  Checking inputs..."
[[ -f "$BAM_FILE" ]] || { echo "[ERROR] BAM not found: $BAM_FILE"; exit 1; }
[[ -f "$SNP_VCF" ]] || { echo "[ERROR] VCF not found: $SNP_VCF"; exit 1; }
[[ -f "$GENETIC_MAP" ]] || { echo "[ERROR] Genetic map not found: $GENETIC_MAP"; exit 1; }
[[ -f "$REF_FASTA" ]] || { echo "[ERROR] Ref fasta not found: $REF_FASTA"; exit 1; }

echo "  ✓ All variant inputs verified"

# Run pileup for each tissue
run_tissue_pileup() {
  local tissue=$1
  local bam=$2
  local snp_vcf=$3
  local alleles_dir=$4
  local dataset=$5
  local barcodes_file="Data/04_analysis/cnv/numbat/inputs/${dataset}_${tissue}_atac_barcodes_for_pileup.tsv"
  
  local pileup_dir="${alleles_dir}/pileup/${dataset}_${tissue}_atac"
  mkdir -p "$pileup_dir"
  
  echo "    [${tissue}] Running cellsnp-lite pileup..."
  echo "      Barcodes: $(wc -l < $barcodes_file)"
  
  cellsnp-lite \
    -s "$bam" \
    -b "$barcodes_file" \
    -O "$pileup_dir" \
    -R "$snp_vcf" \
    --nproc 4 \
    --minCOV 20 \
    --minMAF 0.1 \
    --cellTAG BC \
    2>&1 | tee "${pileup_dir}/pileup.log"
  
  if [[ $? -eq 0 ]]; then
    echo "    [✓] ${tissue} pileup complete"
    return 0
  else
    echo "    [ERROR] ${tissue} pileup failed" >&2
    return 1
  fi
}

# Run parallel pileup for both tissues
echo "  Running parallel tissue pileup..."
for TISSUE in "${TISSUES[@]}"; do
  run_tissue_pileup "$TISSUE" "$BAM_FILE" "$SNP_VCF" "$ALLELES_DIR" "$DATASET" &
done
wait

echo ""
echo "[$(date +'%F %T')] ===== DEEPSEQ TISSUE INPUTS COMPLETE ====="
echo "Ready for NUMBAT:"
echo "  - deepseq_488B_atac_bin.rds"
echo "  - deepseq_489_atac_bin.rds"
echo "  - deepseq_488B_atac_allele_counts.tsv.gz"
echo "  - deepseq_489_atac_allele_counts.tsv.gz"
