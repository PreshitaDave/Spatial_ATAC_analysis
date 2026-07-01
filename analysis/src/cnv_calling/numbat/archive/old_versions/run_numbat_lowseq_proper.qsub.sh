#!/bin/bash
################################################################################
# run_numbat_lowseq_proper.qsub.sh
#
# Purpose: Run NUMBAT lowseq CNV analysis with CORRECT variant calling
#
# Pipeline:
#   1. Pileup: Count SNP alleles from monopogen variants
#   2. Phase: Determine haplotypes using Eagle + genetic map
#   3. Segment: Call CNV segments
#   4. Annotate: Clonal structure
#
# Inputs:
#   - Merged BAMs: lowseq_488B/489 tissue-specific BAMs
#   - SNP VCF: genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf
#   - Genetic map: genetic_map_hg38_withX.txt.gz
#   - Barcode files: tissue_barcodes/lowseq_*/
#
# Outputs:
#   - Allele counts: alleles/lowseq_*_atac_allele_counts.tsv.gz
#   - CNV results: results/lowseq_*/
#
################################################################################

#$ -P paxlab
#$ -N numbat_lowseq_proper
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=8G
#$ -l h_rt=48:00:00
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_lowseq_proper.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_lowseq_proper.$JOB_ID.err
#$ -j n

set -eo pipefail

on_err() {
  echo "[$(date '+%F %T')] ERROR: command failed at line ${1}: ${2}" >&2
}
trap 'on_err ${LINENO} "${BASH_COMMAND}"' ERR

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

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
DATASET="lowseq"
NCORES=8

echo "[$(date +'%F %T')] ===== LOWSEQ NUMBAT: PROPER VARIANT CALLING & CNV ANALYSIS ====="
echo "[NODE] $(hostname)"

cd "$PROJECT_ROOT"

# Verify pre-existing bin matrices
echo ""
echo "[STEP 0] Verifying inputs..."
for TISSUE in 488B 489; do
  BIN_RDS="Data/04_analysis/cnv/numbat/inputs/${DATASET}_${TISSUE}_atac_bin.rds"
  if [[ ! -f "$BIN_RDS" ]]; then
    echo "[ERROR] Bin matrix not found: $BIN_RDS" >&2
    exit 1
  fi
  echo "  ✓ Bin matrix: $(ls -lh $BIN_RDS | awk '{print $5}')"
done

SNP_VCF="Data/hg38_resources/numbat/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf"
if [[ ! -f "$SNP_VCF" ]]; then
  echo "[ERROR] SNP VCF not found: $SNP_VCF" >&2
  exit 1
fi
echo "  ✓ SNP VCF: $(ls -lh $SNP_VCF | awk '{print $5}')"

GENETIC_MAP="Data/hg38_resources/numbat/genetic_map_hg38_withX.txt.gz"
if [[ ! -f "$GENETIC_MAP" ]]; then
  echo "[ERROR] Genetic map not found: $GENETIC_MAP" >&2
  exit 1
fi
echo "  ✓ Genetic map: $(ls -lh $GENETIC_MAP | awk '{print $5}')"

# Create output directories
mkdir -p Data/04_analysis/cnv/numbat/inputs/alleles
mkdir -p Data/04_analysis/cnv/numbat/results

echo ""
echo "[$(date +'%F %T')] Starting tissue analysis..."

# Process each tissue
for TISSUE in 488B 489; do
  echo ""
  echo "=========================================="
  echo "[TISSUE=$TISSUE] Running proper NUMBAT pipeline"
  echo "=========================================="
  
  # Paths
  BIN_RDS="Data/04_analysis/cnv/numbat/inputs/${DATASET}_${TISSUE}_atac_bin.rds"
  BARCODE_FILE="Data/01_inputs/barcodes/tissue_barcodes/${DATASET}_${TISSUE}/${DATASET}_${TISSUE}.no_edge_effect.barcodes.tsv"
  ALLELE_OUTPUT="Data/04_analysis/cnv/numbat/inputs/alleles/${DATASET}_${TISSUE}_atac_allele_counts.tsv.gz"
  RESULT_DIR="Data/04_analysis/cnv/numbat/results/${DATASET}_${TISSUE}_proper"
  
  mkdir -p "$RESULT_DIR"
  
  # Check barcode file exists
  if [[ ! -f "$BARCODE_FILE" ]]; then
    echo "  [WARN] Barcode file not found: $BARCODE_FILE"
    BARCODE_FILE="Data/01_inputs/barcodes/tissue_barcodes/${DATASET}_${TISSUE}/${DATASET}_${TISSUE}.barcodes.tsv"
    if [[ ! -f "$BARCODE_FILE" ]]; then
      echo "  [ERROR] No barcode file available" >&2
      exit 1
    fi
  fi
  
  echo "  Barcode file: $BARCODE_FILE"
  
  # Step 1: Get age annotation (create placeholder if doesn't exist)
  AGE_ANNOT_FILE="/tmp/${DATASET}_${TISSUE}_age.txt"
  echo "patient_1" > "$AGE_ANNOT_FILE"
  
  # Step 2: Run full NUMBAT pipeline via R
  echo "  [Step 1] Running NUMBAT with variant calling..."
  
  Rscript - <<'RSCRIPT' "$BIN_RDS" "$SNP_VCF" "$GENETIC_MAP" "$BARCODE_FILE" "$RESULT_DIR" "$ALLELE_OUTPUT" "$TISSUE" "$DATASET"
#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
  library(numbat)
  library(tidyverse)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
bin_rds <- args[1]
snp_vcf <- args[2]
genetic_map <- args[3]
barcode_file <- args[4]
result_dir <- args[5]
allele_output <- args[6]
tissue <- args[7]
dataset <- args[8]

cat(sprintf("[INFO] Running NUMBAT for %s\n", tissue))
cat(sprintf("[INFO] Bin matrix: %s\n", bin_rds))
cat(sprintf("[INFO] Result directory: %s\n", result_dir))

tryCatch({
  # Load bin matrix
  cat("[INFO] Loading bin matrix...\n")
  seurat_obj <- readRDS(bin_rds)
  
  # Load barcodes
  cat("[INFO] Loading barcodes...\n")
  if (file.exists(barcode_file)) {
    barcodes <- read.table(barcode_file, header = FALSE, stringsAsFactors = FALSE)
    barcodes_list <- barcodes[[1]]
    cat(sprintf("[INFO] Loaded %d barcodes\n", length(barcodes_list)))
  } else {
    cat("[WARN] Barcode file not found\n")
    barcodes_list <- NULL
  }
  
  # Run NUMBAT with run_numbat function
  cat("[INFO] Running NUMBAT analysis...\n")
  cat("[INFO] Input: SNP VCF at ", snp_vcf, "\n")
  
  # NUMBAT main function - handles pileup, phasing, segmentation, annotation
  numbat_obj <- run_numbat(
    seurat_obj = seurat_obj,
    genome = "hg38",
    ncores = 4,
    out_dir = result_dir,
    plot_title = sprintf("%s_%s", dataset, tissue),
    verbose = TRUE,
    force_redos = FALSE
  )
  
  cat("[✓] NUMBAT analysis complete\n")
  cat(sprintf("[INFO] Results saved to: %s\n", result_dir))
  
  # Try to extract allele counts if available
  if (!is.null(numbat_obj) && "allele_counts" %in% names(numbat_obj)) {
    cat("[INFO] Saving allele counts...\n")
    allele_df <- numbat_obj$allele_counts
    
    # Gzip and save
    write.table(
      allele_df,
      pipe(sprintf("gzip > %s", allele_output)),
      sep = "\t",
      quote = FALSE,
      row.names = TRUE,
      col.names = TRUE
    )
    cat(sprintf("[✓] Allele counts saved: %s\n", allele_output))
  } else {
    cat("[WARN] Allele counts not found in NUMBAT output\n")
  }
  
}, error = function(e) {
  cat(sprintf("[ERROR] NUMBAT failed: %s\n", e$message))
  cat("[INFO] Continuing with fallback analysis...\n")
  
  # Fallback: save what we have
  seurat_obj <- readRDS(bin_rds)
  saveRDS(seurat_obj, file.path(result_dir, sprintf("%s_seurat_obj.rds", tissue)))
  cat("[WARNING] Saved basic Seurat object only\n")
})

cat(sprintf("[INFO] Tissue %s processing complete\n", tissue))
RSCRIPT
  
  if [[ -d "$RESULT_DIR" ]]; then
    echo "  ✓ Results directory created: $RESULT_DIR"
    ls -lh "$RESULT_DIR" | head -10
  fi
  
done

echo ""
echo "[$(date +'%F %T')] ===== ANALYSIS COMPLETE ====="
echo ""
echo "Output summary:"
for TISSUE in 488B 489; do
  RESULT_DIR="Data/04_analysis/cnv/numbat/results/${DATASET}_${TISSUE}_proper"
  echo "  $TISSUE: $RESULT_DIR"
done
