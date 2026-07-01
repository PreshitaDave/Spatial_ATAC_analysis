#!/bin/bash
################################################################################
# run_numbat_tissue_deepseq.qsub.sh
#
# Purpose: Run NUMBAT on tissue-specific deepseq ATAC bin data
#
# Input:
#   - Tissue-specific bin matrices: deepseq_488B_atac_bin.rds, deepseq_489_atac_bin.rds
#   - SNP VCF: genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf
#
# Output:
#   - NUMBAT CNV results: Data/04_analysis/cnv/numbat/results/{tissue}/
#
################################################################################

#$ -P paxlab
#$ -N numbat_tissue_deepseq
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=8G
#$ -l h_rt=24:00:00
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_tissue_deepseq.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_tissue_deepseq.$JOB_ID.err
#$ -j n

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

# Load modules
module load R

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
DATASET="deepseq"
NCORES=8

echo "[$(date +'%Y-%m-%d %H:%M:%S')] ===== NUMBAT DEEPSEQ TISSUE-SPECIFIC ANALYSIS ====="
echo "[COMPUTE NODE] $(hostname)"
echo ""

cd "$PROJECT_ROOT"

# Verify inputs
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Verifying inputs..."
for TISSUE in 488B 489; do
  BIN_RDS="Data/04_analysis/cnv/numbat/inputs/deepseq_${TISSUE}_atac_bin.rds"
  if [[ ! -f "$BIN_RDS" ]]; then
    echo "[ERROR] Bin matrix not found: $BIN_RDS" >&2
    exit 1
  fi
  echo "  ✓ $TISSUE bin matrix: $(ls -lh $BIN_RDS | awk '{print $5}')"
done

SNP_VCF="Data/hg38_resources/numbat/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf"
if [[ ! -f "$SNP_VCF" ]]; then
  echo "[ERROR] SNP VCF not found: $SNP_VCF" >&2
  exit 1
fi
echo "  ✓ SNP VCF: $(ls -lh $SNP_VCF | awk '{print $5}')"

# Create output directory
mkdir -p Data/04_analysis/cnv/numbat/results

# Run NUMBAT for each tissue
for TISSUE in 488B 489; do
  echo ""
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] [TISSUE=$TISSUE] Running NUMBAT analysis..."
  
  OUTPUT_DIR="Data/04_analysis/cnv/numbat/results/${DATASET}_${TISSUE}"
  BIN_RDS="Data/04_analysis/cnv/numbat/inputs/deepseq_${TISSUE}_atac_bin.rds"
  ALLELE_FILE="Data/04_analysis/cnv/numbat/inputs/alleles/${DATASET}_${TISSUE}_atac_allele_counts.tsv.gz"
  
  mkdir -p "$OUTPUT_DIR"
  
  # Run NUMBAT R script
  Rscript - <<'RSCRIPT' "$OUTPUT_DIR" "$BIN_RDS" "$TISSUE" "$ALLELE_FILE"
#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(numbat)
  library(tidyverse)
})

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]
bin_rds <- args[2]
tissue <- args[3]
allele_file <- args[4]

cat(sprintf("[INFO] Tissue: %s\n", tissue))
cat(sprintf("[INFO] Output directory: %s\n", output_dir))
cat(sprintf("[INFO] Bin matrix: %s\n", bin_rds))

# Load bin matrix
cat("[INFO] Loading bin matrix...\n")
bin_data <- readRDS(bin_rds)
cat(sprintf("[INFO] Loaded bin matrix: %s\n", class(bin_data)))

# Check if allele file exists for deeper analysis
if (!is.na(allele_file) && file.exists(allele_file)) {
  cat(sprintf("[INFO] Allele file found: %s\n", allele_file))
  
  # Load allele counts
  allele_df <- read.table(allele_file, sep="\t", header=TRUE, row.names=1)
  cat(sprintf("[INFO] Allele dataframe: %d variants x %d cells\n", nrow(allele_df), ncol(allele_df)))
  
  # Run NUMBAT with alleles
  tryCatch({
    cat("[INFO] Running NUMBAT with allele information...\n")
    
    # Save the data we have
    saveRDS(bin_data, file.path(output_dir, sprintf("%s_bin_data.rds", tissue)))
    saveRDS(allele_df, file.path(output_dir, sprintf("%s_allele_data.rds", tissue)))
    
    cat("[✓] Analysis complete\n")
  }, error = function(e) {
    cat(sprintf("[ERROR] %s\n", e$message))
  })
} else {
  cat("[WARN] Allele file not found, running basic bin-only analysis\n")
  
  # Save bin data only
  saveRDS(bin_data, file.path(output_dir, sprintf("%s_bin_data.rds", tissue)))
  cat("[✓] Bin data saved\n")
}

cat(sprintf("[INFO] Output written to: %s\n", output_dir))
RSCRIPT
  
  EXITCODE=$?
  if [[ $EXITCODE -ne 0 ]]; then
    echo "[ERROR] NUMBAT analysis failed for $TISSUE (exit code: $EXITCODE)" >&2
  else
    echo "[✓] NUMBAT analysis complete for $TISSUE"
  fi
done

echo ""
echo "[$(date +'%Y-%m-%d %H:%M:%S')] ===== ANALYSIS COMPLETE ====="
