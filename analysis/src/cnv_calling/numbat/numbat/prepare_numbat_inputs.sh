#!/bin/bash -l
#$ -P paxlab
#$ -N numbat_low_489
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=8G
#$ -l h_rt=48:00:00
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_low_489.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_low_489.$JOB_ID.err
#$ -j n

################################################################################
# prepare_numbat_inputs.sh
#
# PURPOSE: Prepare NUMBAT ATAC-bin mode inputs with proper variant calling
#
# INPUTS:
#   - BAM file (merged tissue-specific BAM with cell barcodes)
#   - Fragment file (ATAC fragments TSV/BED format)
#   - Barcode file (tissue-specific barcodes)
#   - SNP VCF (genome1K variants)
#   - Genetic map (for phasing)
#
# OUTPUTS:
#   - Allele count matrix (TSV.GZ)
#   - Binned ATAC matrix (RDS)
#   - Reference matrix (RDS)
#
# USAGE:
#   bash prepare_numbat_inputs.sh <dataset> <tissue>
#   Example: bash prepare_numbat_inputs.sh lowseq 488B
#
################################################################################

set -eo pipefail

# CRITICAL: Initialize module system before using 'module' command
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

module load R
module load samtools

# cellsnp-lite lives in the calicost conda env; eagle binary is in external/
export PATH="/projectnb/paxlab/presh/env/calicost_env/bin:/projectnb/paxlab/presh/software/external/Eagle_v2.4.1:${PATH}"

# Configuration
PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
DATASET="${1:-lowseq}"  # lowseq or deepsep
TISSUE="${2:-488B}"     # 488B or 489
TISSUE="489" # HARDCODED FOR TESTING - CHANGE BACK TO ARGUMENT IN PRODUCTION
NCORES=8

# Paths
NUMBAT_BIN="/projectnb/paxlab/presh/Rlibs/4.5/numbat/bin"
NUMBAT_EXTDATA="/projectnb/paxlab/presh/Rlibs/4.5/numbat/extdata"

# Input paths (using absolute paths for better compatibility)
# BAM files may be in different locations depending on dataset
if [[ "$DATASET" == "deepseq" ]]; then
  # Deepseq BAMs are in bam_merged/ subdirectory
  BAM_FILE="$PROJECT_ROOT/Data/04_analysis/cnv/numbat/inputs/bam_merged/${DATASET}_${TISSUE}/bam/${DATASET}_${TISSUE}_merged_for_numbat.bam"
else
  # Lowseq BAMs are directly in inputs/ folder
  BAM_FILE="$PROJECT_ROOT/Data/01_inputs/bam/${DATASET}_${TISSUE}.bam"
fi

BARCODE_FILE="$PROJECT_ROOT/Data/01_inputs/barcodes/tissue_barcodes/${DATASET}_${TISSUE}/${DATASET}_${TISSUE}.barcodes.tsv"
FRAGMENT_FILE="$PROJECT_ROOT/Data/01_inputs/fragments/${DATASET}_${TISSUE}/${DATASET}_${TISSUE}.fragments.sort.filtered_std.bed.gz"
SNP_VCF="$PROJECT_ROOT/Data/02_references/genome/hg38_resources/numbat/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
GENETIC_MAP="$PROJECT_ROOT/Data/02_references/genome/hg38_resources/numbat/genetic_map_hg38_withX.txt.gz"
# BINGR="$NUMBAT_EXTDATA/var220kb.rds"
BINGR="/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/cnv_calling/numbat/numbat/var1Mb.rds"

# Output paths
OUTPUT_DIR="$PROJECT_ROOT/Data/04_analysis/cnv/numbat/inputs/${DATASET}_${TISSUE}"
ALLELE_OUTPUT="$OUTPUT_DIR/${DATASET}_${TISSUE}_atac_allele_counts.tsv.gz"
ATAC_BIN_OUTPUT="$PROJECT_ROOT/Data/04_analysis/cnv/numbat/inputs/${DATASET}_${TISSUE}/${DATASET}_${TISSUE}_atac_bin.rds"
RESULT_DIR="$PROJECT_ROOT/Data/04_analysis/cnv/numbat/results/${DATASET}_${TISSUE}/atac_bin_run3"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$RESULT_DIR"

cd "$PROJECT_ROOT"

echo "[STAGE 1] Verifying inputs for $DATASET $TISSUE..."
echo ""

# CRITICAL VALIDATION: Ensure correct tissue-specific inputs
echo "  Tissue Assignment:"
echo "    DATASET: $DATASET"
echo "    TISSUE: $TISSUE"
echo ""

# Check all required inputs exist
INPUTS_VALID=true
for input in "$BAM_FILE" "$BARCODE_FILE" "$SNP_VCF" "$GENETIC_MAP" "$BINGR"; do
  if [[ ! -f "$input" ]]; then
    echo "[ERROR] Missing: $input"
    INPUTS_VALID=false
  else
    SIZE=$(du -sh "$input" | awk '{print $1}')
    echo "  ✓ $(basename $input): $SIZE"
  fi
done

if [[ "$INPUTS_VALID" != true ]]; then
  echo "[ERROR] Some inputs missing. Aborting."
  exit 1
fi

# TISSUE-SPECIFIC VALIDATION
echo ""
echo "  Barcode Count Validation:"
BARCODE_COUNT=$(wc -l < "$BARCODE_FILE")
case "$TISSUE" in
  488B)
    EXPECTED_BARCODES=11645
    ;;
  489)
    EXPECTED_BARCODES=4671
    ;;
  *)
    echo "[ERROR] Unknown tissue: $TISSUE"
    exit 1
    ;;
esac

echo "    Expected barcodes for $TISSUE: $EXPECTED_BARCODES"
echo "    Actual barcodes in file: $BARCODE_COUNT"

if [[ $BARCODE_COUNT -ne $EXPECTED_BARCODES ]]; then
  echo "[ERROR] TISSUE MISMATCH! Got $BARCODE_COUNT but expected $EXPECTED_BARCODES for $TISSUE"
  echo "[ERROR] This suggests wrong barcode file was provided"
  exit 1
fi
echo "    ✓ Barcode count CORRECT for $TISSUE"
echo ""

echo ""
echo "[STAGE 2] Running pileup and phase (variant calling)..."


# TISSUE ISOLATION: Remove any old pileup data for this tissue to prevent contamination
PILEUP_DIR="$OUTPUT_DIR/pileup_phase"
# if [[ -d "$PILEUP_DIR" ]]; then
#     log_warn "Removing old pileup directory: $PILEUP_DIR"
#     rm -rf "$PILEUP_DIR"
#     log_info "✓ Removed (ensures fresh start)"
# fi

# mkdir -p "$PILEUP_DIR"
# log_info "✓ Created pileup directory: $PILEUP_DIR"


# Use pre-built reference panel for Eagle phasing
PANEL_DIR="$PROJECT_ROOT/Data/02_references/genome/hg38_resources/numbat/1000G_hg38"
if [[ ! -d "$PANEL_DIR" ]]; then
  echo "[ERROR] Panel BCF directory not found: $PANEL_DIR"
  # exit 1
fi

# Load R module
module load R 2>/dev/null || true

# Step 1: Pileup and phase to get allele counts
echo "  Running: pileup_and_phase.R"
Rscript "$NUMBAT_BIN/pileup_and_phase.R" \
  --label "${DATASET}_${TISSUE}" \
  --samples "${DATASET}_${TISSUE}" \
  --bams "$BAM_FILE" \
  --barcodes "$BARCODE_FILE" \
  --snpvcf "$SNP_VCF" \
  --gmap "$GENETIC_MAP" \
  --eagle "/projectnb/paxlab/presh/software/external/Eagle_v2.4.1/eagle" \
  --paneldir "$PANEL_DIR" \
  --ncores $NCORES \
  --cellTAG CB \
  --UMItag None \
  --outdir "$PILEUP_DIR" 


echo ""
echo "[STAGE 3] Generating binned ATAC matrix..."
echo "  Running: get_binned_atac.R"

cd "$OUTPUT_DIR"
Rscript "$NUMBAT_BIN/get_binned_atac.R" \
  --CB "$BARCODE_FILE" \
  --frag "$FRAGMENT_FILE" \
  --binGR "$BINGR" \
  --outFile "$ATAC_BIN_OUTPUT"  


echo "  ✓ Binned ATAC matrix: $ATAC_BIN_OUTPUT"

# echo ""
# echo "[STAGE 4] Generate reference matrix..."

REF_DIR="$PROJECT_ROOT/Data/04_analysis/cnv/numbat/inputs/${DATASET}_${TISSUE}/reference"


# R_LIBS_USER="/projectnb/paxlab/presh/Rlibs/4.5":$R_LIBS_USER


echo ""
echo "[STAGE 5] Running numbat multiome mode for CNV calling..."
parL="/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/cnv/numbat/reference/par_numbatm_run2.rds"

Rscript "$NUMBAT_BIN/run_numbat_multiome.R"  \
            --countmat "$OUTPUT_DIR/${DATASET}_${TISSUE}_atac_bin.rds" \
            --alleledf "$PILEUP_DIR/${DATASET}_${TISSUE}_allele_counts.tsv.gz" \
            --out_dir "$RESULT_DIR" \
            --ref "$REF_DIR/lambdas_ATAC_bincnt.rds" \
            --gtf "$BINGR" --parL "$parL" > "$RESULT_DIR/numbat_run.log" 2>&1

