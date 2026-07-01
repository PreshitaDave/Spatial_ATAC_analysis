#!/bin/bash
set -eo pipefail

################################################################################
# merge_bam_for_numbat.qsub.sh
# 
# Purpose: Merge chromosome-split BAM files from monopogen variant calling
#          into a single merged BAM for NUMBAT analysis with correct BC16 format
#
# Input:   Chromosome-split BAMs (deepseq_chr1.filter.bam, etc.) from
#          Data/04_analysis/03_intermediate/variant_calling/monopogen/variant_calling/deepseq/Bam/
#
# Output:  Merged BAM: Data/04_analysis/cnv/numbat/inputs/deepseq_merged_for_numbat_bc16.bam
#          Index:      Data/04_analysis/cnv/numbat/inputs/deepseq_merged_for_numbat_bc16.bam.bai
#
# Author:  Spatial ATAC Pipeline
# Date:    May 16, 2026
#
# SGE CONFIGURATION:
#$ -N merge_bam_numbat
#$ -pe omp 8
#$ -l h_rt=48:00:00
#$ -l mem_per_core=8G
#$ -P paxlab
#$ -e analysis/qsub_logs/merge_bam_numbat.err
#$ -o analysis/qsub_logs/merge_bam_numbat.out
#
################################################################################

# Initialize module system
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

# Load required modules
module load samtools/1.14 || module load samtools

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
NCORES=${NSLOTS:-8}

# Logging function
log() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" >&2
}

run_cmd() {
  log "EXEC: $*"
  "$@" || { log "ERROR: Command failed: $*"; exit 1; }
}

# Verify we're on a compute node (not login node)
HOSTNAME=$(hostname)
if [[ "$HOSTNAME" == scc1* ]]; then
  log "ERROR: Running on login node ($HOSTNAME). Use qrsh to get compute node."
  exit 1
fi
log "Running on compute node: $HOSTNAME"

################################################################################
# STEP 1: Verify Input Files
################################################################################

log "========== STEP 1: VERIFY INPUT FILES =========="

BAM_SOURCE_DIR="${PROJECT_ROOT}/Data/04_analysis/03_intermediate/variant_calling/monopogen/variant_calling/deepseq/Bam"

if [[ ! -d "$BAM_SOURCE_DIR" ]]; then
  log "ERROR: BAM source directory not found: $BAM_SOURCE_DIR"
  exit 1
fi

# Count chromosome BAMs
CHR_BAM_COUNT=$(find "$BAM_SOURCE_DIR" -name "deepseq_chr*.filter.bam" -type f | wc -l)
log "Found $CHR_BAM_COUNT chromosome BAMs in $BAM_SOURCE_DIR"

if [[ $CHR_BAM_COUNT -ne 22 ]]; then
  log "WARNING: Expected 22 chromosome BAMs, found $CHR_BAM_COUNT"
fi

# List BAMs that will be merged (sorted by chromosome number)
log "BAMs to merge:"
for chr in {1..22}; do
  bam_file="${BAM_SOURCE_DIR}/deepseq_chr${chr}.filter.bam"
  if [[ -f "$bam_file" ]]; then
    size_gb=$(du -h "$bam_file" | cut -f1)
    log "  chr${chr}: $size_gb"
  else
    log "  ERROR: Missing deepseq_chr${chr}.filter.bam"
    exit 1
  fi
done

################################################################################
# STEP 2: Prepare Output Directory
################################################################################

log "========== STEP 2: PREPARE OUTPUT DIRECTORY =========="

OUTPUT_DIR="${PROJECT_ROOT}/Data/04_analysis/cnv/numbat/inputs"
OUTPUT_BAM="${OUTPUT_DIR}/deepseq_merged_for_numbat_bc16.bam"
OUTPUT_BAI="${OUTPUT_DIR}/deepseq_merged_for_numbat_bc16.bam.bai"

# Create output directory if needed
if [[ ! -d "$OUTPUT_DIR" ]]; then
  log "Creating output directory: $OUTPUT_DIR"
  mkdir -p "$OUTPUT_DIR" || { log "ERROR: Failed to create $OUTPUT_DIR"; exit 1; }
fi

log "Output BAM will be written to: $OUTPUT_BAM"

# Check if output already exists
if [[ -f "$OUTPUT_BAM" ]]; then
  log "WARNING: Output BAM already exists ($(du -h "$OUTPUT_BAM" | cut -f1))"
  log "Backup existing file before merge"
  BACKUP_BAM="${OUTPUT_BAM}.backup.$(date +%s)"
  run_cmd mv "$OUTPUT_BAM" "$BACKUP_BAM"
  log "Backed up to: $BACKUP_BAM"
fi

################################################################################
# STEP 3: Merge Chromosome BAMs
################################################################################

log "========== STEP 3: MERGE CHROMOSOME BAMS =========="

# Build list of BAMs in order (chr1, chr2, ..., chr22)
BAM_LIST=()
for chr in {1..22}; do
  BAM_LIST+=("${BAM_SOURCE_DIR}/deepseq_chr${chr}.filter.bam")
done

log "Merging ${#BAM_LIST[@]} chromosome BAMs..."
log "Output BAM: $OUTPUT_BAM"
log "Using $NCORES cores"

# Run samtools merge with threading
START_TIME=$(date +%s)
run_cmd samtools merge \
  -f \
  -h "${BAM_LIST[0]}" \
  -@ $((NCORES - 1)) \
  "$OUTPUT_BAM" \
  "${BAM_LIST[@]}"

END_TIME=$(date +%s)
MERGE_TIME=$((END_TIME - START_TIME))
log "BAM merge completed in ${MERGE_TIME}s ($(( MERGE_TIME / 60 ))m)"

################################################################################
# STEP 4: Verify Merged BAM
################################################################################

log "========== STEP 4: VERIFY MERGED BAM =========="

if [[ ! -f "$OUTPUT_BAM" ]]; then
  log "ERROR: Merged BAM file not found: $OUTPUT_BAM"
  exit 1
fi

OUTPUT_SIZE_GB=$(du -h "$OUTPUT_BAM" | cut -f1)
log "Merged BAM size: $OUTPUT_SIZE_GB"

# Index the merged BAM
log "Indexing merged BAM..."
run_cmd samtools index -@ $((NCORES - 1)) "$OUTPUT_BAM" "$OUTPUT_BAI"

if [[ ! -f "$OUTPUT_BAI" ]]; then
  log "ERROR: BAM index not created: $OUTPUT_BAI"
  exit 1
fi

INDEX_SIZE=$(du -h "$OUTPUT_BAI" | cut -f1)
log "BAM index size: $INDEX_SIZE"

################################################################################
# STEP 5: Sample and Report CB Tags
################################################################################

log "========== STEP 5: SAMPLE AND REPORT CB TAGS =========="

# Sample first few reads to check CB tag format
log "Sampling first 5 reads with CB tags:"
SAMPLE_COUNT=0
samtools view "$OUTPUT_BAM" 2>/dev/null | while read -r line && [[ $SAMPLE_COUNT -lt 5 ]]; do
  CB_TAG=$(echo "$line" | tr '\t' '\n' | grep "^CB:")
  if [[ -n "$CB_TAG" ]]; then
    log "  CB tag: $CB_TAG"
    ((SAMPLE_COUNT++))
  fi
done

# Get unique CB tags count
log "Counting unique CB tags..."
UNIQUE_CB_COUNT=$(samtools view "$OUTPUT_BAM" 2>/dev/null | awk '{for (i=1; i<=NF; i++) if ($i ~ /^CB:/) print $i}' | sort -u | wc -l)
log "Unique CB tags in merged BAM: $UNIQUE_CB_COUNT"

################################################################################
# STEP 6: Summary
################################################################################

log "========== STEP 6: SUMMARY =========="
log "✓ BAM merge completed successfully"
log "  Input:   22 chromosome BAMs from ${BAM_SOURCE_DIR}/"
log "  Output:  ${OUTPUT_BAM}"
log "  Size:    $OUTPUT_SIZE_GB"
log "  Index:   $OUTPUT_BAI"
log "  CB tags: $UNIQUE_CB_COUNT unique barcodes"
log ""
log "Next step: Extract BC16 barcodes and tissue filters"
log "  Script:  analysis/src/cnv_calling/numbat/extract_barcodes_from_bam.sh"
log ""

log "MERGE COMPLETE"
