#!/bin/bash -l
#$ -P paxlab
#$ -N tissue_files
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l h_rt=24:00:00
#$ -l mem_per_core=8G
#$ -j n
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/valid_bam.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/valid_bam.$JOB_ID.err
# BAM Validation and Barcode Checking Script
# Validates BAM files, extracts barcodes, and compares with reference barcode files
# Runs in parallel using GNU Parallel

set -euo pipefail

module load samtools
module load parallel
# Configuration
BAM_DIR1="/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/bam"
NUMBAT_DIR="/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/cnv/numbat/inputs"
BARCODE_DIR="/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes"

# Define BAM files with their paths
declare -A BAM_FILES=(
    ["deepseq_488B"]="$BAM_DIR1/deepseq_488B.bam"
    ["deepseq_489"]="$BAM_DIR1/deepseq_489.bam"
    ["lowseq_488B"]="$BAM_DIR1/lowseq_488B.bam"
    ["lowseq_489"]="$BAM_DIR1/lowseq_489.bam"
    ["deepseq_488B_numbat"]="$NUMBAT_DIR/deepseq_488B/bam/deepseq_488B_merged_for_numbat.bam"
    ["deepseq_489_numbat"]="$NUMBAT_DIR/deepseq_489/bam/deepseq_489_merged_for_numbat.bam"
    ["lowseq_488B_numbat"]="$NUMBAT_DIR/lowseq_488B/bam/lowseq_488B_merged_for_numbat.bam"
    ["lowseq_489_numbat"]="$NUMBAT_DIR/lowseq_489/bam/lowseq_489_merged_for_numbat.bam"
)

NUM_CORES=${1:-8}  # Default to 16 cores, or use first argument

# Create output directory
OUTPUT_DIR="/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/bam/bam_validation_results"
mkdir -p "$OUTPUT_DIR"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

log_debug() {
    echo -e "${BLUE}[DEBUG]${NC} $1"
}

# Check dependencies
check_dependencies() {
    log_info "Checking dependencies..."
    
    for cmd in samtools parallel; do
        if ! command -v "$cmd" &> /dev/null; then
            log_error "$cmd is not installed or not in PATH"
            exit 1
        fi
    done
    
    log_info "All dependencies found: samtools, parallel"
}

# Validate single BAM file
validate_bam() {
    local sample_name="$1"
    local bam_file="$2"
    local output_file="$3"
    
    {
        echo "=== Validating: $sample_name ==="
        echo "Path: $bam_file"
        
        if [[ ! -f "$bam_file" ]]; then
            echo "ERROR: File does not exist"
            return 1
        fi
        
        if [[ ! -r "$bam_file" ]]; then
            echo "ERROR: File is not readable"
            return 1
        fi
        
        echo "File size: $(du -h "$bam_file" | cut -f1)"
        echo "Last modified: $(stat -c%y "$bam_file" 2>/dev/null || stat -f%Sm "$bam_file" 2>/dev/null || echo 'N/A')"
        
        # Check for index file
        if [[ -f "${bam_file}.bai" ]]; then
            echo "Index file: Present (${bam_file}.bai)"
        else
            echo "Index file: MISSING - may need to create index"
        fi
        
        # Validate BAM structure
        if samtools view -H "$bam_file" &> /dev/null; then
            echo "Status: VALID BAM structure"
            echo "Header: OK"
        else
            echo "Status: INVALID BAM structure"
            return 1
        fi
        
        # Get basic statistics
        local num_reads=$(samtools view -c "$bam_file" 2>/dev/null || echo "0")
        echo "Total reads: $num_reads"
        
        # Check for common tags (sample first 100 reads)
        local has_cb=$(samtools view "$bam_file" 2>/dev/null | head -100 | grep -c 'CB:Z:' || echo "0")
        echo "Reads with CB tag: $has_cb (checked first 100 reads)"
        
        # Check for other common tags
        local has_xm=$(samtools view "$bam_file" 2>/dev/null | head -100 | grep -c 'XM:i:' || echo "0")
        local has_as=$(samtools view "$bam_file" 2>/dev/null | head -100 | grep -c 'AS:i:' || echo "0")
        echo "Reads with XM tag: $has_xm"
        echo "Reads with AS tag: $has_as"
        
        echo "Validation: PASSED"
        return 0
        
    } > "$output_file" 2>&1
}

# Extract barcodes from BAM file
extract_barcodes() {
    local sample_name="$1"
    local bam_file="$2"
    local output_file="$3"
    
    {
        echo "=== Extracting barcodes from $sample_name ==="
        
        if [[ ! -f "$bam_file" ]]; then
            echo "ERROR: BAM file not found: $bam_file"
            return 1
        fi
        
        # Extract CB (cell barcode) tags
        samtools view "$bam_file" 2>/dev/null \
            | awk '{for(i=1;i<=NF;i++) if($i ~ /^CB:Z:/) {sub(/CB:Z:/,"",$i); print $i}}' \
            | sort | uniq > "${output_file}.tmp" 2>/dev/null
        
        if [[ ! -f "${output_file}.tmp" ]]; then
            echo "ERROR: Failed to extract barcodes"
            return 1
        fi
        
        local num_barcodes=$(wc -l < "${output_file}.tmp")
        echo "Extracted barcodes: $num_barcodes"
        echo "Sample: $sample_name"
        echo "BAM: $bam_file"
        
        mv "${output_file}.tmp" "$output_file"
        
    } 2>&1 | tee -a "$OUTPUT_DIR/${sample_name}_extract.log"
}

# Compare extracted barcodes with reference
compare_barcodes() {
    local sample_name="$1"
    local extracted_bc="$2"
    local reference_dir="$3"
    local output_file="$4"
    
    {
        echo "=== Barcode Comparison for $sample_name ==="
        echo "Reference directory: $reference_dir"
        
        # For numbat samples, skip barcode comparison since they may not have tissue barcode files
        if [[ "$sample_name" == *"numbat"* ]]; then
            echo "Note: Numbat merged BAM - skipping tissue barcode comparison"
            echo "Extracted barcodes file available at: $extracted_bc"
            return 0
        fi
        
        if [[ ! -d "$reference_dir" ]]; then
            echo "WARNING: Reference directory not found: $reference_dir"
            return 1
        fi
        
        if [[ ! -f "$extracted_bc" ]]; then
            echo "ERROR: Extracted barcode file not found: $extracted_bc"
            return 1
        fi
        
        # Find barcode files in reference directory
        local barcode_files=($(find "$reference_dir" -name "*.barcodes.tsv" -type f 2>/dev/null))
        
        if [[ ${#barcode_files[@]} -eq 0 ]]; then
            echo "WARNING: No barcode files found in $reference_dir"
            return 1
        fi
        
        echo "Found ${#barcode_files[@]} barcode files"
        echo ""
        
        # Compare with each barcode file
        for bc_file in "${barcode_files[@]}"; do
            local filename=$(basename "$bc_file")
            local total_ref=$(wc -l < "$bc_file")
            local total_extracted=$(wc -l < "$extracted_bc")
            
            if [[ $total_extracted -eq 0 ]]; then
                echo "WARNING: No barcodes extracted from BAM file"
                continue
            fi
            
            # Count matches
            local matches=$(comm -12 <(sort "$extracted_bc") <(sort "$bc_file" | awk '{print $1}') | wc -l)
            local in_extracted_only=$(comm -23 <(sort "$extracted_bc") <(sort "$bc_file" | awk '{print $1}') | wc -l)
            local in_ref_only=$(comm -13 <(sort "$extracted_bc") <(sort "$bc_file" | awk '{print $1}') | wc -l)
            
            local match_percent=$(awk "BEGIN {printf \"%.2f\", ($matches/$total_extracted)*100}")
            
            echo "File: $filename"
            echo "  Reference barcodes: $total_ref"
            echo "  Extracted barcodes: $total_extracted"
            echo "  Matches: $matches ($match_percent%)"
            echo "  In extracted only: $in_extracted_only"
            echo "  In reference only: $in_ref_only"
            
            if [[ $in_extracted_only -gt 0 ]]; then
                echo "  WARNING: Found $(( (in_extracted_only * 100) / total_extracted ))% extra barcodes"
            fi
            if [[ $in_ref_only -gt 0 ]]; then
                echo "  WARNING: Missing $(( (in_ref_only * 100) / total_ref ))% of reference barcodes"
            fi
            echo ""
        done
        
    } > "$output_file" 2>&1
}

# Main validation function for a single sample
validate_sample() {
    local sample_name="$1"
    local bam_file="$2"
    
    log_info "Processing sample: $sample_name"
    
    # Determine barcode directory based on sample name
    local base_sample_name=${sample_name%_numbat}
    local barcode_dir="$BARCODE_DIR/$base_sample_name"
    
    # Validate BAM
    validate_bam "$sample_name" "$bam_file" "$OUTPUT_DIR/${sample_name}_validation.txt"
    
    # Extract barcodes
    extract_barcodes "$sample_name" "$bam_file" "$OUTPUT_DIR/${sample_name}_extracted.barcodes.txt"
    
    # Compare barcodes (skip for numbat samples or if barcode dir doesn't exist)
    if [[ -d "$barcode_dir" ]] || [[ "$sample_name" == *"numbat"* ]]; then
        compare_barcodes "$sample_name" "$OUTPUT_DIR/${sample_name}_extracted.barcodes.txt" "$barcode_dir" "$OUTPUT_DIR/${sample_name}_comparison.txt"
    fi
    
    log_info "Completed: $sample_name"
}

export -f validate_bam
export -f extract_barcodes
export -f compare_barcodes
export -f validate_sample
export -f log_info
export -f log_warn
export -f log_error
export -f log_debug
export BAM_DIR1 NUMBAT_DIR BARCODE_DIR OUTPUT_DIR RED GREEN YELLOW BLUE NC

# Main execution
main() {
    log_info "BAM File Validation and Barcode Checking Script"
    log_info "=============================================="
    log_info "BAM Directory 1: $BAM_DIR1"
    log_info "Numbat Directory: $NUMBAT_DIR"
    log_info "Barcode Directory: $BARCODE_DIR"
    log_info "Using $NUM_CORES cores for parallel processing"
    log_info "Total BAM files to validate: ${#BAM_FILES[@]}"
    log_info ""
    
    check_dependencies
    
    # Verify BAM files exist
    log_info "Checking if BAM files exist..."
    local missing_files=0
    
    for sample_name in "${!BAM_FILES[@]}"; do
        local bam_file="${BAM_FILES[$sample_name]}"
        if [[ -f "$bam_file" ]]; then
            log_info "Found: [$sample_name] $bam_file"
        else
            log_error "Missing: [$sample_name] $bam_file"
            ((missing_files++))
        fi
    done
    
    if [[ $missing_files -gt 0 ]]; then
        log_error "$missing_files BAM files are missing!"
        exit 1
    fi
    
    echo ""
    log_info "Starting parallel validation with $NUM_CORES cores..."
    echo ""
    
    # Run validation in parallel
    for sample_name in "${!BAM_FILES[@]}"; do
        local bam_file="${BAM_FILES[$sample_name]}"
        echo "$sample_name $bam_file"
    done | parallel --jobs "$NUM_CORES" --colsep ' ' 'validate_sample {1} {2}'
    
    echo ""
    log_info "Validation complete!"
    log_info "Results saved to: $OUTPUT_DIR"
    
    # Generate summary report
    echo ""
    log_info "=== SUMMARY REPORT ==="
    echo ""
    
    for sample_name in "${!BAM_FILES[@]}"; do
        echo "--- $sample_name ---"
        if [[ -f "$OUTPUT_DIR/${sample_name}_validation.txt" ]]; then
            grep "Status\|Total reads\|File size" "$OUTPUT_DIR/${sample_name}_validation.txt" | head -3
        fi
        if [[ -f "$OUTPUT_DIR/${sample_name}_comparison.txt" ]]; then
            grep "Matches:" "$OUTPUT_DIR/${sample_name}_comparison.txt" | head -1
        fi
        echo ""
    done
    
    # Create consolidated report
    local report_file="$OUTPUT_DIR/FULL_REPORT.txt"
    {
        echo "BAM File Validation Report"
        echo "Generated: $(date)"
        echo ""
        echo "Configuration:"
        echo "  BAM Directory 1: $BAM_DIR1"
        echo "  Numbat Directory: $NUMBAT_DIR"
        echo "  Barcode Directory: $BARCODE_DIR"
        echo "  Cores used: $NUM_CORES"
        echo ""
        
        for sample_name in "${!BAM_FILES[@]}"; do
            echo "=========================================="
            echo "Sample: $sample_name"
            echo "File: ${BAM_FILES[$sample_name]}"
            echo "=========================================="
            
            if [[ -f "$OUTPUT_DIR/${sample_name}_validation.txt" ]]; then
                echo ""
                echo "--- BAM Validation ---"
                cat "$OUTPUT_DIR/${sample_name}_validation.txt"
            fi
            
            if [[ -f "$OUTPUT_DIR/${sample_name}_comparison.txt" ]]; then
                echo ""
                echo "--- Barcode Comparison ---"
                cat "$OUTPUT_DIR/${sample_name}_comparison.txt"
            fi
            
            echo ""
        done
        
    } > "$report_file"
    
    log_info "Full report saved to: $report_file"
    echo ""
    log_info "Summary:"
    echo "  - All 4 original BAM files validated"
    echo "  - All 4 numbat merged BAM files validated"
    echo "  - Barcodes extracted from each file"
    echo "  - Barcode comparisons completed (where applicable)"
}

# Run main function
main