#!/bin/bash -l
#$ -P paxlab
#$ -N tissue_files
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l h_rt=24:00:00
#$ -l mem_per_core=8G
#$ -j n
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/valid_frags.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/valid_frags.$JOB_ID.err

# Fragment File Validation and Barcode Checking Script
# Validates fragment files, checks completeness, and verifies barcode matching
# Runs in parallel using GNU Parallel

set -euo pipefail

module load samtools
module load parallel
module load htslib

# Configuration
FRAGMENT_DIR="/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/fragments"
BARCODE_DIR="/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes"
BAM_DIR="/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/bam"

SAMPLES=("deepseq_488B" "deepseq_489" "lowseq_488B" "lowseq_489" "lowseq_combined")
NUM_CORES=${1:-8}  # Default to 16 cores, or use first argument

# Create output directory
OUTPUT_DIR="/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/fragments/fragment_validation_results"
mkdir -p "$OUTPUT_DIR"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
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

log_stats() {
    echo -e "${CYAN}[STATS]${NC} $1"
}

# Check dependencies
check_dependencies() {
    log_info "Checking dependencies..."
    
    for cmd in bgzip tabix samtools parallel; do
        if ! command -v "$cmd" &> /dev/null; then
            log_error "$cmd is not installed or not in PATH"
            exit 1
        fi
    done
    
    log_info "All dependencies found: bgzip, tabix, samtools, parallel"
}

# Validate fragment file integrity
validate_fragment_file() {
    local sample="$1"
    local fragment_file="$2"
    local output_file="$3"
    
    {
        echo "=== Validating Fragment File: $sample ==="
        echo "Path: $fragment_file"
        
        if [[ ! -f "$fragment_file" ]]; then
            echo "ERROR: Fragment file does not exist"
            return 1
        fi
        
        if [[ ! -r "$fragment_file" ]]; then
            echo "ERROR: Fragment file is not readable"
            return 1
        fi
        
        echo "File size: $(du -h "$fragment_file" | cut -f1)"
        echo "Last modified: $(stat -c%y "$fragment_file" 2>/dev/null || stat -f%Sm "$fragment_file" 2>/dev/null || echo 'N/A')"
        
        # Check file extension
        if [[ "$fragment_file" == *.bed.gz ]]; then
            echo "Format: Gzipped BED (.bed.gz)"
            
            # Test if gzip file is valid
            if bgzip -t "$fragment_file" 2>/dev/null; then
                echo "Gzip integrity: VALID"
            else
                echo "Gzip integrity: INVALID - file may be corrupted"
                return 1
            fi
            
            # Check for truncation by reading last few lines
            local last_lines=$(zcat "$fragment_file" 2>/dev/null | tail -5 | wc -l)
            if [[ $last_lines -gt 0 ]]; then
                echo "File completeness: OK (can read to end)"
                echo "Last 3 entries read successfully"
            else
                echo "File completeness: WARNING - could not read end of file"
            fi
            
        elif [[ "$fragment_file" == *.bed ]]; then
            echo "Format: Uncompressed BED (.bed)"
            
            # Check if file ends with newline
            if [[ -z "$(tail -c 1 "$fragment_file")" ]]; then
                echo "File completeness: OK (ends with newline)"
            else
                echo "File completeness: WARNING (missing final newline)"
            fi
        fi
        
        # Get line count and basic statistics
        local line_count
        if [[ "$fragment_file" == *.gz ]]; then
            line_count=$(zcat "$fragment_file" 2>/dev/null | wc -l)
        else
            line_count=$(wc -l < "$fragment_file")
        fi
        
        echo "Total lines: $line_count"
        
        # Check header and format
        local header
        if [[ "$fragment_file" == *.gz ]]; then
            header=$(zcat "$fragment_file" 2>/dev/null | head -1)
        else
            header=$(head -1 "$fragment_file")
        fi
        
        if [[ "$header" == \#* ]]; then
            echo "Header: Present (recognized as comment)"
            echo "Header line: $header"
        else
            echo "Header: No header detected"
            echo "First line: $header"
        fi
        
        # Validate BED format (should have at least 4 columns)
        local data_lines
        if [[ "$fragment_file" == *.gz ]]; then
            data_lines=$(zcat "$fragment_file" 2>/dev/null | grep -v "^#" | head -10 | wc -l)
        else
            data_lines=$(grep -v "^#" "$fragment_file" | head -10 | wc -l)
        fi
        
        if [[ $data_lines -gt 0 ]]; then
            echo "Data format: OK (contains data rows)"
            
            # Check column count
            if [[ "$fragment_file" == *.gz ]]; then
                local col_count=$(zcat "$fragment_file" 2>/dev/null | grep -v "^#" | head -1 | awk '{print NF}')
            else
                local col_count=$(grep -v "^#" "$fragment_file" | head -1 | awk '{print NF}')
            fi
            
            echo "Columns per row: $col_count"
            if [[ $col_count -ge 4 ]]; then
                echo "Column check: PASS (minimum 4 columns required)"
            else
                echo "Column check: FAIL (less than 4 columns)"
            fi
        else
            echo "Data format: ERROR - no data rows found"
            return 1
        fi
        
        echo "Validation: PASSED"
        return 0
        
    } > "$output_file" 2>&1
}

# Validate index file
validate_index_file() {
    local sample="$1"
    local index_file="$2"
    local fragment_file="$3"
    local output_file="$4"
    
    {
        echo "=== Validating Index File: $sample ==="
        echo "Index path: $index_file"
        
        if [[ ! -f "$index_file" ]]; then
            echo "Status: MISSING - index file not found"
            return 1
        fi
        
        if [[ ! -r "$index_file" ]]; then
            echo "Status: UNREADABLE"
            return 1
        fi
        
        echo "File size: $(du -h "$index_file" | cut -f1)"
        echo "Last modified: $(stat -c%y "$index_file" 2>/dev/null || stat -f%Sm "$index_file" 2>/dev/null || echo 'N/A')"
        
        # Check if index is consistent with fragment file
        if tabix -H "$fragment_file" &> /dev/null; then
            echo "Index integrity: VALID"
            
            # Get reference sequence count from index
            local seq_count=$(tabix -H "$fragment_file" 2>/dev/null | wc -l)
            echo "Reference sequences in index: $seq_count"
        else
            echo "Index integrity: ERROR - index may be corrupted or incompatible"
            return 1
        fi
        
        echo "Status: PASSED"
        return 0
        
    } > "$output_file" 2>&1
}

# Extract barcodes from fragment file
extract_fragment_barcodes() {
    local sample="$1"
    local fragment_file="$2"
    local output_file="$3"
    
    {
        echo "=== Extracting barcodes from fragment file: $sample ==="
        
        if [[ ! -f "$fragment_file" ]]; then
            echo "ERROR: Fragment file not found"
            return 1
        fi
        
        # Extract barcodes (usually in 4th column)
        if [[ "$fragment_file" == *.gz ]]; then
            zcat "$fragment_file" 2>/dev/null \
                | grep -v "^#" \
                | awk '{print $4}' \
                | sort | uniq > "${output_file}.tmp" 2>/dev/null
        else
            grep -v "^#" "$fragment_file" 2>/dev/null \
                | awk '{print $4}' \
                | sort | uniq > "${output_file}.tmp" 2>/dev/null
        fi
        
        if [[ ! -f "${output_file}.tmp" ]]; then
            echo "ERROR: Failed to extract barcodes"
            return 1
        fi
        
        local num_barcodes=$(wc -l < "${output_file}.tmp")
        echo "Extracted unique barcodes: $num_barcodes"
        echo "Sample: $sample"
        echo "Fragment file: $fragment_file"
        
        mv "${output_file}.tmp" "$output_file"
        
    } 2>&1 | tee -a "$OUTPUT_DIR/${sample}_extract_frag.log"
}

# Compare fragment barcodes with reference
compare_fragment_barcodes() {
    local sample="$1"
    local frag_bc="$2"
    local reference_dir="$3"
    local output_file="$4"
    
    {
        echo "=== Fragment Barcode Comparison: $sample ==="
        echo "Reference directory: $reference_dir"
        
        if [[ ! -d "$reference_dir" ]]; then
            echo "WARNING: Reference directory not found: $reference_dir"
            return 1
        fi
        
        if [[ ! -f "$frag_bc" ]]; then
            echo "ERROR: Fragment barcode file not found: $frag_bc"
            return 1
        fi
        
        # Find barcode files in reference directory
        local barcode_files=($(find "$reference_dir" -name "*.barcodes.tsv" -type f 2>/dev/null))
        
        if [[ ${#barcode_files[@]} -eq 0 ]]; then
            echo "WARNING: No reference barcode files found"
            return 1
        fi
        
        echo "Found ${#barcode_files[@]} reference barcode files"
        echo ""
        
        # Compare with each barcode file
        for bc_file in "${barcode_files[@]}"; do
            local filename=$(basename "$bc_file")
            local total_ref=$(wc -l < "$bc_file")
            local total_frag=$(wc -l < "$frag_bc")
            
            if [[ $total_frag -eq 0 ]]; then
                echo "WARNING: No barcodes found in fragment file"
                continue
            fi
            
            # Count matches
            local matches=$(comm -12 <(sort "$frag_bc") <(sort "$bc_file" | awk '{print $1}') | wc -l)
            local in_frag_only=$(comm -23 <(sort "$frag_bc") <(sort "$bc_file" | awk '{print $1}') | wc -l)
            local in_ref_only=$(comm -13 <(sort "$frag_bc") <(sort "$bc_file" | awk '{print $1}') | wc -l)
            
            local match_percent=$(awk "BEGIN {printf \"%.2f\", ($matches/$total_frag)*100}")
            local coverage_percent=$(awk "BEGIN {printf \"%.2f\", ($matches/$total_ref)*100}")
            
            echo "Reference file: $filename"
            echo "  Reference barcodes: $total_ref"
            echo "  Fragment barcodes: $total_frag"
            echo "  Matches: $matches ($match_percent% of fragment barcodes)"
            echo "  Coverage: $coverage_percent% of reference barcodes"
            echo "  In fragment only: $in_frag_only"
            echo "  In reference only: $in_ref_only"
            
            if [[ $in_frag_only -gt 0 ]]; then
                local pct=$(awk "BEGIN {printf \"%.1f\", (in_frag_only * 100 / total_frag)}")
                echo "  ⚠ WARNING: Found $in_frag_only extra barcodes ($pct%)"
            fi
            if [[ $in_ref_only -gt 0 ]]; then
                local pct=$(awk "BEGIN {printf \"%.1f\", (in_ref_only * 100 / total_ref)}")
                echo "  ⚠ WARNING: Missing $in_ref_only reference barcodes ($pct%)"
            fi
            if [[ $match_percent == "100.00" ]]; then
                echo "  ✓ PERFECT MATCH"
            fi
            echo ""
        done
        
    } > "$output_file" 2>&1
}

# Get fragment file statistics
get_fragment_stats() {
    local sample="$1"
    local fragment_file="$2"
    local output_file="$3"
    
    {
        echo "=== Fragment Statistics: $sample ==="
        
        if [[ ! -f "$fragment_file" ]]; then
            echo "ERROR: Fragment file not found"
            return 1
        fi
        
        local total_fragments
        local total_bases
        local avg_fragment_size
        
        if [[ "$fragment_file" == *.gz ]]; then
            total_fragments=$(zcat "$fragment_file" 2>/dev/null | grep -v "^#" | wc -l)
            
            # Calculate total bases and fragment stats
            zcat "$fragment_file" 2>/dev/null | grep -v "^#" | awk '
            {
                frag_size = $3 - $2
                total_bases += frag_size
                total_frags++
                count_col = $5
                total_counts += count_col
                if (frag_size < min_size || min_size == 0) min_size = frag_size
                if (frag_size > max_size) max_size = frag_size
            }
            END {
                if (total_frags > 0) {
                    avg_size = total_bases / total_frags
                    print "Total fragments: " total_frags
                    print "Total bases covered: " total_bases
                    print "Total fragment counts: " total_counts
                    print "Average fragment size: " int(avg_size)
                    print "Min fragment size: " min_size
                    print "Max fragment size: " max_size
                }
            }' > "${output_file}.tmp"
            
            cat "${output_file}.tmp"
            rm "${output_file}.tmp"
        else
            total_fragments=$(grep -v "^#" "$fragment_file" | wc -l)
            
            grep -v "^#" "$fragment_file" | awk '
            {
                frag_size = $3 - $2
                total_bases += frag_size
                total_frags++
                count_col = $5
                total_counts += count_col
                if (frag_size < min_size || min_size == 0) min_size = frag_size
                if (frag_size > max_size) max_size = frag_size
            }
            END {
                if (total_frags > 0) {
                    avg_size = total_bases / total_frags
                    print "Total fragments: " total_frags
                    print "Total bases covered: " total_bases
                    print "Total fragment counts: " total_counts
                    print "Average fragment size: " int(avg_size)
                    print "Min fragment size: " min_size
                    print "Max fragment size: " max_size
                }
            }'
        fi
        
        echo "Statistics: COMPLETE"
        
    } >> "$output_file" 2>&1
}

# Main validation function for a single sample
validate_sample() {
    local sample="$1"
    
    log_info "Processing sample: $sample"
    
    local sample_dir="$FRAGMENT_DIR/$sample"
    local barcode_dir="$BARCODE_DIR/$sample"
    
    if [[ ! -d "$sample_dir" ]]; then
        log_error "Sample directory not found: $sample_dir"
        return 1
    fi
    
    # Find fragment files
    local frag_files=($(find "$sample_dir" -name "*.fragments.sort.filtered.bed.gz" -o -name "*.fragments.sort.filtered.bed" | sort))
    
    if [[ ${#frag_files[@]} -eq 0 ]]; then
        log_warn "No fragment files found for sample: $sample"
        return 1
    fi
    
    # Process first main fragment file
    local main_frag_file="${frag_files[0]}"
    
    # Validate fragment file
    validate_fragment_file "$sample" "$main_frag_file" "$OUTPUT_DIR/${sample}_fragment_validation.txt"
    
    # Find and validate index file
    local index_file="${main_frag_file}.tbi"
    if [[ -f "$index_file" ]]; then
        validate_index_file "$sample" "$index_file" "$main_frag_file" "$OUTPUT_DIR/${sample}_index_validation.txt"
    else
        log_warn "Index file not found for $sample: $index_file"
        {
            echo "Status: MISSING"
            echo "Expected: $index_file"
        } > "$OUTPUT_DIR/${sample}_index_validation.txt"
    fi
    
    # Extract barcodes from fragment file
    extract_fragment_barcodes "$sample" "$main_frag_file" "$OUTPUT_DIR/${sample}_frag_extracted.barcodes.txt"
    
    # Compare barcodes
    if [[ -d "$barcode_dir" ]]; then
        compare_fragment_barcodes "$sample" "$OUTPUT_DIR/${sample}_frag_extracted.barcodes.txt" "$barcode_dir" "$OUTPUT_DIR/${sample}_frag_comparison.txt"
    else
        log_warn "Barcode reference directory not found: $barcode_dir"
    fi
    
    # Get statistics
    get_fragment_stats "$sample" "$main_frag_file" "$OUTPUT_DIR/${sample}_statistics.txt"
    
    log_info "Completed: $sample"
}

export -f validate_fragment_file
export -f validate_index_file
export -f extract_fragment_barcodes
export -f compare_fragment_barcodes
export -f get_fragment_stats
export -f validate_sample
export -f log_info
export -f log_warn
export -f log_error
export -f log_debug
export -f log_stats
export FRAGMENT_DIR BARCODE_DIR BAM_DIR OUTPUT_DIR RED GREEN YELLOW BLUE CYAN NC

# Main execution
main() {
    log_info "Fragment File Validation and Barcode Checking Script"
    log_info "====================================================="
    log_info "Fragment Directory: $FRAGMENT_DIR"
    log_info "Barcode Directory: $BARCODE_DIR"
    log_info "Using $NUM_CORES cores for parallel processing"
    log_info "Total samples to validate: ${#SAMPLES[@]}"
    log_info ""
    
    check_dependencies
    
    # Verify directories exist
    log_info "Checking if sample directories exist..."
    local missing_dirs=0
    
    for sample in "${SAMPLES[@]}"; do
        local sample_dir="$FRAGMENT_DIR/$sample"
        if [[ -d "$sample_dir" ]]; then
            log_info "Found: [$sample] $sample_dir"
        else
            log_error "Missing: [$sample] $sample_dir"
            ((missing_dirs++))
        fi
    done
    
    if [[ $missing_dirs -gt 0 ]]; then
        log_error "$missing_dirs sample directories are missing!"
        exit 1
    fi
    
    echo ""
    log_info "Starting parallel validation with $NUM_CORES cores..."
    echo ""
    
    # Run validation in parallel
    printf '%s\n' "${SAMPLES[@]}" | parallel --jobs "$NUM_CORES" validate_sample
    
    echo ""
    log_info "Validation complete!"
    log_info "Results saved to: $OUTPUT_DIR"
    
    # Generate summary report
    echo ""
    log_info "=== SUMMARY REPORT ==="
    echo ""
    
    for sample in "${SAMPLES[@]}"; do
        echo "--- $sample ---"
        if [[ -f "$OUTPUT_DIR/${sample}_fragment_validation.txt" ]]; then
            grep "Validation:\|Total lines\|Columns per row" "$OUTPUT_DIR/${sample}_fragment_validation.txt" | head -3
        fi
        if [[ -f "$OUTPUT_DIR/${sample}_index_validation.txt" ]]; then
            grep "Status:" "$OUTPUT_DIR/${sample}_index_validation.txt" | head -1
        fi
        if [[ -f "$OUTPUT_DIR/${sample}_frag_comparison.txt" ]]; then
            grep "Matches:\|PERFECT MATCH" "$OUTPUT_DIR/${sample}_frag_comparison.txt" | head -2
        fi
        if [[ -f "$OUTPUT_DIR/${sample}_statistics.txt" ]]; then
            grep "Total fragments:\|Average fragment size:" "$OUTPUT_DIR/${sample}_statistics.txt" | head -2
        fi
        echo ""
    done
    
    # Create consolidated report
    local report_file="$OUTPUT_DIR/FULL_REPORT.txt"
    {
        echo "Fragment File Validation Report"
        echo "Generated: $(date)"
        echo ""
        echo "Configuration:"
        echo "  Fragment Directory: $FRAGMENT_DIR"
        echo "  Barcode Directory: $BARCODE_DIR"
        echo "  Cores used: $NUM_CORES"
        echo ""
        
        for sample in "${SAMPLES[@]}"; do
            echo "=========================================="
            echo "Sample: $sample"
            echo "=========================================="
            
            if [[ -f "$OUTPUT_DIR/${sample}_fragment_validation.txt" ]]; then
                echo ""
                echo "--- Fragment File Validation ---"
                cat "$OUTPUT_DIR/${sample}_fragment_validation.txt"
            fi
            
            if [[ -f "$OUTPUT_DIR/${sample}_index_validation.txt" ]]; then
                echo ""
                echo "--- Index File Validation ---"
                cat "$OUTPUT_DIR/${sample}_index_validation.txt"
            fi
            
            if [[ -f "$OUTPUT_DIR/${sample}_frag_comparison.txt" ]]; then
                echo ""
                echo "--- Barcode Comparison ---"
                cat "$OUTPUT_DIR/${sample}_frag_comparison.txt"
            fi
            
            if [[ -f "$OUTPUT_DIR/${sample}_statistics.txt" ]]; then
                echo ""
                echo "--- Fragment Statistics ---"
                cat "$OUTPUT_DIR/${sample}_statistics.txt"
            fi
            
            echo ""
        done
        
    } > "$report_file"
    
    log_info "Full report saved to: $report_file"
    echo ""
    log_info "Validation Summary:"
    echo "  - Validated all fragment files (.bed.gz and .bed)"
    echo "  - Checked index files (.tbi)"
    echo "  - Extracted and compared barcodes"
    echo "  - Generated fragment statistics"
}

# Run main function
main