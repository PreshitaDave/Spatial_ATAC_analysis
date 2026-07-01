# NUMBAT Library - Reusable Utilities

This folder contains helper scripts, utilities, and core analysis functions used by NUMBAT workflows.

## Data Preparation Scripts

### BAM Processing
- **extract_barcodes_from_bam.sh** - Extract cell barcodes from BAM file headers
- **extract_barcodes_from_bam.qsub.sh** - SGE job wrapper for barcode extraction
- **merge_bam_for_numbat.qsub.sh** - Merge BAM files across chromosomes/regions
- **merge_lowseq_fragments.qsub.sh** - Merge fragment files from multiple lowseq runs

### Data Preparation
- **prepare_numbat_inputs.sh** - Generic input preparation (used by tissue-specific scripts)
- **prepare_deepseq_tissue_inputs.sh** - Deepseq-specific input organization
- **prepare_numbat_atac_inputs.sh** - Prepare ATAC-specific inputs

## Analysis Scripts (R)

### Core Analysis
- **run_numbat_multiome.R** - Main NUMBAT CNV analysis engine
  - Input: allele counts + ATAC binning
  - Output: CNV profiles, tumor classification
  - Called by: `run_numbat_analysis_{tissue}.qsub.sh`

- **run_numbat_combined.R** - Combined multi-tissue analysis (legacy)
- **run_numbat_refhca.R** - NUMBAT with HCA reference panel (experimental)

### Data Processing
- **get_binned_atac_fixed.R** ⭐ **IMPORTANT**
  - Generates ATAC bin matrices for CNV detection
  - Handles barcode suffix mismatches (required for lowseq)
  - Output: adata_atac.rds
  - Use this instead of standard binning scripts

- **validate_numbat_inputs.R** ⭐ **CRITICAL**
  - Validates input file consistency
  - Checks: barcode overlap, cell counts, coverage
  - Run this before analysis if anything seems off
  - Usage: `Rscript validate_numbat_inputs.R {tissue}`

- **postprocess_numbat_results.R** - Post-processing CNV calls
  - Optional: enhances CNV output with annotations

- **combine_lowseq_archr_projects.R** - Combine ArchR projects
  - Integrates with ArchR analysis pipeline (lowseq-specific)

### Testing/Debug
- **test_numbat_load.R** - Test NUMBAT loading and basic functionality
- **test_numbat_minimal.R** - Minimal NUMBAT execution for debugging

## Reference Management

- **run_reference_generation.sh** - Generate tissue-specific NUMBAT references
  - Creates: lambdas_ATAC_bincnt.rds
  - Runtime: ~30 minutes per tissue
  - Must run BEFORE first analysis
  - Usage: `qrsh` → `bash lib/run_reference_generation.sh`

- **setup_numbat_references.sh** - Initialize reference file structure
  - Sets up directory: `Data/04_analysis/cnv/numbat/reference/`
  - One-time setup script

- **regenerate_atac_full_barcodes.qsub.sh** - Regenerate ATAC matrices
  - Use if barcode file changes or mismatch detected
  - Alternative: `regenerate_atac_with_correct_barcodes.qsub.sh`

## Orchestration & Templates

- **submit_all_tissue_numbat_orchestrated.sh** - Run all 4 tissues in sequence
  - Submits: prepare → analyze for each tissue
  - Handles job dependencies
  - Usage: `bash lib/submit_all_tissue_numbat_orchestrated.sh`

- **submit_deepseq_prepare_sequential.sh** - Prepare deepseq tissues only
  - Prepares 488B and 489 in sequence

- **run_numbat_analysis_with_validation_TEMPLATE.qsub.sh** ⭐ **TEMPLATE**
  - Template for creating new tissue-specific analysis scripts
  - Includes: path validation, module loading, error handling
  - Copy and customize for new tissues/parameters

## Common Utilities

- **numbat_common.sh** - Bash utility functions
  - Shared across all scripts
  - Contains: path checks, logging, error handling

- **setup_numbat_tools.sh** - Initialize tools and environment
  - Loads R modules
  - Configures environment variables

## Usage Patterns

### Standard Workflow
```bash
# 1. Check inputs
Rscript lib/validate_numbat_inputs.R deepseq_488B

# 2. Prepare
qsub production/deepseq_488B/prepare_numbat_inputs_deepseq_488B.qsub.sh

# 3. Analyze
qsub production/deepseq_488B/run_numbat_analysis_deepseq_488B.qsub.sh
```

### Troubleshooting
```bash
# 1. Validate inputs first
Rscript lib/validate_numbat_inputs.R {tissue}

# 2. Check ATAC binning
Rscript lib/get_binned_atac_fixed.R {tissue}

# 3. Review template
less lib/run_numbat_analysis_with_validation_TEMPLATE.qsub.sh
```

### Reference Generation (First Time)
```bash
qrsh -l h_rt=02:00:00 -pe omp 8 -P paxlab
cd analysis/src/numbat
bash lib/run_reference_generation.sh
```

### All Tissues At Once
```bash
bash lib/submit_all_tissue_numbat_orchestrated.sh
```

## Key Functions in Scripts

### get_binned_atac_fixed.R
```R
# Main function: Generate ATAC bin matrix
# - Strips "-1"/"-2" suffix from fragment barcodes
# - Matches to cell barcodes without suffix
# - Fixes: "No fragments matched" error
```

### validate_numbat_inputs.R
```R
# Main function: Check consistency
# - Compares ATAC and allele cell counts
# - Identifies barcode overlap
# - Reports pass/fail verdict
```

### run_numbat_multiome.R
```R
# Main function: CNV analysis
# Input: alleles.csv, adata_atac.rds
# Output: Seurat object with CNV calls
```

## File Locations

### Inputs (read from)
- `Data/01_inputs/bam/*.bam.lnk` - BAM symlinks
- `Data/01_inputs/fragments/{tissue}/` - Fragment files
- `Data/01_inputs/barcodes/tissue_barcodes/{tissue}/` - Barcode files

### References (read from)
- `Data/04_analysis/cnv/numbat/reference/` - Generated references

### Outputs (written to)
- `Data/04_analysis/cnv/numbat/results/{tissue}/` - Analysis outputs

## Notes

⭐ **IMPORTANT**: 
- Always use `get_binned_atac_fixed.R` (not original get_binned_atac.R)
- Run `validate_numbat_inputs.R` if anything seems wrong
- Reference files must be generated before first analysis

📋 **TEMPLATE**:
- Use `run_numbat_analysis_with_validation_TEMPLATE.qsub.sh` as starting point
- Copy, customize, test, then move to production

See Also:
- `../ORGANIZATION.md` - Full folder structure
- `../production/*/README.md` - Tissue-specific workflows
- `../archive/README.md` - Removed scripts and history
