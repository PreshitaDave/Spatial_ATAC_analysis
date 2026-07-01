# NUMBAT Analysis Folder Organization

**Last Updated**: June 8, 2026

## Folder Structure

```
numbat/
├── production/              # ← PRODUCTION SCRIPTS (use these!)
│   ├── deepseq_488B/       # Tissue-specific scripts
│   ├── deepseq_489/
│   ├── lowseq_488B/
│   └── lowseq_489/
├── lib/                     # Reusable utilities and helpers
├── archive/                 # Archived files (old/test scripts)
│   ├── test_scripts/        # Test/debug scripts (removed from production)
│   └── old_versions/        # Deprecated/alternative versions
├── numbat/                  # Working scripts directory (subdirectory with shared utilities)
├── README.md                # Original documentation
└── REGENERATE_ATAC_WORKFLOW.md  # ATAC reference generation guide
```

## Quick Start

### For Deepseq Analysis

```bash
cd production/deepseq_488B
# 1. Prepare inputs
qsub prepare_numbat_inputs_deepseq_488B.qsub.sh

# 2. Run analysis
qsub run_numbat_analysis_deepseq_488B.qsub.sh
```

### For Lowseq Analysis (ATAC-bin mode)

```bash
cd production/lowseq_488B
# 1. Prepare inputs
qsub prepare_numbat_inputs_lowseq_488B.qsub.sh

# 2. Run pileup/phasing (if needed)
# Reference: lib/numbat/run_numbat_atac_bin.sh

# 3. Run ATAC-bin analysis
qsub run_numbat_atac_bin_lowseq_488B.qsub.sh

# 4. OR: Run standard analysis
qsub run_numbat_analysis_lowseq_488B.qsub.sh
```

## Production Scripts

### Tissue-Specific Folders

Each tissue folder contains THREE types of scripts:

1. **prepare_numbat_inputs_{tissue}.qsub.sh**
   - Prepares allele counts and variant calls (pileup stage)
   - Generates alleles.csv for CNV calling
   - Must run BEFORE analysis

2. **run_numbat_analysis_{tissue}.qsub.sh**
   - Main NUMBAT CNV analysis
   - Calls run_numbat_multiome.R with tissue-specific parameters
   - Produces CNV profiles and tumor/normal classification

3. **run_numbat_atac_bin_{tissue}.qsub.sh** (lowseq only)
   - Alternative workflow for ATAC-bin mode
   - Uses ATAC-seq-based binning for CNV detection
   - Better for lowseq data with limited SNP coverage

### Supported Tissues

- `deepseq_488B` - Deep-sequenced 488B tissue (~12K cells)
- `deepseq_489` - Deep-sequenced 489 tissue (~5K cells)  
- `lowseq_488B` - Low-sequenced 488B tissue (~12K cells)
- `lowseq_489` - Low-sequenced 489 tissue (~5K cells)

## Library Scripts (lib/)

Reusable utilities and core functions:

### Data Preparation
- `extract_barcodes_from_bam.sh` - Extract cell barcodes from BAM
- `merge_bam_for_numbat.qsub.sh` - Merge BAM files across regions
- `merge_lowseq_fragments.qsub.sh` - Merge fragment files

### Analysis Utilities
- `get_binned_atac_fixed.R` - Generate ATAC-bin matrices (FIXED for barcode mismatch)
- `validate_numbat_inputs.R` - Validate input consistency
- `run_numbat_multiome.R` - Main NUMBAT analysis R script
- `postprocess_numbat_results.R` - Post-process CNV calls

### Reference Generation
- `run_reference_generation.sh` - Generate NUMBAT reference matrices
- `setup_numbat_references.sh` - Initialize reference files
- `regenerate_atac_full_barcodes.qsub.sh` - Regenerate ATAC matrices with correct barcodes

### Orchestration & Templates
- `submit_all_tissue_numbat_orchestrated.sh` - Run all tissues in sequence
- `submit_deepseq_prepare_sequential.sh` - Prepare all deepseq tissues
- `run_numbat_analysis_with_validation_TEMPLATE.qsub.sh` - Template with validation

## Archive

### test_scripts/ (9 scripts)
**Removed**: Test/debug scripts and preliminary experiments
- `check_params.qsub.sh`
- `test_deepseq*.sh`
- `generate_*_references.test.sh`
- `test_run_numbat*.qsub.sh`, `test_run_numbat*.R`

### old_versions/ (23 scripts)
**Removed**: Deprecated/alternative workflows and experiments
- `*_refhca*.qsub.sh` - Reference panel variations (mostly experimental)
- `*_v2.qsub.sh` - Older versions of tissue scripts
- `*_combined.qsub.sh` - Combined tissue analyses (deprecated)
- `*_multiome_*.qsub.sh` - Early multiome integration attempts
- `run_numbat_*_analysis_only.qsub.sh` - Analysis-only variants

## Critical Paths & Files

### Input Locations
```
Data/01_inputs/barcodes/tissue_barcodes/{tissue}/
  - {tissue}.barcodes.tsv              # Full barcode set
  - {tissue}.no_edge_effect.barcodes.tsv  # Filtered barcodes

Data/01_inputs/bam/
  - *.bam.lnk                          # Symlinks to BAM files

Data/01_inputs/fragments/{tissue}/
  - *.bed.gz                           # Fragment files
```

### Reference Files
```
Data/04_analysis/cnv/numbat/reference/
  - lambdas_ATAC_bincnt.rds            # ATAC aggregated reference
  - var220kb.rds                       # Genomic bins (220kb)
  - par_numbatm.rds                    # NUMBAT parameters
  - phased_panel_bcf_links/            # Eagle phasing panel
```

### Output Locations
```
Data/04_analysis/cnv/numbat/results/{tissue}/
  - alleles.csv                        # From pileup stage
  - adata_atac.rds                     # ATAC bin matrix
  - numbat_seurat_obj.RDS              # CNV calls & tumor profiles
```

## Best Practices

### 1. Always Validate Before Running
```bash
# Check input consistency
Rscript lib/validate_numbat_inputs.R {tissue}
```

### 2. Use Correct Barcode Files
- Pileup stage: Full barcode set (`*.barcodes.tsv`)
- ATAC binning: Must use SAME barcode file as pileup
- Mismatch causes 0-coverage cells → no CNV calls

### 3. Monitor Jobs
```bash
# Check status
qstat -j {JOBID}

# Watch logs in real-time
tail -f analysis/qsub_logs/{name}_{JOBID}.log
```

### 4. Reference Generation (First Time Setup)
```bash
# In interactive session
qrsh -l h_rt=02:00:00 -pe omp 8 -P paxlab
cd analysis/src/numbat
bash lib/run_reference_generation.sh
```

## Reorganization Summary

**Date**: June 8, 2026

### What Changed

1. **Cleaned up root directory**
   - Removed 32 test/debug scripts → `archive/test_scripts/`
   - Removed 23 deprecated scripts → `archive/old_versions/`
   - Kept only 2 documentation files at root

2. **Organized by tissue**
   - Created `production/{tissue}/` folders
   - Each tissue has 2-3 focused scripts
   - Clearer production workflow

3. **Centralized utilities**
   - All helper scripts in `lib/`
   - All R analysis functions in `lib/`
   - Templates and orchestration scripts in `lib/`

4. **Improved discoverability**
   - Clear separation: production vs. utilities vs. archived
   - Self-documenting folder structure
   - Easier to find the right script

### Files Removed (Now in Archive)

- **9 test scripts**: Various debug and test runs
- **23 old versions**: Experimental workflows, parallel attempts, deprecated solutions

See `archive/` for complete history.

## Troubleshooting

### Scripts Not Found
- Check: Are you in `analysis/src/numbat/production/{tissue}/`?
- All production scripts are tissue-specific and located there

### Input Path Errors
- Verify: `Data/01_inputs/fragments/{tissue}/`, `Data/01_inputs/barcodes/tissue_barcodes/{tissue}/`
- Use absolute paths in scripts

### Barcode Mismatch Issues
- Run: `Rscript lib/validate_numbat_inputs.R {tissue}`
- Fix: `Rscript lib/get_binned_atac_fixed.R` (uses corrected barcode handling)

### Reference Not Found
- Regenerate: `qrsh` → `bash lib/run_reference_generation.sh`
- Check: `Data/04_analysis/cnv/numbat/reference/` for lambdas_*.rds files

## Questions?

Refer to:
- `README.md` - Original workflow documentation
- `REGENERATE_ATAC_WORKFLOW.md` - Reference generation guide
- `lib/run_numbat_analysis_with_validation_TEMPLATE.qsub.sh` - Template with inline documentation
