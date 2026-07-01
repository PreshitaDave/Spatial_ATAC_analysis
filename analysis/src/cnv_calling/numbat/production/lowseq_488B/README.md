# NUMBAT Analysis: lowseq_488B

**Tissue**: 488B tissue section  
**Sequencing**: Low-sequenced ATAC (~8M reads)  
**Cells**: ~12,000 barcodes (11,612 post-edge-filter)

## Scripts in This Folder

1. **prepare_numbat_inputs_lowseq_488B.qsub.sh**
   - Prepares inputs for NUMBAT analysis
   - Runtime: ~1-2 hours
   - Output: `Data/04_analysis/cnv/numbat/results/lowseq_488B/`

2. **run_numbat_analysis_lowseq_488B.qsub.sh** (Standard workflow)
   - Standard NUMBAT CNV analysis
   - Use this if pileup has sufficient variant coverage (>100K variants)
   - Runtime: ~1-2 hours

3. **run_numbat_atac_bin_lowseq_488B.qsub.sh** (Recommended for lowseq)
   - ATAC-bin mode: Uses ATAC-seq features instead of SNP coverage
   - Better for low-coverage data with few variants
   - Use this if standard pileup has sparse SNP calls
   - Runtime: ~1-2 hours

## Quick Start

### Option A: ATAC-bin Mode (Recommended)

```bash
# Step 1: Prepare
qsub prepare_numbat_inputs_lowseq_488B.qsub.sh

# Step 2: Run ATAC-bin analysis
qsub run_numbat_atac_bin_lowseq_488B.qsub.sh
```

### Option B: Standard Mode

```bash
# Step 1: Prepare
qsub prepare_numbat_inputs_lowseq_488B.qsub.sh

# Step 2: Run standard analysis
qsub run_numbat_analysis_lowseq_488B.qsub.sh
```

## Which Workflow to Use?

**Use ATAC-bin mode if**:
- Pileup creates <100K variants
- SNP coverage is sparse across cells
- You want to leverage ATAC-seq signal for CNV calling

**Use standard mode if**:
- Pileup creates >1M variants
- SNP coverage is good across cells
- You trust the allele count data

## Key Files

### Inputs
- BAM: `Data/01_inputs/bam/lowseq_488B.bam.lnk`
- Barcodes: `Data/01_inputs/barcodes/tissue_barcodes/lowseq_488B/lowseq_488B.barcodes.tsv`
- Fragments: `Data/01_inputs/fragments/lowseq_488B/`
- Reference: `Data/04_analysis/cnv/numbat/reference/lambdas_*.rds`

### Outputs
- Standard: `Data/04_analysis/cnv/numbat/results/lowseq_488B/numbat_seurat_obj.RDS`
- ATAC-bin: `Data/04_analysis/cnv/numbat/results/lowseq_488B_atac_bin/numbat_seurat_obj.RDS`

## Troubleshooting

### Check pileup coverage
```bash
wc -l Data/04_analysis/cnv/numbat/results/lowseq_488B/alleles.csv
```
- >1M lines → use standard mode
- <100K lines → use ATAC-bin mode

### Validate inputs
```bash
cd ../..
Rscript lib/validate_numbat_inputs.R lowseq_488B
```

## Notes

- Low-seq is challenging for SNP-based CNV detection (sparse variants)
- ATAC-bin mode is an alternative that leverages chromatin accessibility patterns
- If both workflows fail, check: (1) BAM symlink validity, (2) barcode matching, (3) reference files

## See Also
- `../ORGANIZATION.md` - Full folder structure
- `../../lib/get_binned_atac_fixed.R` - ATAC binning (handles barcode mismatches)
- `../lowseq_489/README.md` - Similar workflow for 489
