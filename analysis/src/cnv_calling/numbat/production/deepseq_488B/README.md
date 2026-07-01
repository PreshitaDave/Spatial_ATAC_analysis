# NUMBAT Analysis: deepseq_488B

**Tissue**: 488B tissue section  
**Sequencing**: Deep-sequenced ATAC (~190M reads)  
**Cells**: ~12,000 barcodes (11,467 post-edge-filter)

## Scripts in This Folder

1. **prepare_numbat_inputs_deepseq_488B.qsub.sh**
   - Prepares allele counts from BAM file via pileup
   - Creates `alleles.csv` (variant calls per barcode)
   - Runtime: ~2-4 hours
   - Output: `Data/04_analysis/cnv/numbat/results/deepseq_488B/alleles.csv`

2. **run_numbat_analysis_deepseq_488B.qsub.sh**
   - Main NUMBAT CNV analysis
   - Matches allele counts to ATAC binning
   - Infers CNV profiles and tumor classification
   - Runtime: ~1-2 hours
   - Output: `Data/04_analysis/cnv/numbat/results/deepseq_488B/numbat_seurat_obj.RDS`

## Quick Start

```bash
# From this folder:
qsub prepare_numbat_inputs_deepseq_488B.qsub.sh

# After pileup completes (check qstat):
qsub run_numbat_analysis_deepseq_488B.qsub.sh
```

## Key Files

### Inputs (must exist before running)
- BAM: `Data/01_inputs/bam/deepseq_488B.bam.lnk`
- Barcodes: `Data/01_inputs/barcodes/tissue_barcodes/deepseq_488B/deepseq_488B.barcodes.tsv`
- Fragments: `Data/01_inputs/fragments/deepseq_488B/`
- Reference: `Data/04_analysis/cnv/numbat/reference/lambdas_*.rds`

### Outputs (after completion)
- Alleles: `Data/04_analysis/cnv/numbat/results/deepseq_488B/alleles.csv`
- ATAC matrix: `Data/04_analysis/cnv/numbat/results/deepseq_488B/adata_atac.rds`
- CNV calls: `Data/04_analysis/cnv/numbat/results/deepseq_488B/numbat_seurat_obj.RDS`

## Troubleshooting

### Pileup fails with "BAM file not found"
```bash
# Verify BAM symlink
ls -lh Data/01_inputs/bam/deepseq_488B.bam.lnk
file Data/01_inputs/bam/deepseq_488B.bam.lnk  # should show: GZIP compressed data
```

### Analysis fails with "Barcode mismatch"
```bash
# Validate input consistency
cd ../..
Rscript lib/validate_numbat_inputs.R deepseq_488B
```

### CNV calls are sparse or missing
- Check pileup output: `wc -l Data/04_analysis/cnv/numbat/results/deepseq_488B/alleles.csv`
  - Should be >1M lines (variants × cells)
  - If <100K: pileup may have failed silently
- Verify reference files exist: `ls Data/04_analysis/cnv/numbat/reference/`

## Notes

- Deep-seq has high SNP coverage, good for standard NUMBAT workflow
- Consider ATAC-bin mode if standard analysis shows sparse CNV calls
- Run validation script if anything seems off

## See Also
- `../ORGANIZATION.md` - Full folder structure
- `../../lib/validate_numbat_inputs.R` - Input validation
- `../../lib/run_numbat_multiome.R` - Analysis script source
