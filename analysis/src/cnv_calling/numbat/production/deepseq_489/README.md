# NUMBAT Analysis: deepseq_489

**Tissue**: 489 tissue section  
**Sequencing**: Deep-sequenced ATAC (~50M reads)  
**Cells**: ~5,000 barcodes (4,622 post-edge-filter)

## Scripts in This Folder

1. **prepare_numbat_inputs_deepseq_489.qsub.sh**
   - Prepares allele counts from BAM file via pileup
   - Creates `alleles.csv` (variant calls per barcode)
   - Runtime: ~1-2 hours
   - Output: `Data/04_analysis/cnv/numbat/results/deepseq_489/alleles.csv`

2. **run_numbat_analysis_deepseq_489.qsub.sh**
   - Main NUMBAT CNV analysis
   - Runtime: ~30-60 minutes (smaller tissue)
   - Output: `Data/04_analysis/cnv/numbat/results/deepseq_489/numbat_seurat_obj.RDS`

## Quick Start

```bash
qsub prepare_numbat_inputs_deepseq_489.qsub.sh
# Wait for pileup to complete
qsub run_numbat_analysis_deepseq_489.qsub.sh
```

## Key Files

### Inputs
- BAM: `Data/01_inputs/bam/deepseq_489.bam.lnk`
- Barcodes: `Data/01_inputs/barcodes/tissue_barcodes/deepseq_489/deepseq_489.barcodes.tsv`
- Reference: `Data/04_analysis/cnv/numbat/reference/lambdas_*.rds`

### Outputs
- Alleles: `Data/04_analysis/cnv/numbat/results/deepseq_489/alleles.csv`
- CNV: `Data/04_analysis/cnv/numbat/results/deepseq_489/numbat_seurat_obj.RDS`

## See Also
- `../ORGANIZATION.md` - Full documentation
- `../deepseq_488B/README.md` - Similar workflow for 488B
