# NUMBAT Analysis: lowseq_489

**Tissue**: 489 tissue section  
**Sequencing**: Low-sequenced ATAC (~8M reads)  
**Cells**: ~5,000 barcodes (4,622 post-edge-filter)

## Scripts in This Folder

1. **prepare_numbat_inputs_lowseq_489.qsub.sh**
   - Prepares inputs (allele counts, ATAC binning)
   - Runtime: ~30-60 minutes (smaller tissue)

2. **run_numbat_analysis_lowseq_489.qsub.sh** (Standard)
   - Standard SNP-based CNV analysis

3. **run_numbat_atac_bin_lowseq_489.qsub.sh** (Recommended)
   - ATAC-bin mode for low-coverage data
   - Use this if pileup has sparse SNP calls

## Quick Start (ATAC-bin - Recommended)

```bash
qsub prepare_numbat_inputs_lowseq_489.qsub.sh
qsub run_numbat_atac_bin_lowseq_489.qsub.sh
```

## Key Files

### Inputs
- BAM: `Data/01_inputs/bam/lowseq_489.bam.lnk`
- Barcodes: `Data/01_inputs/barcodes/tissue_barcodes/lowseq_489/lowseq_489.barcodes.tsv`

### Outputs
- ATAC-bin: `Data/04_analysis/cnv/numbat/results/lowseq_489_atac_bin/numbat_seurat_obj.RDS`
- Standard: `Data/04_analysis/cnv/numbat/results/lowseq_489/numbat_seurat_obj.RDS`

## See Also
- `../lowseq_488B/README.md` - Detailed workflow guide
- `../ORGANIZATION.md` - Full documentation
