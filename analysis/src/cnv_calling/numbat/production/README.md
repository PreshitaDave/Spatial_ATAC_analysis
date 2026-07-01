# NUMBAT Production Scripts by Tissue

This folder contains production-ready NUMBAT analysis scripts, organized by tissue type and sequencing depth.

## Quick Reference

| Tissue | Sequencing | Cells | Scripts |
|--------|-----------|-------|---------|
| **deepseq_488B** | Deep (~190M reads) | ~12K | prepare → analyze |
| **deepseq_489** | Deep (~50M reads) | ~5K | prepare → analyze |
| **lowseq_488B** | Low (~8M reads) | ~12K | prepare → atac_bin OR analyze |
| **lowseq_489** | Low (~8M reads) | ~5K | prepare → atac_bin OR analyze |

## Tissue Folders

Each tissue folder contains the scripts needed for that analysis:

- `prepare_numbat_inputs_{tissue}.qsub.sh` - Data prep (allele counts)
- `run_numbat_analysis_{tissue}.qsub.sh` - Main analysis (all tissues)
- `run_numbat_atac_bin_{tissue}.qsub.sh` - ATAC-bin workflow (lowseq only)

## Workflow

### Standard Workflow (All Tissues)

```bash
cd {tissue}

# Step 1: Prepare
qsub prepare_numbat_inputs_{tissue}.qsub.sh

# Step 2: Analyze
qsub run_numbat_analysis_{tissue}.qsub.sh

# Results: Data/04_analysis/cnv/numbat/results/{tissue}/
```

### ATAC-bin Workflow (Lowseq Only)

Use ATAC-seq binning instead of SNP-based for better performance on low-coverage data:

```bash
cd lowseq_{488B|489}

# Step 1: Prepare
qsub prepare_numbat_inputs_lowseq_{488B|489}.qsub.sh

# Step 2: ATAC-bin Analysis
qsub run_numbat_atac_bin_lowseq_{488B|489}.qsub.sh

# Results: Data/04_analysis/cnv/numbat/results/lowseq_{488B|489}_atac_bin/
```

## See Also

- `../lib/` - Helper scripts and utilities
- `../ORGANIZATION.md` - Complete folder structure and documentation
- `../README.md` - Original NUMBAT workflow guide
