# NUMBAT & Alleloscope Complete Analysis Guide

## Overview

This document summarizes the complete workflow for:
1. **NUMBAT**: Copy Number Variation (CNV) detection from scATAC-seq data
2. **Alleloscope**: Single-cell haplotype inference from variant data

Both pipelines have been tested and validated with proper error handling for barcode format compatibility issues.

---

## Critical Issue: 10X Barcode Format Mismatch

### The Problem
- **Fragment files** from 10X Cell Ranger pipeline include barcode suffix: `TGGCTTCAAGCCATGC-1`
- **Cell barcode files** contain raw barcodes without suffix: `TGGCTTCAAGCCATGC`
- **Result without fix**: barcode matching returns 0 fragments → empty output files (182 bytes)
- **Example failure**: ATAC bin matrix should be ~100M, but only 182 bytes

### The Solution
Use the patched script: `analysis/src/cnv_calling/numbat/get_binned_atac_fixed.R`

**Key fix**: Strip `-1` or `-2` suffix before barcode matching
```R
barcodes_clean <- sub("-[12]$", "", barcodes_from_fragments)
# Now matches work: TGGCTTCAAGCCATGC-1 → TGGCTTCAAGCCATGC
```

### Diagnostic Checks
Always run these before starting analysis:
```bash
# Check fragment barcodes (should have -1 suffix)
zcat Data/01_inputs/fragments/{TISSUE}/*.bed.gz | head -5 | cut -f4 | sort -u
# Output: XXXXXX-1

# Check cell barcode file (should NOT have suffix)
head Data/01_inputs/barcodes/tissue_barcodes/{TISSUE}/*.tsv | cut -f1 | sort -u
# Output: XXXXXX

# Check ATAC bin matrix output size (should be >100M, NOT <1K)
ls -lh Data/04_analysis/cnv/numbat/inputs/*_atac_bin.rds
# Should show: -rw-r-- ... 98M ... (NOT 182 bytes)
```

---

## NUMBAT Workflow

### Prerequisite: Understanding NUMBAT
NUMBAT combines:
1. **Variant calling** (pileup stage) - copies variants from BAM file
2. **Long-range phasing** - uses genetic maps + reference panels
3. **ATAC binning** - aggregates ATAC reads into 220kb windows
4. **Reference comparison** - compares against normal cell baseline
5. **CNV inference** - predicts copy number per cell

### Full Workflow Steps

#### Step 1: Generate Lambda References (First Time Only)
```bash
# On login node (scc1):
cd /projectnb/paxlab/presh/projects/spatial_atac

# Allocate 8 cores for 2.5 hours
qrsh -l h_rt=02:30:00 -pe omp 8 -P paxlab -l mem_per_core=8G

# Inside compute node, create/attach tmux session
tmux new-session -s numbat_refs
# OR: tmux attach -t spatial_atac_work

# Verify compute node (should NOT be scc1)
hostname
# Output: scc-tb3, scc-tc3, etc.

# Run setup script
cd /projectnb/paxlab/presh/projects/spatial_atac
bash analysis/setup_numbat_env.sh

# Generate all tissue references (creates tissue-specific lambdas files)
bash analysis/src/cnv_calling/numbat/generate_all_tissue_references.test.sh

# Verify output
ls -lh Data/04_analysis/cnv/numbat/reference/lambdas_*_ATAC_bincnt.rds
# Should show 4 files: lowseq_488B, lowseq_489, deepseq_488B, deepseq_489

# Detach from tmux (keeps session alive)
# Press: Ctrl+B then D
```

#### Step 2: Generate Pileup & Phasing (Variant Calling)
```bash
# From login node, submit for each tissue:
qsub analysis/src/cnv_calling/numbat/run_numbat_pileup_phase_lowseq_488B.qsub.sh
qsub analysis/src/cnv_calling/numbat/run_numbat_pileup_phase_lowseq_489.qsub.sh
qsub analysis/src/cnv_calling/numbat/run_numbat_pileup_phase_deepseq_488B.qsub.sh
qsub analysis/src/cnv_calling/numbat/run_numbat_pileup_phase_deepseq_489.qsub.sh

# Monitor progress
qstat -u preshita
tail -f analysis/qsub_logs/numbat/run_numbat_pileup_phase_lowseq_489_*.log

# Expected output:
# Data/04_analysis/cnv/numbat/inputs/lowseq_489_comb_allele_counts.tsv.gz (~50M)
# Data/04_analysis/cnv/numbat/inputs/lowseq_489_atac_bin.rds (~100M)
```

#### Step 3: Run CNV Analysis
```bash
# Only after pileup jobs complete successfully, submit analysis:
qsub analysis/src/cnv_calling/numbat/run_numbat_analysis_lowseq_488B.qsub.sh
qsub analysis/src/cnv_calling/numbat/run_numbat_analysis_lowseq_489.qsub.sh
qsub analysis/src/cnv_calling/numbat/run_numbat_analysis_deepseq_488B.qsub.sh
qsub analysis/src/cnv_calling/numbat/run_numbat_analysis_deepseq_489.qsub.sh

# Monitor results
qstat -u preshita
tail -f analysis/qsub_logs/numbat/run_numbat_analysis_lowseq_489_*.log

# Check outputs once complete
ls -lh Data/04_analysis/cnv/numbat/results/lowseq_489/
# Should contain: seurat_obj_updated.rds, plot_*.pdf, etc.
```

### Input Requirements

| File | Location | Format | Notes |
|------|----------|--------|-------|
| BAM | `Data/01_inputs/bam/{tissue}.bam` | BAM | Indexed, all barcodes before filtering |
| Barcodes | `Data/01_inputs/barcodes/tissue_barcodes/{tissue}/{tissue}.barcodes.tsv` | TSV | Raw format, no -1 suffix |
| Fragments | `Data/01_inputs/fragments/{tissue}/{tissue}.fragments.sort.filtered.bed.gz` | BED.GZ | Cell Ranger format, has -1 suffix |
| Reference Genome | `Data/02_references/` | FASTA | hg38 required |

### Output Locations

| Output | Location | Size | Description |
|--------|----------|------|-------------|
| Lambda Reference | `Data/04_analysis/cnv/numbat/reference/lambdas_{tissue}_ATAC_bincnt.rds` | ~180K | Aggregated ATAC baseline |
| Allele Counts | `Data/04_analysis/cnv/numbat/inputs/alleles/{tissue}_comb_allele_counts.tsv.gz` | ~50M | Variant counts per cell |
| ATAC Matrix | `Data/04_analysis/cnv/numbat/inputs/{tissue}_atac_bin.rds` | ~100M | Cell × 220kb bin matrix |
| CNV Results | `Data/04_analysis/cnv/numbat/results/{tissue}/seurat_obj_updated.rds` | ~200M | Final CNV calls + metadata |
| Plots | `Data/04_analysis/cnv/numbat/results/{tissue}/plot_*.pdf` | Various | Publication-ready visualizations |

### Key Resources

- **Lambda Reference File**: `lambdas_{tissue}_ATAC_bincnt.rds` (178K per tissue)
  - Generated from that tissue's normal ATAC data
  - Used as baseline for CNV detection
  - Must exist BEFORE running CNV analysis

- **Genomic Bins**: `var220kb.rds` (80K, shared across all tissues)
  - 220 kilobase windows covering entire genome
  - Location: `Data/04_analysis/cnv/numbat/reference/`

- **Phasing Panel**: `phased_panel_bcf_links/` (shared directory)
  - Eagle phasing reference panel
  - Pre-downloaded and available

---

## Alleloscope Workflow

### Overview
Alleloscope infers single-cell haplotypes from variant data. Useful for:
- Validating CNV calls from NUMBAT
- Phase-aware analysis of tumor clones
- Understanding allele-specific copy number changes

### Workflow Steps

#### Step 1: Prepare Variant Inputs
```bash
# From login node, submit for each tissue:
qsub analysis/qsub/alleloscope/lowseq/prepare_alleloscope_lowseq_488B.qsub.sh
qsub analysis/qsub/alleloscope/lowseq/prepare_alleloscope_lowseq_489.qsub.sh
qsub analysis/qsub/alleloscope/deepseq/prepare_alleloscope_deepseq_488B.qsub.sh
qsub analysis/qsub/alleloscope/deepseq/prepare_alleloscope_deepseq_489.qsub.sh

# Monitor
qstat -u preshita
tail -f analysis/qsub_logs/alleloscope/*.log

# Expected outputs:
# Data/04_analysis/cnv/alleloscope/inputs/{tissue}_variant_matrix.tsv
```

#### Step 2: Run Haplotype Inference
```bash
# After preparation completes:
qsub analysis/qsub/alleloscope/lowseq/run_alleloscope_lowseq_488B.qsub.sh
qsub analysis/qsub/alleloscope/lowseq/run_alleloscope_lowseq_489.qsub.sh
qsub analysis/qsub/alleloscope/deepseq/run_alleloscope_deepseq_488B.qsub.sh
qsub analysis/qsub/alleloscope/deepseq/run_alleloscope_deepseq_489.qsub.sh

# Check results
ls -lh Data/04_analysis/cnv/alleloscope/results/lowseq_489/
# Should contain: haplotype_assignments.csv, phase_probabilities.csv, etc.
```

### Key Points for Alleloscope

- **Use tumor tissue ONLY** - reference samples may skew results
- **Consistent barcode format** - ensure all inputs use same format (no mixing of -1 suffixed vs raw)
- **Fragment data aggregation** - aggregate variants across spatial barcodes before input
- **Validation** - cross-reference with NUMBAT results for clonal CNV patterns

---

## Important Workflow Rules

### Always On Compute Node (NEVER Login Node)
```bash
# WRONG (on login node scc1):
bash analysis/src/cnv_calling/numbat/generate_all_tissue_references.test.sh

# CORRECT (on compute node via qrsh):
qrsh -l h_rt=02:30:00 -pe omp 8 -P paxlab -l mem_per_core=8G
# Then inside compute node:
bash analysis/src/cnv_calling/numbat/generate_all_tissue_references.test.sh
```

### Always Use Absolute Paths
```bash
# WRONG (relative):
/path/to/analysis/src/script.R --frag ./Data/fragments.bed.gz

# CORRECT (absolute):
/path/to/analysis/src/script.R --frag /projectnb/paxlab/presh/projects/spatial_atac/Data/fragments.bed.gz
```

### Always Verify Inputs Before Running
```bash
# Before any submission, verify these exist:
ls -lh Data/01_inputs/fragments/{TISSUE}/*.bed.gz  # Fragment file
ls -lh Data/01_inputs/barcodes/tissue_barcodes/{TISSUE}/*.tsv  # Barcode file
ls -lh Data/01_inputs/bam/{TISSUE}.bam  # BAM file
ls -lh Data/04_analysis/cnv/numbat/reference/lambdas_*ATAC_bincnt.rds  # References
```

### Use tmux for Long-Running Sessions
```bash
# Create new session
tmux new-session -s numbat_refs

# OR attach to existing
tmux list-sessions
tmux attach -t spatial_atac_work

# Detach (keeps running): Ctrl+B then D
# Kill session: tmux kill-session -t session_name
```

### Always Include Module Initialization
```bash
#!/bin/bash
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

module load R
which Rscript  # Verify it worked
```

---

## Troubleshooting

### Problem: "0 fragments matched"
**Cause**: Barcode format mismatch (fragment -1 suffix vs raw barcode)
```bash
# Check formats
zcat Data/01_inputs/fragments/{TISSUE}/*.bed.gz | head -3 | cut -f4
# Should show: XXXXXX-1 (fragment files)

head Data/01_inputs/barcodes/tissue_barcodes/{TISSUE}/*.tsv | cut -f1
# Should show: XXXXXX (barcode files)

# If mismatch exists, ensure script uses get_binned_atac_fixed.R
ls -l analysis/src/cnv_calling/numbat/get_binned_atac_fixed.R
```

### Problem: ATAC matrix output is tiny (<10MB)
**Cause**: Usually 0 fragments matched or R error
```bash
# Check the generation log
tail -100 Data/04_analysis/cnv/numbat/reference/{TISSUE}_reference_generation.log

# Look for error messages or 0 fragment count
grep -i "error\|matched\|warning" Data/04_analysis/cnv/numbat/reference/*.log
```

### Problem: Job submission hangs
**Cause**: Usually waiting for node allocation
```bash
# Check queue
qstat -u preshita

# Check if already have nodes
qstat -l -u preshita | grep -v Error

# Kill stuck job if needed
qdel {JOBID}
```

### Problem: Module not found on compute node
**Cause**: Module initialization script not sourced
```bash
# Inside compute node, manually test:
. /etc/profile
module load R
which Rscript
Rscript --version

# All should work without errors
```

---

## Quick Reference Commands

```bash
# Verify compute node allocation
hostname  # Should NOT start with scc1

# Create/attach tmux session
tmux new-session -s numbat_refs
tmux attach -t numbat_refs

# Check all reference files exist
ls -lh Data/04_analysis/cnv/numbat/reference/lambdas_*_ATAC_bincnt.rds

# Generate references for all tissues
bash analysis/src/cnv_calling/numbat/generate_all_tissue_references.test.sh

# Submit NUMBAT analysis jobs
qsub analysis/src/cnv_calling/numbat/run_numbat_analysis_*.qsub.sh

# Monitor job progress
qstat -u preshita
tail -f analysis/qsub_logs/numbat/*.log

# Check final outputs
ls -lh Data/04_analysis/cnv/numbat/results/lowseq_489/seurat_obj_updated.rds
```

---

## File Organization Summary

```
Data/
├── 01_inputs/
│   ├── fragments/
│   │   ├── lowseq_488B/
│   │   │   └── lowseq_488B.fragments.sort.filtered.bed.gz
│   │   ├── lowseq_489/
│   │   └── ...
│   ├── barcodes/tissue_barcodes/
│   │   ├── lowseq_488B/
│   │   │   └── lowseq_488B.barcodes.tsv
│   │   └── ...
│   ├── bam/
│   │   ├── lowseq_488B.bam
│   │   └── ...
│   └── archive/
│       └── test_files_archive_YYYYMMDD/
│
├── 02_references/
│   └── [genome reference files]
│
└── 04_analysis/
    └── cnv/
        ├── numbat/
        │   ├── reference/  (Generated outputs)
        │   │   ├── lambdas_lowseq_488B_ATAC_bincnt.rds
        │   │   ├── lambdas_lowseq_489_ATAC_bincnt.rds
        │   │   ├── lambdas_deepseq_488B_ATAC_bincnt.rds
        │   │   ├── lambdas_deepseq_489_ATAC_bincnt.rds
        │   │   ├── var220kb.rds
        │   │   ├── par_numbatm.rds
        │   │   └── phased_panel_bcf_links/
        │   │
        │   ├── inputs/
        │   │   ├── lowseq_489_atac_bin.rds
        │   │   ├── lowseq_489_comb_allele_counts.tsv.gz
        │   │   └── alleles/
        │   │
        │   └── results/
        │       ├── lowseq_489/
        │       │   ├── seurat_obj_updated.rds
        │       │   └── plot_*.pdf
        │       └── ...
        │
        └── alleloscope/
            ├── inputs/
            │   └── {tissue}_variant_matrix.tsv
            │
            └── results/
                └── {tissue}/
                    ├── haplotype_assignments.csv
                    └── phase_probabilities.csv

analysis/
├── src/cnv_calling/numbat/
│   ├── get_binned_atac_fixed.R  (PATCHED - always use this)
│   ├── generate_all_tissue_references.test.sh  (Run in qrsh)
│   ├── run_numbat_pileup_phase_*.qsub.sh
│   └── run_numbat_analysis_*.qsub.sh
│
├── setup_numbat_env.sh  (Environment setup)
├── NUMBAT_REFERENCE_GENERATION_WORKFLOW.md  (Detailed guide)
└── COMPLETE_ANALYSIS_GUIDE.md  (This file)
```

---

## References

- **NUMBAT Paper**: Jiang et al. 2024 - "Single-cell CNV detection using ATAC-seq"
- **NUMBAT GitHub**: https://github.com/kharchenkolab/numbat
- **Alleloscope GitHub**: https://github.com/camplab/alleloscope
- **10X Cell Ranger**: https://support.10xgenomics.com/single-cell-atac/software

---

## Questions or Issues?

See the following documents:
1. `.github/copilot-instructions.md` - Comprehensive project guidelines
2. `analysis/NUMBAT_REFERENCE_GENERATION_WORKFLOW.md` - Detailed reference generation steps
3. `analysis/setup_numbat_env.sh` - Environment verification script
4. `Data/ORGANIZATION_SUMMARY.md` - File organization reference

Last updated: May 18, 2026
