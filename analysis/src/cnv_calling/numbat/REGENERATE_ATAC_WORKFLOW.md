# NUMBAT ATAC Matrix Regeneration Workflow

**Last Updated**: May 19, 2026  
**Document Path**: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/cnv_calling/numbat/REGENERATE_ATAC_WORKFLOW.md`

---

## Overview

This document describes the ATAC matrix regeneration workflow for ensuring barcode consistency in NUMBAT CNV analysis. This workflow is required when different barcode files were used in different pipeline stages, causing cell count mismatches.

### Why This Matters

NUMBAT requires **exact cell-by-cell correspondence** between:
- **ATAC bin-by-cell matrix** (from fragment counts)
- **Allele counts** (from variant pileup)

If these inputs use different barcode sets:
- Result: 460+ cells in one file but not the other
- Effect: NUMBAT filters all mismatched cells to 0 coverage
- Outcome: Analysis fails with **no CNV calls**

### The Fix

Regenerate the ATAC matrix using the **same barcode file as the pileup stage**, ensuring 100% cell-by-cell correspondence.

---

## Quick Start (Single Command)

For any tissue, use the universal regeneration script:

```bash
# Submit regeneration job
qsub analysis/src/cnv_calling/numbat/regenerate_atac_full_barcodes.qsub.sh {TISSUE}

# Example for lowseq_489
qsub analysis/src/cnv_calling/numbat/regenerate_atac_full_barcodes.qsub.sh lowseq_489

# Example for deepseq_488B  
qsub analysis/src/cnv_calling/numbat/regenerate_atac_full_barcodes.qsub.sh deepseq_488B
```

The script automatically:
1. Finds all input files for the specified tissue
2. Uses the FULL barcode file (not filtered)
3. Generates ATAC matrix with ~4,600-4,700 cells (depending on tissue)
4. Validates output size
5. Creates backups

---

## Detailed Workflow

### Step 1: Understand the Barcode Hierarchy

Each tissue has multiple barcode files in `Data/01_inputs/barcodes/tissue_barcodes/{tissue}/`:

| File | Cells | Used For | Description |
|------|-------|----------|-------------|
| `{tissue}.barcodes.tsv` | ~4,671 | **Pileup stage** (allele generation) | ALL cells before filtering |
| `{tissue}.no_edge_effect.barcodes.tsv` | ~4,622 | Old ATAC (caused mismatch) | Cells after edge-effect removal |
| `{tissue}.edge_effect.barcodes.tsv` | ~49 | Reference/validation | Cells removed as edge effects |
| `{tissue}_nFrags_from_fragments.tsv.gz` | ~4,162 | Validation | Actual cells with fragments |

**KEY INSIGHT**: Must use the same barcode file for both pileup AND ATAC stages.

### Step 2: Verify Pileup Stage Used Correct Barcodes

Check which barcode file was used in the pileup job:

```bash
# Check pileup job submission (look at past qsub logs)
grep -r "barcodes.tsv" analysis/qsub_logs/*pileup*

# Or check the pileup results for cell count
wc -l Data/04_analysis/cnv/numbat/results/{tissue}/alleles.csv
# Should show ~4,671 rows (4,670 variants + 1 header)
```

### Step 3: Check Current ATAC Matrix Cell Count

```bash
# Load and inspect current ATAC matrix dimensions
cd /projectnb/paxlab/presh/projects/spatial_atac
Rscript << 'EOF'
library(Matrix)
atac <- readRDS("Data/04_analysis/cnv/numbat/inputs/{tissue}/atac_bin/{tissue}_atac_bin.rds")
cat("Current ATAC cells:", ncol(atac), "\n")
cat("Current ATAC bins:", nrow(atac), "\n")
EOF

# If cells < 4,600, regeneration is needed
```

### Step 4: Submit Regeneration Job

```bash
# For lowseq_489
qsub -N "regenerate_atac_lowseq_489" \
  analysis/src/cnv_calling/numbat/regenerate_atac_full_barcodes.qsub.sh lowseq_489

# For deepseq_489
qsub -N "regenerate_atac_deepseq_489" \
  analysis/src/cnv_calling/numbat/regenerate_atac_full_barcodes.qsub.sh deepseq_489

# For lowseq_488B
qsub -N "regenerate_atac_lowseq_488B" \
  analysis/src/cnv_calling/numbat/regenerate_atac_full_barcodes.qsub.sh lowseq_488B
```

### Step 5: Monitor Job Progress

```bash
# Check job status (while running)
qstat -j {JOB_ID}

# Monitor logs in real-time
tail -f analysis/qsub_logs/regenerate_atac_{JOB_ID}.out
tail -f analysis/qsub_logs/regenerate_atac_{JOB_ID}.err

# Check when done (should exit queue)
qstat -u preshita | grep regenerate
```

### Step 6: Validate Regeneration

Once the job completes, validate barcode consistency:

```bash
# Verify ATAC matrix was regenerated with correct cell count
Rscript analysis/src/cnv_calling/numbat/validate_numbat_inputs.R {tissue}

# Example output (if PASS):
# ────────────────────────────────────────────
# ✓ PASS: Barcode files are CONSISTENT
# ATAC cells: 4671
# Allele cells: 4671
# Overlap: 100%
# ────────────────────────────────────────────
```

### Step 7: Run NUMBAT Analysis

If validation **PASSES**, submit the analysis:

```bash
qsub analysis/src/cnv_calling/numbat/run_numbat_analysis_{tissue}_refhca.qsub.sh
```

---

## Script Details: `regenerate_atac_full_barcodes.qsub.sh`

### Input Parameters

The script takes ONE parameter:

```
Usage: qsub regenerate_atac_full_barcodes.qsub.sh {TISSUE}
       where {TISSUE} = lowseq_489, deepseq_489, lowseq_488B, deepseq_488B, etc.
```

The script automatically determines:
- Dataset type (lowseq vs deepseq) from tissue name
- All file paths based on tissue name
- Cell counts from barcode file

### Resource Allocation

```
Cores:     8 (OMP)
Memory:    8 GB per core (64 GB total)
Walltime:  4 hours (240 min)
Project:   paxlab
```

**Expected Duration**: 45-120 minutes (depends on fragment file size)

### Files It Creates/Modifies

| File | Purpose | Action |
|------|---------|--------|
| `Data/04_analysis/cnv/numbat/inputs/{tissue}/atac_bin/{tissue}_atac_bin.rds` | Active ATAC matrix | **REPLACED** with regenerated version |
| `Data/04_analysis/cnv/numbat/inputs/{tissue}/atac_bin/{tissue}_atac_bin_ORIGINAL_BACKUP.rds` | Backup | Created (for recovery if needed) |
| `analysis/qsub_logs/regenerate_atac_{JOB_ID}.out` | Standard output | Created (shows progress) |
| `analysis/qsub_logs/regenerate_atac_{JOB_ID}.err` | Error log | Created (should be small if success) |

### Output File Validation

After completion, verify:

```bash
# File should exist and be >50M
ls -lh Data/04_analysis/cnv/numbat/inputs/{tissue}/atac_bin/{tissue}_atac_bin.rds

# Load and inspect dimensions
Rscript << 'EOF'
library(Matrix)
atac <- readRDS("Data/04_analysis/cnv/numbat/inputs/{tissue}/atac_bin/{tissue}_atac_bin.rds")
cat("ATAC cells:", ncol(atac), "\n")  # Should be ~4,671
cat("ATAC bins:", nrow(atac), "\n")    # Should be ~12,145
EOF
```

---

## Reproducibility Checklist for Any New Tissue

Use this checklist when regenerating ATAC for a new tissue:

- [ ] **Tissue Identified**: Tissue name follows pattern: `{lowseq|deepseq}_{patient}`
- [ ] **Barcode Files Exist**: 
  - [ ] `Data/01_inputs/barcodes/tissue_barcodes/{tissue}/{tissue}.barcodes.tsv` exists
  - [ ] Check cell count: `wc -l Data/01_inputs/barcodes/tissue_barcodes/{tissue}/{tissue}.barcodes.tsv`
- [ ] **Fragment Files Exist**:
  - [ ] `Data/01_inputs/fragments/{tissue}/{tissue}.fragments.sort.filtered.bed.gz` exists  
  - [ ] Check file size: `ls -lh Data/01_inputs/fragments/{tissue}/`
- [ ] **Genomic Bins Exist**:
  - [ ] `Data/04_analysis/cnv/numbat/reference/var220kb.rds` exists (80K, shared across all tissues)
- [ ] **Output Directory Ready**:
  - [ ] `Data/04_analysis/cnv/numbat/inputs/{tissue}/atac_bin/` directory exists
  - [ ] Test writability: `touch Data/04_analysis/cnv/numbat/inputs/{tissue}/atac_bin/.test && rm .test`
- [ ] **Job Submission**:
  - [ ] Job queued without errors: `qsub ... {tissue}` returns Job ID
- [ ] **Job Completion**:
  - [ ] Check logs: `tail analysis/qsub_logs/regenerate_atac_{JOB_ID}.{out,err}`
  - [ ] Check file size: Output >50M, <2G (reasonable for ATAC matrix)
  - [ ] Verify cells: Use R to load and check `ncol(readRDS(...))` ≈ 4,600-4,700
- [ ] **Validation**:
  - [ ] Run: `Rscript validate_numbat_inputs.R {tissue}`
  - [ ] Expected: PASS with 100% barcode overlap
- [ ] **Analysis Ready**:
  - [ ] Submit: `qsub run_numbat_analysis_{tissue}_refhca.qsub.sh`

---

## Troubleshooting

### Symptom: Validation Still Shows Mismatch

**Cause**: Wrong barcode file used (script used filtered instead of full)

**Fix**:
1. Check which barcode file was actually used in pileup stage
2. Look at pileup job submission/logs for clues
3. Count allele file cells: `wc -l Data/04_analysis/cnv/numbat/results/{tissue}/alleles.csv`
4. Re-examine barcode files and use the matching one

### Symptom: Job Fails with "binGR not found"

**Cause**: Script invoked without `--binGR` parameter (happens if using old version)

**Fix**: Use the updated script: `regenerate_atac_full_barcodes.qsub.sh`

### Symptom: Output File <1M (Corrupted)

**Cause**: Barcode mismatch between fragments and barcode file

**Check**:
```bash
# Compare barcode formats
zcat Data/01_inputs/fragments/{tissue}/{tissue}.fragments.sort.filtered.bed.gz | \
  head -5 | cut -f4 | sort -u
# Should show: {16bp}-1, {16bp}-2, etc (with suffix)

head Data/01_inputs/barcodes/tissue_barcodes/{tissue}/{tissue}.barcodes.tsv | sort -u
# Should show: {16bp} (no suffix)

# The script handles this automatically with gsub("-[0-9]+$", "", ...)
```

### Symptom: Job Takes >4 Hours

**Cause**: Large fragment file or slow filesystem

**Solution**: Resubmit with more time:
```bash
qsub -l h_rt=06:00:00 analysis/src/cnv_calling/numbat/regenerate_atac_full_barcodes.qsub.sh {tissue}
```

---

## For Reference: Key Paths

```
# Inputs (by tissue)
Data/01_inputs/fragments/{tissue}/{tissue}.fragments.sort.filtered.bed.gz
Data/01_inputs/barcodes/tissue_barcodes/{tissue}/{tissue}.barcodes.tsv

# Shared reference
Data/04_analysis/cnv/numbat/reference/var220kb.rds

# Output (by tissue)  
Data/04_analysis/cnv/numbat/inputs/{tissue}/atac_bin/{tissue}_atac_bin.rds

# NUMBAT pipeline inputs (by tissue)
Data/04_analysis/cnv/numbat/inputs/{tissue}/alleles/
Data/04_analysis/cnv/numbat/inputs/{tissue}/barcodes/
Data/04_analysis/cnv/numbat/inputs/{tissue}/atac_bin/

# NUMBAT results (by tissue)
Data/04_analysis/cnv/numbat/results/{tissue}/
```

---

## Related Documents

- **NUMBAT Pipeline**: `.github/copilot-instructions.md` → "CNV Analysis Workflows" section
- **Barcode Consistency Issue**: `.github/copilot-instructions.md` → "CRITICAL ISSUE 2: Barcode File Consistency"
- **Validation Script**: `analysis/src/cnv_calling/numbat/validate_numbat_inputs.R`
- **NUMBAT Analysis**: `analysis/src/cnv_calling/numbat/run_numbat_analysis_*_refhca.qsub.sh`

---

## Change Log

| Date | Tissue | Job ID | Status | Notes |
|------|--------|--------|--------|-------|
| 2026-05-19 | lowseq_489 | 5723875 | SUBMITTED | ATAC regeneration with FULL barcodes (4,671 cells) |
| | | | | See: This workflow document |

