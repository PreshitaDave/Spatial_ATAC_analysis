# NUMBAT ATAC-bin Workflow for Spatial ATAC-seq

## Overview

This folder contains the complete pipeline for running NUMBAT (clonal copy-number variation analysis) in **ATAC-bin mode** on spatial ATAC-seq data. The workflow follows the [NUMBAT multiome guide](https://github.com/compbio-UofT/numbat/wiki/Running-NUMBAT-on-Multiome-data) and processes **lowseq** and **deepseq** datasets separately for tissues 488B and 489.

### Why ATAC-bin Mode?

- Uses chromatin accessibility patterns (ATAC-seq) instead of RNA expression
- Requires SNP allele counts + accessibility bins for robust CNV inference
- CRITICAL: Must include full pileup/phasing pipeline (not bin-only analysis)

---

## Scripts (In Order of Execution)

### 1. `test_numbat_inputs.test.sh` - **VALIDATE INPUTS (RUN FIRST)**

Comprehensive input validation before full analysis.

**Usage:**
```bash
bash test_numbat_inputs.test.sh <dataset> <tissue>
```

**Examples:**
```bash
bash test_numbat_inputs.test.sh lowseq 488B
bash test_numbat_inputs.test.sh deepseq 489
```

**What it tests:**
- ✓ BAM files exist and are readable (`samtools view -H`)
- ✓ Barcode files in correct format (one per line)
- ✓ SNP VCF accessible with variant count validation
- ✓ Genetic map present (gzip format)
- ✓ Genomic bins (var220kb.rds) RDS file readable
- ✓ Reference directory structure and permissions

**Runtime:** 2-5 minutes  
**Required before:** prepare_numbat_inputs.sh  
**On Failure:** Debug issues and re-run before proceeding

---

### 2. `prepare_numbat_inputs.sh` - **PREPARE INPUTS (MAIN PIPELINE)**

Generates all required files:
1. **Pileup & Phase** → SNP allele counts per cell
2. **ATAC Bin Count** → 220kb bin × cell matrix
3. **Generate Reference** → Aggregated reference from normal cells

**Usage:**
```bash
bash prepare_numbat_inputs.sh <dataset> <tissue>
```

**Examples:**
```bash
# Prepare lowseq tissue 488B
bash prepare_numbat_inputs.sh lowseq 488B

# Prepare deepseq tissue 489
bash prepare_numbat_inputs.sh deepseq 489
```

**Outputs Generated:**
- `Data/04_analysis/cnv/numbat/inputs/alleles/${DATASET}_${TISSUE}_comb_allele_counts.tsv.gz` (SNP alleles)
- `Data/04_analysis/cnv/numbat/inputs/${DATASET}_${TISSUE}_atac_bin.rds` (bin matrix)
- `Data/02_references/numbat/lambdas_ATAC_bincnt.rds` (reference, shared)

**Runtime per tissue:** 2-4 hours  
**Parallelization:** Run tissues 488B and 489 simultaneously (separate SGE jobs)  
**NCORES:** 8 (configured internally)

**Prerequisites:**
- ✅ Merged BAM files: `Data/04_analysis/cnv/numbat/inputs/bam_merged/`
- ✅ Barcode files: `Data/01_inputs/barcodes/tissue_barcodes/`
- ✅ SNP VCF: `Data/02_references/numbat/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf`
- ✅ Genetic map: `Data/02_references/numbat/genetic_map_hg38_withX.txt.gz`
- ✅ Bins: `Data/02_references/numbat/var220kb.rds`

---

### 3. `run_numbat_analysis_atac.sh` - **RUN NUMBAT ANALYSIS**

Performs CNV inference on prepared inputs.

**Usage:**
```bash
bash run_numbat_analysis_atac.sh <dataset> <tissue>
```

**Outputs:**
- `Data/04_analysis/cnv/numbat/results/${DATASET}_${TISSUE}_atac/cnv_calls.rds`
- `Data/04_analysis/cnv/numbat/results/${DATASET}_${TISSUE}_atac/phylogeny.png`
- `Data/04_analysis/cnv/numbat/results/${DATASET}_${TISSUE}_atac/numbat_${DATASET}_${TISSUE}.log`

**Runtime per tissue:** 4-8 hours  
**Prerequisites:** All outputs from prepare_numbat_inputs.sh

---

## Complete Workflow

### For Lowseq Dataset

```bash
cd /projectnb/paxlab/presh/projects/spatial_atac

# STEP 1: Validate inputs for 488B
bash analysis/src/cnv_calling/numbat/test_numbat_inputs.test.sh lowseq 488B
# Expected: All tests PASS

# STEP 2: Validate inputs for 489
bash analysis/src/cnv_calling/numbat/test_numbat_inputs.test.sh lowseq 489
# Expected: All tests PASS

# STEP 3: Prepare inputs (run in parallel)
# Terminal 1:
bash analysis/src/cnv_calling/numbat/prepare_numbat_inputs.sh lowseq 488B

# Terminal 2 (simultaneously):
bash analysis/src/cnv_calling/numbat/prepare_numbat_inputs.sh lowseq 489

# STEP 4: Once preparation complete, run NUMBAT analysis (sequential)
bash analysis/src/cnv_calling/numbat/run_numbat_analysis_atac.sh lowseq 488B
bash analysis/src/cnv_calling/numbat/run_numbat_analysis_atac.sh lowseq 489
```

### For Deepseq Dataset

Identical workflow, replace `lowseq` with `deepseq`:

```bash
bash analysis/src/cnv_calling/numbat/test_numbat_inputs.test.sh deepseq 488B
bash analysis/src/cnv_calling/numbat/prepare_numbat_inputs.sh deepseq 488B
bash analysis/src/cnv_calling/numbat/run_numbat_analysis_atac.sh deepseq 488B
```

---

## File Organization

**Input Files:**
```
Data/
├── 01_inputs/barcodes/tissue_barcodes/
│   ├── lowseq_488B/{sample}.barcodes.tsv
│   ├── lowseq_489/{sample}.barcodes.tsv
│   └── [deepseq_* tissue directories]
├── 02_references/numbat/
│   ├── genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf (3.76M SNPs)
│   ├── genetic_map_hg38_withX.txt.gz (54M, symlink)
│   └── var220kb.rds (220kb genomic bins)
└── 04_analysis/cnv/numbat/inputs/
    ├── bam_merged/
    │   ├── lowseq_488B_merged_for_numbat.bam (indexed)
    │   └── lowseq_489_merged_for_numbat.bam (indexed)
    └── alleles/ (GENERATED)
        ├── lowseq_488B_comb_allele_counts.tsv.gz
        └── lowseq_489_comb_allele_counts.tsv.gz
```

**Output Files:**
```
Data/04_analysis/cnv/numbat/
├── inputs/ (intermediate)
│   ├── lowseq_488B_atac_bin.rds (307M)
│   └── lowseq_489_atac_bin.rds (96M)
├── results/ (final)
│   ├── lowseq_488B_atac/
│   │   ├── cnv_calls.rds (clone CNVs)
│   │   ├── phylogeny.png (clone tree)
│   │   └── numbat_lowseq_488B.log
│   └── lowseq_489_atac/
│       ├── cnv_calls.rds
│       ├── phylogeny.png
│       └── numbat_lowseq_489.log
└── references/ (shared)
    └── lambdas_ATAC_bincnt.rds (ATAC reference)
```

---

## Key Parameters

| Parameter | Value | Purpose |
|-----------|-------|---------|
| **Mode** | ATAC-bin | Chromatin-only (no RNA) |
| **Resolution** | 220kb | Genomic bin size |
| **Min Cells/Clone** | 5 | Minimum clone size threshold |
| **NCORES** | 8 | Parallelization threads |
| **Genetic Map** | hg38_withX | For phasing |
| **Ref VCF** | 1KG Phase 3 (3.76M SNPs) | Variant annotation |

---

## Troubleshooting

### Test Script Fails

**Problem:** Test script reports missing files

**Solution:**
```bash
# Check BAM files exist
ls -lh Data/04_analysis/cnv/numbat/inputs/bam_merged/lowseq_*.bam

# Check barcodes exist
ls -lh Data/01_inputs/barcodes/tissue_barcodes/lowseq_*/*.tsv

# Check reference files
ls -lh Data/02_references/numbat/
```

### Pileup/Phasing Fails

**Problem:** Step 1 of prepare_numbat_inputs.sh fails with barcode/BAM error

**Solution:**
- Verify barcode file format: `head Data/01_inputs/barcodes/tissue_barcodes/lowseq_488B/*`
- Check BAM index: `samtools index Data/04_analysis/cnv/numbat/inputs/bam_merged/lowseq_488B_merged_for_numbat.bam`
- Verify BAM readability: `samtools view -c Data/04_analysis/cnv/numbat/inputs/bam_merged/lowseq_488B_merged_for_numbat.bam`

### Low Variant Counts

**Problem:** Allele counts file has few variants

**This is expected** if tissue coverage is sparse. NUMBAT handles sparse variant data.

### Reference Generation Fails

**Problem:** Reference generation in step 3 fails

**Solution:**
- Verify bin matrix RDS readability: Run test script
- Check NUMBAT bin function: `Rscript -e "require(numbat); ?get_binned_atac"`

---

## References

- **NUMBAT GitHub:** https://github.com/compbio-UofT/numbat
- **NUMBAT Multiome Guide:** https://github.com/compbio-UofT/numbat/wiki/Running-NUMBAT-on-Multiome-data
- **Project Root:** `/projectnb/paxlab/presh/projects/spatial_atac`

---

## Change History

| Date | Change | Status |
|------|--------|--------|
| 2026-05-17 | Initial lowseq analysis (INCORRECT) | ❌ Archived - missing pileup/phasing |
| 2026-05-17 | Corrected workflow per NUMBAT guide | ✅ Ready for production |

**Notes:**
- Previous incorrect analysis archived at: `Data/01_inputs/archive/incorrect_numbat_results_20260517/`
- Reason: Missing pileup/phasing step, incomplete variant calling
- All paths use absolute paths (required for SGE jobs)
- Validation is MANDATORY before full production runs
- Parallelization implemented at tissue level and core level
