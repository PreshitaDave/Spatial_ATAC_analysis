# NUMBERED GUIDE: NUMBAT ATAC-bin Analysis Pipeline

## Project Context
- **Project:** Spatial ATAC-seq Analysis
- **Datasets:** Deepseq vs Lowseq for patients 448B and 489
- **Objective:** Clone inference via CNV calling in ATAC-bin mode
- **Workflow:** NUMBAT multiome ATAC-bin pipeline

---

## Phase 1: Preparation (Weeks 1-2)

### 1. Input Data Organization

**Files Required:**
1. Merged BAM files (tissue-specific, all barcodes)
   - Location: `Data/04_analysis/cnv/numbat/inputs/bam_merged/`
   - Format: `{DATASET}_{TISSUE}_merged_for_numbat.bam` (indexed)
   - Size: 8-15 GB per tissue

2. Barcode files (one per line, tissue-specific)
   - Location: `Data/01_inputs/barcodes/tissue_barcodes/{DATASET}_{TISSUE}/`
   - Format: `.barcodes.tsv`
   - Content: One barcode per line (no header)

3. Reference VCF (SNPs for pileup)
   - Location: `Data/02_references/numbat/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf`
   - Size: 764 MB uncompressed
   - Content: 3.76M high-frequency SNPs from 1000 Genomes Phase 3

4. Genetic map (for phasing)
   - Location: `Data/02_references/numbat/genetic_map_hg38_withX.txt.gz`
   - Size: 54 MB (symlinked from Eagle)
   - Content: Recombination rates per megabase

5. Genomic bins (220kb resolution)
   - Location: `Data/02_references/numbat/var220kb.rds`
   - Format: GRanges object (RDS)
   - Content: 220kb non-overlapping genomic bins

**Action:**
```bash
# Verify all input files exist
cd /projectnb/paxlab/presh/projects/spatial_atac
ls -lh Data/04_analysis/cnv/numbat/inputs/bam_merged/lowseq_*.bam
ls -lh Data/01_inputs/barcodes/tissue_barcodes/lowseq_*/
ls -lh Data/02_references/numbat/
```

---

### 2. Input Validation

**Script:** `test_numbat_inputs.test.sh`  
**Purpose:** Verify all inputs before full processing  
**Runtime:** 2-5 minutes

**Action:**
```bash
cd /projectnb/paxlab/presh/projects/spatial_atac
bash analysis/src/cnv_calling/numbat/test_numbat_inputs.test.sh lowseq 488B
bash analysis/src/cnv_calling/numbat/test_numbat_inputs.test.sh lowseq 489
```

**Expected Output:**
```
✓ BAM file readable
✓ Barcode file format correct
✓ VCF accessible (XXXX variants)
✓ Genetic map present
✓ Bins (var220kb.rds) readable
✓ Reference directory ready
```

**On Failure:** Stop and debug (do NOT proceed to step 3)

---

## Phase 2: Input Preparation (Weeks 2-3)

### 3. Generate Pileup & Allele Counts

**Script:** `prepare_numbat_inputs.sh` (Step 1)  
**Purpose:** Extract SNP alleles at each cell  
**Output:** `${DATASET}_${TISSUE}_comb_allele_counts.tsv.gz` (~50-200 MB)  
**Runtime:** 2-4 hours per tissue

**Key Function:**
```bash
# Internally calls:
Rscript ${NUMBAT_REPO}/inst/bin/pileup_and_phase.R \
  --label lowseq_488B \
  --samples lowseq_488B \
  --bams Data/04_analysis/cnv/numbat/inputs/bam_merged/lowseq_488B_merged_for_numbat.bam \
  --barcodes Data/01_inputs/barcodes/tissue_barcodes/lowseq_488B/lowseq_488B.barcodes.tsv \
  --snpvcf Data/02_references/numbat/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
  --gmap Data/02_references/numbat/genetic_map_hg38_withX.txt.gz \
  --ncores 8 \
  --cellTAG CB \
  --UMItag None
```

**Action (Parallel Submission):**
```bash
# Terminal 1:
bash analysis/src/cnv_calling/numbat/prepare_numbat_inputs.sh lowseq 488B

# Terminal 2 (simultaneous):
bash analysis/src/cnv_calling/numbat/prepare_numbat_inputs.sh lowseq 489
```

**OR as SGE jobs:**
```bash
cd /projectnb/paxlab/presh/projects/spatial_atac
qsub analysis/src/cnv_calling/numbat/prepare_numbat_inputs_lowseq_488B.qsub.sh &
qsub analysis/src/cnv_calling/numbat/prepare_numbat_inputs_lowseq_489.qsub.sh &
wait
```

**Monitor Progress:**
```bash
# Check if allele counts files exist
watch -n 5 'ls -lh Data/04_analysis/cnv/numbat/inputs/alleles/'
```

**Expected Output:**
- `Data/04_analysis/cnv/numbat/inputs/alleles/lowseq_488B_comb_allele_counts.tsv.gz` (100-200 MB)
- `Data/04_analysis/cnv/numbat/inputs/alleles/lowseq_489_comb_allele_counts.tsv.gz` (50-100 MB)

---

### 4. Generate ATAC Bin Matrices

**Script:** `prepare_numbat_inputs.sh` (Step 2)  
**Purpose:** Create 220kb bin × cell count matrix  
**Output:** `${DATASET}_${TISSUE}_atac_bin.rds`  
**Runtime:** 1-2 hours per tissue

**Prerequisites:** Step 3 complete (allele counts exist)

**Action:** Already included in prepare_numbat_inputs.sh - no separate step needed

**Expected Output:**
- `Data/04_analysis/cnv/numbat/inputs/lowseq_488B_atac_bin.rds` (307 MB)
- `Data/04_analysis/cnv/numbat/inputs/lowseq_489_atac_bin.rds` (96 MB)

---

### 5. Generate ATAC Reference

**Script:** `prepare_numbat_inputs.sh` (Step 3)  
**Purpose:** Create aggregated reference from normal/diploid cells  
**Output:** `lambdas_ATAC_bincnt.rds` (shared across all tissues)  
**Runtime:** 30-60 minutes

**Prerequisites:** Step 4 complete (bin matrices exist)

**Expected Output:**
- `Data/02_references/numbat/lambdas_ATAC_bincnt.rds` (50-100 MB)

**Note:** Generated once, used for all tissues

---

## Phase 3: Analysis (Weeks 3-4)

### 6. Run NUMBAT ATAC-bin Analysis

**Script:** `run_numbat_analysis_atac.sh`  
**Purpose:** Call CNVs and infer clonal structure  
**Output:** CNV calls, phylogeny, clone assignments  
**Runtime:** 4-8 hours per tissue

**Prerequisites:** Steps 3-5 complete (all inputs ready)

**Action:**
```bash
# After allele counts + bin matrices + reference are ready:
bash analysis/src/cnv_calling/numbat/run_numbat_analysis_atac.sh lowseq 488B
bash analysis/src/cnv_calling/numbat/run_numbat_analysis_atac.sh lowseq 489
```

**Expected Outputs:**
```
Data/04_analysis/cnv/numbat/results/lowseq_488B_atac/
├── cnv_calls.rds (primary output)
├── phylogeny.png (clone tree visualization)
├── clones.rds (clone assignments)
└── numbat_lowseq_488B.log (processing log)
```

**Key Results:**
- Clone phylogeny tree (PNG image)
- Per-cell clone assignments
- CNV segments and copy numbers
- Clonal fractions

---

### 7. Validate NUMBAT Results

**Checks:**
```bash
# Verify output files exist
ls -lh Data/04_analysis/cnv/numbat/results/lowseq_488B_atac/

# Check file sizes (should be >0)
du -h Data/04_analysis/cnv/numbat/results/lowseq_488B_atac/*

# View phylogeny image
eog Data/04_analysis/cnv/numbat/results/lowseq_488B_atac/phylogeny.png

# Inspect log for errors
tail -50 Data/04_analysis/cnv/numbat/results/lowseq_488B_atac/numbat_lowseq_488B.log
```

---

## Phase 4: Repeat for Other Datasets

### 8. Run Deepseq Equivalent

Identical workflow, replace `lowseq` with `deepseq`:

```bash
# Validate inputs
bash analysis/src/cnv_calling/numbat/test_numbat_inputs.test.sh deepseq 488B
bash analysis/src/cnv_calling/numbat/test_numbat_inputs.test.sh deepseq 489

# Prepare inputs (parallel)
bash analysis/src/cnv_calling/numbat/prepare_numbat_inputs.sh deepseq 488B &
bash analysis/src/cnv_calling/numbat/prepare_numbat_inputs.sh deepseq 489 &
wait

# Run analysis
bash analysis/src/cnv_calling/numbat/run_numbat_analysis_atac.sh deepseq 488B
bash analysis/src/cnv_calling/numbat/run_numbat_analysis_atac.sh deepseq 489
```

**Expected Outputs:**
```
Data/04_analysis/cnv/numbat/results/
├── lowseq_488B_atac/
├── lowseq_489_atac/
├── deepseq_488B_atac/
└── deepseq_489_atac/
```

---

## Phase 5: Advanced Analysis (Weeks 4+)

### 9. Cross-Dataset Comparison

**Goal:** Compare lowseq vs deepseq clonal structures for same patient

**Steps:**
1. Load both NUMBAT results (`cnv_calls.rds` files)
2. Compare CNV patterns (should be similar - same tumor)
3. Check clone frequencies (may differ due to sampling bias)
4. Validate with somatic_chr SNV comparison

---

### 10. Integrated Multi-Omic Analysis

**Downstream Analyses:**
- Alleloscope (CNV refinement with RNA integration)
- pyCistopic (topic modeling)
- Cellwalkr (multi-omic correlation)
- Spatial visualization (CNV heatmaps over tissue)

---

## Troubleshooting Guide

### Issue: "module: command not found"
**Solution:** Scripts include module initialization. Verify:
```bash
. /etc/profile.d/modules.sh
module load R
which Rscript
```

### Issue: "BAM file not found"
**Solution:** Verify BAM path and index:
```bash
ls -lh Data/04_analysis/cnv/numbat/inputs/bam_merged/lowseq_488B_merged_for_numbat.bam*
samtools index Data/04_analysis/cnv/numbat/inputs/bam_merged/lowseq_488B_merged_for_numbat.bam
```

### Issue: "Barcode mismatch"
**Solution:** Verify barcode format (one per line, no header):
```bash
head -5 Data/01_inputs/barcodes/tissue_barcodes/lowseq_488B/*.tsv
wc -l Data/01_inputs/barcodes/tissue_barcodes/lowseq_488B/*.tsv
```

### Issue: "Low variant count in allele counts"
**This is expected** - sparse coverage is normal. NUMBAT handles it.

### Issue: "Reference generation fails"
**Solution:** Verify bin matrix readability:
```bash
Rscript -e "readRDS('Data/04_analysis/cnv/numbat/inputs/lowseq_488B_atac_bin.rds'); cat('OK\n')"
```

---

## Key Files Reference

| File | Purpose | Size |
|------|---------|------|
| SNP VCF | Variant calling | 764 MB |
| Genetic Map | Phasing | 54 MB |
| Genomic Bins | Feature mapping | 5-10 MB |
| BAM (merged) | Read counts | 8-15 GB |
| Barcodes | Cell identifiers | <1 MB |
| Allele Counts | SNP alleles/cell | 50-200 MB |
| ATAC Bins | Bin×Cell matrix | 50-300 MB |
| Reference | Aggregated bins | 50-100 MB |
| CNV Calls | Output | 50-100 MB |

---

## Timeline Estimate

| Phase | Duration | Tasks |
|-------|----------|-------|
| 1. Prep | 1 day | Organize inputs, validate |
| 2. Generate | 7-10 days | Allele counts, bins, reference (parallel) |
| 3. Analyze | 3-5 days | NUMBAT for lowseq + deepseq |
| 4. Validate | 1-2 days | QC, comparisons |
| 5. Integration | 3-5 days | Multi-omic analysis, visualization |
| **Total** | **2-3 weeks** | Full pipeline completion |

---

## Best Practices

1. **Always validate inputs first** (Step 2) before proceeding
2. **Run tissues in parallel** where possible (488B and 489)
3. **Monitor large jobs** (Steps 3-6) with `qstat -u preshita`
4. **Archive old results** if re-running (prevents file conflicts)
5. **Keep logs** for troubleshooting and reproducibility
6. **Document any customizations** made to scripts

---

## References

- NUMBAT Multiome Guide: https://github.com/compbio-UofT/numbat/wiki/Running-NUMBAT-on-Multiome-data
- Project Root: `/projectnb/paxlab/presh/projects/spatial_atac`
- Scripts Location: `analysis/src/cnv_calling/numbat/`
- Data Location: `Data/04_analysis/cnv/numbat/`

---

**Document Version:** 1.0  
**Last Updated:** 2026-05-17  
**Status:** Production Ready  
**Maintained by:** Paxlab Spatial ATAC Team
