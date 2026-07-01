# NUMBAT Reference Generation Workflow

## Quick Start: Generate References for All Tissues

This guide shows how to generate lambda reference files for all tissues using qrsh and tmux.

### Prerequisites
- All fragment files in: `Data/01_inputs/fragments/{tissue}/`
- All barcode files in: `Data/01_inputs/barcodes/tissue_barcodes/{tissue}/`
- Fixed binning script: `analysis/src/cnv_calling/numbat/get_binned_atac_fixed.R`
- Reference generation script: `analysis/src/cnv_calling/numbat/generate_all_tissue_references.test.sh`

### Step-by-Step Workflow

#### 1. Check Existing tmux Sessions
```bash
# From login node (scc1), list all sessions
tmux list-sessions

# Output should show existing sessions like:
# spatial_atac_work: 1 windows
# spatial_atac_numbat: 1 windows
```

#### 2. Allocate Compute Node with qrsh
```bash
# Request 8 cores, 2.5 hours, 8GB RAM per core
qrsh -l h_rt=02:30:00 -pe omp 8 -P paxlab -l mem_per_core=8G

# You should see output like:
# scc-tb3  (or similar compute node)
```

#### 3. Create/Attach to Persistent tmux Session
```bash
# Create new session (if needed)
tmux new-session -s numbat_refs

# OR attach to existing session
tmux attach -t spatial_atac_work

# You are now in tmux - all commands persist even if connection drops
```

#### 4. Verify Compute Node Allocation (Inside tmux)
```bash
hostname
# Should output: scc-tb3, scc-tc3, or similar (NOT scc1)

pwd
# Should be: /projectnb/paxlab/presh/projects/spatial_atac
```

#### 5. Run Reference Generation Script
```bash
# Run the comprehensive reference generation for all tissues
bash analysis/src/cnv_calling/numbat/generate_all_tissue_references.test.sh

# Expected output:
# ===============================================
# NUMBAT Reference Generation for All Tissues
# ===============================================
# Start Time: [timestamp]
# 
# =========== Processing: lowseq_488B ===========
# [STEP 1] Verifying inputs for lowseq_488B...
#   ✓ Barcodes: lowseq_488B.barcodes.tsv
#   ✓ Fragments: lowseq_488B.fragments.sort.filtered.bed.gz
#   ✓ Bins: var220kb.rds
# 
# [STEP 2] Generating aggregated reference for lowseq_488B...
#   Output: lambdas_lowseq_488B_ATAC_bincnt.rds
# 
# [... similar for other tissues ...]
# 
# ===============================================
# ✓ Reference generation complete!
#   Completed: [timestamp]
# ===============================================
#
# Generated references:
# -rw-rw-r-- 1 preshita paxlab 178K ... lambdas_lowseq_488B_ATAC_bincnt.rds
# -rw-rw-r-- 1 preshita paxlab 175K ... lambdas_lowseq_489_ATAC_bincnt.rds
# -rw-rw-r-- 1 preshita paxlab 189K ... lambdas_deepseq_488B_ATAC_bincnt.rds
# -rw-rw-r-- 1 preshita paxlab 192K ... lambdas_deepseq_489_ATAC_bincnt.rds
```

#### 6. Monitor Progress (Optional)
If you need to monitor from another terminal:
```bash
# From another login node terminal, attach to tmux session
tmux attach -t numbat_refs

# To detach without stopping the process: Ctrl+B then D
```

#### 7. Verify Generated References
```bash
# After script completes, verify all references were created
ls -lh Data/04_analysis/cnv/numbat/reference/lambdas_*ATAC_bincnt.rds

# Should show:
# -rw-rw-r-- 1 preshita paxlab 178K May 18 XX:XX lambdas_lowseq_488B_ATAC_bincnt.rds
# -rw-rw-r-- 1 preshita paxlab 175K May 18 XX:XX lambdas_lowseq_489_ATAC_bincnt.rds
# -rw-rw-r-- 1 preshita paxlab 189K May 18 XX:XX lambdas_deepseq_488B_ATAC_bincnt.rds
# -rw-rw-r-- 1 preshita paxlab 192K May 18 XX:XX lambdas_deepseq_489_ATAC_bincnt.rds
```

#### 8. Check Generation Logs
```bash
# Each tissue gets its own log file
ls -lh Data/04_analysis/cnv/numbat/reference/*_reference_generation.log

# View logs to verify no errors
cat Data/04_analysis/cnv/numbat/reference/lowseq_488B_reference_generation.log

# Check for fragment matching counts
grep -i "fragment" Data/04_analysis/cnv/numbat/reference/*_reference_generation.log
# Should show values like: "103M fragments matched" (NOT 0)
```

#### 9. Exit tmux (Optional)
```bash
# After generation completes and verification is done:
exit

# Or detach: Ctrl+B then D (keeps session alive)
```

---

## What Each Reference File Contains

After running this script, you'll have these files in `Data/04_analysis/cnv/numbat/reference/`:

| File | Size | Tissue | Purpose |
|------|------|--------|---------|
| `lambdas_lowseq_488B_ATAC_bincnt.rds` | ~180K | Lowseq 488B | Aggregated ATAC reference for CNV baseline |
| `lambdas_lowseq_489_ATAC_bincnt.rds` | ~180K | Lowseq 489 | Aggregated ATAC reference for CNV baseline |
| `lambdas_deepseq_488B_ATAC_bincnt.rds` | ~190K | Deepseq 488B | Aggregated ATAC reference for CNV baseline |
| `lambdas_deepseq_489_ATAC_bincnt.rds` | ~190K | Deepseq 489 | Aggregated ATAC reference for CNV baseline |
| `lowseq_488B_reference_generation.log` | - | Lowseq 488B | Diagnostic log |
| `lowseq_489_reference_generation.log` | - | Lowseq 489 | Diagnostic log |
| `deepseq_488B_reference_generation.log` | - | Deepseq 488B | Diagnostic log |
| `deepseq_489_reference_generation.log` | - | Deepseq 489 | Diagnostic log |

---

## Common Issues & Troubleshooting

### Issue: "0 fragments matched"
**Problem**: Fragment barcode suffix mismatch
```bash
# Check what's in your fragment files
zcat Data/01_inputs/fragments/{TISSUE}/*.bed.gz | head -3 | cut -f4
# Should show: XXXXXX-1 (with -1 suffix)

# Check barcode file format
head Data/01_inputs/barcodes/tissue_barcodes/{TISSUE}/*.tsv | cut -f1
# Should show: XXXXXX (no suffix)
```
**Solution**: Use `get_binned_atac_fixed.R` (not the original script)

### Issue: "Barcode file not found"
```bash
# Verify barcode files exist
ls -lh Data/01_inputs/barcodes/tissue_barcodes/*/
# Should show subdirectories with .barcodes.tsv files
```

### Issue: "Fragment file not found"
```bash
# Verify fragment files exist and are compressed
ls -lh Data/01_inputs/fragments/*/
# Should show .bed.gz files (NOT uncompressed .bed)
```

### Issue: Script runs but creates tiny output files (<1MB)
```bash
# Check the log for error messages
tail -50 Data/04_analysis/cnv/numbat/reference/{TISSUE}_reference_generation.log

# Common causes:
# 1. No fragments matched (barcode format issue)
# 2. R library not loaded (check module load R)
# 3. Bin file not found (check NUMBAT_EXTDATA path)
```

---

## Next Steps After Reference Generation

Once all references are generated:

1. **Generate Pileup & Phasing** (if not already done):
   ```bash
   qsub analysis/src/cnv_calling/numbat/run_numbat_pileup_phase_lowseq_488B.qsub.sh
   qsub analysis/src/cnv_calling/numbat/run_numbat_pileup_phase_lowseq_489.qsub.sh
   qsub analysis/src/cnv_calling/numbat/run_numbat_pileup_phase_deepseq_488B.qsub.sh
   qsub analysis/src/cnv_calling/numbat/run_numbat_pileup_phase_deepseq_489.qsub.sh
   ```

2. **Run NUMBAT CNV Analysis**:
   ```bash
   qsub analysis/src/cnv_calling/numbat/run_numbat_analysis_lowseq_488B.qsub.sh
   qsub analysis/src/cnv_calling/numbat/run_numbat_analysis_lowseq_489.qsub.sh
   qsub analysis/src/cnv_calling/numbat/run_numbat_analysis_deepseq_488B.qsub.sh
   qsub analysis/src/cnv_calling/numbat/run_numbat_analysis_deepseq_489.qsub.sh
   ```

3. **Archive Test Files** (clean up after testing):
   ```bash
   mkdir -p Data/01_inputs/archive/test_files_archive_$(date +%Y%m%d)
   mv analysis/src/cnv_calling/numbat/generate_all_tissue_references.test.sh \
     Data/01_inputs/archive/test_files_archive_$(date +%Y%m%d)/
   ```

---

## Resource Requirements

- **CPU Cores**: 8 (multi-threaded binning & aggregation)
- **Memory**: 64 GB total (8 cores × 8 GB per core)
- **Wall Time**: 2.5 hours max (typically 30-40 min per tissue)
- **Disk Space**: ~200MB output + ~500MB temporary files

---

## Key Parameters Reference

| Parameter | Value | Description |
|-----------|-------|-------------|
| Fragment File | `.bed.gz` | Compressed BED format with barcodes |
| Barcode List | `.tsv` (1 per line) | Raw barcode format (no -1 suffix) |
| Genomic Bins | `var220kb.rds` | 220 kilobase windows across genome |
| Output Matrix | `adata_atac.rds` | Cell × bin matrix for NUMBAT |
| Output Reference | `lambdas_*ATAC_bincnt.rds` | Aggregated baseline for CNV detection |

---

## Contact & Documentation

- **NUMBAT Paper**: Jiang et al. 2024
- **NUMBAT GitHub**: kharchenkolab/numbat
- **NUMBAT Vignette**: https://kharchenkolab.github.io/numbat/

For questions specific to the barcode format fix, see section "CRITICAL ISSUE: Barcode Format Compatibility" in copilot-instructions.md
