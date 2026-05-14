# Spatial ATAC Project - Reorganization Complete (2026-05-13)

## Quick Navigation

👉 **For Data Organization**: [Data/ORGANIZATION_SUMMARY.md](Data/ORGANIZATION_SUMMARY.md)
👉 **For Analysis Pipelines**: [analysis/ORGANIZATION_SUMMARY.md](analysis/ORGANIZATION_SUMMARY.md)
👉 **For Job Logs**: [analysis/qsub_logs/README.md](analysis/qsub_logs/README.md)
👉 **For Plots**: [analysis/plots/README.md](analysis/plots/README.md)

---

## 📋 What Was Reorganized

### ✅ Data Folder
- **Barcodes**: Organized by object in `Data/01_inputs/barcodes/tissue_barcodes/`
- **Fragments**: Organized by object in `Data/01_inputs/fragments/` with .tbi indices
- **BAM files**: Symlinked in `Data/01_inputs/bam/` with `.lnk` suffix
- **Documentation**: Comprehensive symlink mapping in `Data/file_mapping/`

**Key Metrics**:
| Object | Barcodes Before | After | % Retained |
|--------|-----------------|-------|-----------|
| deepseq_488B | 11,645 | 11,467 | 98.5% |
| deepseq_489 | 4,671 | 4,622 | 98.9% |
| lowseq_488B | 11,645 | 11,612 | 99.7% |
| lowseq_489 | 4,671 | 4,622 | 98.9% |

### ✅ Analysis QSub Logs
- **build_tissue/**: 9 files - Edge-effect filtering jobs
- **diagnostics/**: 14 files - Variant QC and somatic tests
- **archived/**: 12 files - Historical debug runs

### ✅ Analysis Plots
- **edge_effect/**: 12 PNG plots (4 samples × 3 plot types) - ✅ COMPLETE
- **cnv_analysis/**: Ready for Alleloscope and NUMBAT outputs
- **comparison/**: Organized into variant_qc and somatic_characterization

---

## 🚀 Next Phase - Ready for Analysis

### Alleloscope (3 runs planned)
```bash
# Input paths for your scripts:
barcodes <- "Data/01_inputs/barcodes/tissue_barcodes/{object}/{object}.no_edge_effect.barcodes.tsv"
fragments <- "Data/01_inputs/fragments/{object}/{object}.fragments.sort.filtered.bed.gz"

# Output directory:
analysis/plots/cnv_analysis/alleloscope/
```

**Runs planned**:
1. Tissue 488B (deepseq_488B + lowseq_488B)
2. Tissue 489 (deepseq_489 + lowseq_489) ⚠️ *User notes: "tissue 489 is giving me grief"*
3. Combined tissues

### Somatic SNV Comparison
```bash
# Input paths:
bam_deepseq <- "Data/01_inputs/bam/deepseq_{object}.bam.lnk"
bam_lowseq <- "Data/01_inputs/bam/lowseq_{object}.bam.lnk"
fragments <- "Data/01_inputs/fragments/{object}/*.bed.gz"

# Output directory:
analysis/plots/comparison/somatic_characterization/
```

**Purpose**: Compare somatic variants between lowseq and deepseq (should be identical)

---

## 📂 Folder Structure at a Glance

```
Data/01_inputs/
├── archive/              # BAM files archived here (206GB)
├── bam/                  # Symlinks to BAM files (.lnk suffix)
├── barcodes/             # Barcode files organized by object
│   └── tissue_barcodes/
│       ├── deepseq_488B/
│       ├── deepseq_489/
│       ├── lowseq_488B/
│       └── lowseq_489/
├── fragments/            # Fragment files organized by object
│   ├── deepseq_488B/     # *.bed.gz + *.tbi
│   ├── deepseq_489/
│   ├── lowseq_488B/
│   └── lowseq_489/
├── file_mapping/         # Symlink documentation & audit trail
└── [other folders]

analysis/
├── qsub_logs/
│   ├── build_tissue/     # 9 files - Edge-effect filtering
│   ├── diagnostics/      # 14 files - Tests & comparisons
│   ├── archived/         # 12 files - Old debug runs
│   └── README.md         # Job log reference
├── plots/
│   ├── edge_effect/      # 12 PNG plots (4 samples × 3 types)
│   ├── cnv_analysis/     # Ready for Alleloscope/NUMBAT
│   ├── comparison/
│   │   ├── variant_qc/
│   │   ├── somatic_characterization/
│   │   └── somatic_comparison_old/
│   └── README.md         # Plot organization guide
└── ORGANIZATION_SUMMARY.md
```

---

## 📖 Documentation Created

| File | Size | Purpose |
|------|------|---------|
| Data/ORGANIZATION_SUMMARY.md | 6.6K | Complete data folder reference |
| analysis/ORGANIZATION_SUMMARY.md | 9.0K | Analysis pipeline & structure |
| analysis/plots/README.md | 3.6K | Plot organization & naming |
| analysis/qsub_logs/README.md | 5.5K | Job log categories & tracking |
| Data/file_mapping/* | 57KB | Symlink audit trail (BEFORE/AFTER) |

---

## ✨ Key Improvements

✅ **Before**: Files scattered across multiple locations
✅ **After**: Organized by data type and object, with clear naming

✅ **Before**: No barcode filtering applied
✅ **After**: Edge-effect filtered barcodes available (98%+ retention)

✅ **Before**: Job logs mixed in root folder
✅ **After**: Categorized by pipeline stage (build_tissue, diagnostics, archived)

✅ **Before**: Unclear symlink management
✅ **After**: Documented with BEFORE/AFTER mapping files

✅ **Before**: Minimal documentation
✅ **After**: 5 comprehensive markdown guides

---

## 🔧 Useful Commands

### View Documentation
```bash
# Read data organization guide
cat Data/ORGANIZATION_SUMMARY.md

# Read analysis guide
cat analysis/ORGANIZATION_SUMMARY.md

# Read job log guide
cat analysis/qsub_logs/README.md
```

### Check Job Logs
```bash
# View latest successful job
tail -20 analysis/qsub_logs/build_tissue/build_tissue.o5612496

# Check for errors
cat analysis/qsub_logs/build_tissue/build_tissue.e5612496

# Monitor queue
qstat -u preshita
```

### Verify Input Files
```bash
# Check barcodes
ls -lh Data/01_inputs/barcodes/tissue_barcodes/deepseq_488B/

# Check fragments
ls -lh Data/01_inputs/fragments/deepseq_488B/

# Check BAM symlinks
ls -lh Data/01_inputs/bam/
```

---

## ⚠️ Known Issues & Monitoring

1. **Tissue 489 compatibility** (user note)
   - May have Alleloscope issues
   - Keep detailed logs for debugging
   - Compare barcode format with tissue 488B

2. **Somatic variant comparison**
   - Deepseq and lowseq should have identical variants
   - Document any discrepancies
   - Check variant counts match

---

## 🎯 Status Summary

| Task | Status | Details |
|------|--------|---------|
| Data organization | ✅ Complete | 4 objects, barcodes+fragments+BAM |
| Barcode filtering | ✅ Complete | 98-99.7% retention, edge effects removed |
| Job log categorization | ✅ Complete | 35 files organized into 3 categories |
| Plot organization | ✅ Complete | 4 main categories, 12 edge plots ready |
| Documentation | ✅ Complete | 5 markdown guides + symlink audit trail |
| Alleloscope prep | ✅ Ready | Input paths verified, output dirs ready |
| Somatic comp prep | ✅ Ready | BAM symlinks verified, output dirs ready |

---

## 🚀 Next Actions

1. **Begin Alleloscope runs** (3 parallel jobs):
   - Monitor tissue 489 carefully
   - Output to `analysis/plots/cnv_analysis/alleloscope/`

2. **Run somatic comparison**:
   - Compare deepseq vs lowseq variants
   - Output to `analysis/plots/comparison/somatic_characterization/`

3. **Archive old results** (as needed):
   ```bash
   mv analysis/plots/comparison analysis/plots/comparison_$(date +%Y%m%d)
   ```

---

**Last Updated**: 2026-05-13 22:07 UTC
**Status**: ✅ **ALL REORGANIZATION TASKS COMPLETE**
**Next Phase**: Ready for Alleloscope & Somatic Analysis

