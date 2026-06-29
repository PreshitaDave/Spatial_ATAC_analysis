# Arrow Binsize Test Setup - Verification Checklist

## Setup Completed ✓

All files and directories for the arrow binsize testing workflow have been created and verified.

### Created Files

**Main Scripts:**
- ✓ `analysis/src/build_tissue/binsize_test_archr/create_arrow_variants.R` (14 KB)
  - Creates arrow files with specific tile sizes and binarization options
  - Computes sparsity metrics
  
- ✓ `analysis/src/build_tissue/binsize_test_archr/compare_arrow_sparsity.R` (10 KB)
  - Analyzes and compares sparsity metrics across all combinations
  - Generates visualization plots and recommendations
  
- ✓ `analysis/src/build_tissue/binsize_test_archr/submit_all_arrow_variants.sh` (3.5 KB)
  - Master job submission script
  - Submits all 16 combinations in parallel

**Documentation:**
- ✓ `WORKFLOW.md` - Complete workflow documentation with examples
- ✓ `README.md` - Quick start guide
- ✓ `SETUP_CHECKLIST.md` - This file

### Created Directories

- ✓ `Data/01_inputs/arrow/arrow_binarize/` - For binarized arrow files
- ✓ `Data/01_inputs/arrow/arrow_not_binarize/` - For non-binarized arrow files
- ✓ `analysis/binsize_comparison/` - For comparison results (created at first analysis run)

### Available Input Files

**Fragment Files:** ✓
- deepseq_488B
- deepseq_489
- deepseq_combined
- lowseq_488B
- lowseq_489
- lowseq_combined

**Input Arrow Files:** ✓
- Deepseq_488B.arrow (12 GB)
- Deepseq_489.arrow (3.5 GB)
- Lowseq_488B.arrow (4.7 GB)
- Lowseq_489.arrow (1.5 GB)

**Barcode Files:** ✓
Available for: deepseq_488B, deepseq_489, lowseq_488B, lowseq_489

## Quick Start Commands

### 1. Submit all jobs (16 parallel jobs)
```bash
cd /projectnb/paxlab/presh/projects/spatial_atac
bash analysis/src/build_tissue/binsize_test_archr/submit_all_arrow_variants.sh
```

### 2. Monitor progress
```bash
qstat -u preshita | grep create_arrow
tail -f analysis/qsub_logs/create_arrow_*.log
```

### 3. After all jobs complete (~1 hour), run analysis
```bash
module load R
Rscript analysis/src/build_tissue/binsize_test_archr/compare_arrow_sparsity.R
```

### 4. View results
```bash
cat analysis/binsize_comparison/sparsity_recommendations.txt
```

## Testing Setup (Optional)

### Test with single tissue first
```bash
# Allocate compute resources
qrsh -l h_rt=01:00:00 -pe omp 8 -P paxlab -l mem_per_core=8G

# Inside qrsh:
module load R
cd /projectnb/paxlab/presh/projects/spatial_atac

# Create one variant (500bp, not binarized)
Rscript analysis/src/build_tissue/binsize_test_archr/create_arrow_variants.R deepseq_488B 500 FALSE

# Check metrics
cat analysis/binsize_comparison/deepseq_488B_500bp_binarizeFALSE_metrics.txt
```

## Test Matrix

The full test will create and analyze:

| Tissue | Tile Sizes | Binarize Options | Combinations |
|--------|-----------|-----------------|--------------|
| deepseq_488B | 500bp, 5000bp | FALSE, TRUE | 4 jobs |
| deepseq_489 | 500bp, 5000bp | FALSE, TRUE | 4 jobs |
| lowseq_488B | 500bp, 5000bp | FALSE, TRUE | 4 jobs |
| lowseq_489 | 500bp, 5000bp | FALSE, TRUE | 4 jobs |
| **Total** | | | **16 jobs** |

## Expected Outputs

After running all jobs and analysis:

```
analysis/binsize_comparison/
├── sparsity_comparison_table.csv      # Full metrics
├── sparsity_comparison_plots.pdf      # Visualization
├── sparsity_recommendations.txt       # Summary with best params
└── {tissue}_{tilesize}bp_binarize_{TRUE,FALSE}_metrics.txt  # Per-combination metrics
```

## Key Metrics

For each combination, the following will be computed:
- **Sparsity**: Proportion of zeros in tile matrix (0.0-1.0)
- **Density**: Proportion of non-zeros (1 - sparsity)
- **Cell Coverage**: Average tiles per cell
- **Tile Coverage**: Average cells per tile

## File Locations Reference

| Type | Location |
|------|----------|
| Fragment files | `Data/01_inputs/fragments/` |
| Input arrows | `Data/01_inputs/arrow/` |
| Output arrows | `Data/01_inputs/arrow/arrow_{binarize,not_binarize}/` |
| R scripts | `analysis/src/build_tissue/binsize_test_archr/` |
| Results | `analysis/binsize_comparison/` |
| Job logs | `analysis/qsub_logs/create_arrow_*.{log,err}` |

## Resource Allocation

Each job uses:
- **CPU cores**: 8 (via `-pe omp 8`)
- **Memory**: 32 GB total (8G per core via `-l mem_per_core=8G`)
- **Time limit**: 4 hours (via `-l h_rt=04:00:00`)
- **Queue**: paxlab (via `-P paxlab`)

Typical execution time:
- Single job: 20-40 minutes
- All 16 jobs parallel: 40-60 minutes
- Comparison analysis: 10-15 minutes

## Next Steps After Completion

1. **Review recommendations** in `sparsity_recommendations.txt`
2. **Check plots** in `sparsity_comparison_plots.pdf`
3. **Update downstream scripts** to use optimal parameters
4. **Archive results** if satisfied with findings
5. **Document findings** in analysis notes

## Troubleshooting

### Files not found
```bash
# Verify all inputs exist
ls /projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/fragments/*/
ls /projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/arrow/*.arrow
```

### Script won't run
```bash
# Verify R installation on compute node
qrsh -l h_rt=00:30:00 -pe omp 2 -P paxlab
which Rscript && Rscript --version
module load R
which Rscript && Rscript --version
```

### Job submission fails
```bash
# Check if on compute node (should NOT be scc1)
hostname

# If on login node, allocate compute:
qrsh -l h_rt=01:00:00 -pe omp 4 -P paxlab
```

## Important Notes

1. **Run from appropriate location**: Must be on compute node or qsub system, NOT login node
2. **Module loading**: Scripts will load R module automatically in qsub jobs
3. **Data persistence**: Output directories will accumulate results; clean if rerunning
4. **Parallel execution**: All jobs can run simultaneously (no dependencies between them)
5. **QC filtering**: Analysis uses TSS ≥ 3 and nFrags ≥ 1000 filters consistent with main pipeline

## Script Parameters

### create_arrow_variants.R
```bash
Usage: Rscript create_arrow_variants.R <tissue> <tilesize> <binarize>

tissue:    deepseq_488B, deepseq_489, lowseq_488B, lowseq_489, deepseq_combined, lowseq_combined
tilesize:  500 or 5000
binarize:  FALSE or TRUE
```

### compare_arrow_sparsity.R
```bash
Usage: Rscript compare_arrow_sparsity.R [tissue_filter]

tissue_filter: (optional) Filter to single tissue, or omit for all tissues
```

### submit_all_arrow_variants.sh
```bash
Usage: bash submit_all_arrow_variants.sh

No parameters needed. Submits all 16 combinations automatically.
Modify TISSUES array in script to test different tissues.
```

---

**Setup completed**: 2026-06-21  
**Status**: Ready for job submission ✓  
**Next step**: Run `bash analysis/src/build_tissue/binsize_test_archr/submit_all_arrow_variants.sh`
