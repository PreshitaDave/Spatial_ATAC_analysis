# Arrow Binsize Test Workflow
## Comparing Tile Sizes and Binarization Options for ArchR Performance

**STATUS: Arrow files (16 variants) already exist and are valid. This workflow focuses on metrics computation and comparison.**

This workflow tests different combinations of tile sizes (500bp vs 5000bp) and binarization options (TRUE vs FALSE) to determine which parameters produce the best ArchR performance based on sparsity and coverage metrics.

### Overview

The workflow has three phases:

1. **Creation Phase** (DONE): Arrow files (16 variants × 4 tissues) already created from fragment files
   - Each arrow was created with binarization baked into the storage layer
2. **Metrics Phase** (CURRENT): Recompute metrics for all variants from existing arrow files
   - No arrow creation needed; only metrics computation
3. **Analysis Phase**: Use the best-performing parameters for downstream ArchR analysis

### Directory Structure

```
analysis/src/build_tissue/binsize_test_archr/
├── create_arrow_variants.R          # Main R script to create arrow variants
├── compare_arrow_sparsity.R         # Analysis script to compare sparsity metrics
├── submit_all_arrow_variants.sh     # Master submission script for all jobs
├── WORKFLOW.md                      # This file
└── README.md                        # Quick reference guide

Data/01_inputs/arrow/
├── arrow_binarize/                  # Output: binarized arrow files
├── arrow_not_binarize/              # Output: non-binarized arrow files
└── [standard arrow files]           # Input: existing arrow files

analysis/binsize_comparison/         # Output: metrics and plots
├── sparsity_comparison_table.csv    # Results table
├── sparsity_comparison_plots.pdf    # Visualization
└── sparsity_recommendations.txt     # Summary and recommendations
```

### Phase 1: Arrow Creation (ALREADY COMPLETED)

Arrow files exist in:
- `Data/01_inputs/arrow/arrow_binarize/` (8 files, binarized TileMatrix)
- `Data/01_inputs/arrow/arrow_not_binarize/` (8 files, count-based TileMatrix)

Each file is ~1.3–6.7GB and contains all fragments from the original files for one tissue,
with TileMatrix baked in at storage time using the specified tilesize and binarization.

### Phase 2: Metrics Computation and Aggregation (CURRENT)

#### Step 2a: Submit all jobs at once (recommended)

```bash
cd /projectnb/paxlab/presh/projects/spatial_atac
bash analysis/src/build_tissue/binsize_test_archr/submit_metrics_recompute.sh
```

This will submit:
- 16 parallel metrics-recomputation jobs (one per tissue × tilesize × binarize)
  - Each loads one existing arrow file and computes sparsity, density, and coverage metrics
  - Logs to `analysis/qsub_logs/recompute_metrics_*.{log,err}`
  - ~30–60 seconds per job (no arrow creation, just metric computation)
  - Requires only 4 cores, 8GB memory, 1-hour time limit (rightcores-sized for this step)
- 1 aggregation job (auto-submitted with `-hold_jid` dependency on all 16)
  - Reads the 16 metrics files, aggregates into CSV
  - Generates comparison plots and recommendations
  - Logs to `analysis/qsub_logs/aggregate_metrics.{log,err}`
  - ~1-2 minutes total

#### Step 2b: Monitor progress

```bash
# View all jobs
qstat -u preshita

# Watch metrics jobs
tail -f analysis/qsub_logs/recompute_metrics_*.log

# Watch aggregation (starts after metrics jobs finish)
tail -f analysis/qsub_logs/aggregate_metrics.log
```

#### Step 2c: Review output files (after all jobs complete)

```bash
# Check metrics table
cat analysis/binsize_comparison/sparsity_comparison_table.csv

# View recommendations
cat analysis/binsize_comparison/sparsity_recommendations.txt

# View plots
# Download: analysis/binsize_comparison/sparsity_comparison_plots.pdf
```

#### Step 2d: Manual submission (if needed)

For a single tissue or debugging:

```bash
# Test one combination
qrsh -l h_rt=00:30:00 -pe omp 4 -P paxlab -l mem_per_core=8G

# Inside qrsh session:
module load R
cd /projectnb/paxlab/presh/projects/spatial_atac
Rscript analysis/src/build_tissue/binsize_test_archr/create_arrow_variants.R deepseq_488B 500 FALSE
```

Or to manually run the aggregation after metrics are ready:

```bash
module load R
Rscript analysis/src/build_tissue/binsize_test_archr/compare_arrow_sparsity.R

# Or for a specific tissue:
Rscript analysis/src/build_tissue/binsize_test_archr/compare_arrow_sparsity.R deepseq_488B
```

### Phase 3: Use Best Parameters Downstream

Once you've reviewed the metrics and recommendations:

1. **Identify optimal parameters** from `sparsity_recommendations.txt`
   - Look for per-tissue best configurations by sparsity/density
   - Check "tiles per cell" values — aim for ≥100 for robust per-spot CNV calling

2. **Update downstream scripts** (e.g., `0_create_archr_no_doublet.R`, CNV-calling scripts):
   ```R
   # Example: if best is 500bp, binarize=TRUE for deepseq_488B
   tilesize_by_tissue <- list(
     deepseq_488B = 500,
     deepseq_489 = 5000,
     lowseq_488B = 500,
     lowseq_489 = 500
   )
   binarize_by_tissue <- list(
     deepseq_488B = TRUE,
     deepseq_489 = FALSE,
     lowseq_488B = FALSE,
     lowseq_489 = TRUE
   )

   # In your script:
   proj <- addTileMatrix(proj, force = TRUE, tileSize = tilesize_by_tissue[[tissue]])
   ```

3. **For production CNV analysis**:
   Create separate ArchR projects per tissue using the optimal parameters, or use one unified
   parameter set for multi-tissue consistency (may require a compromise if per-tissue optima differ)

### Script Usage Details

#### create_arrow_variants.R

Creates arrow files with specific parameters and computes sparsity metrics.

**Usage:**
```bash
Rscript create_arrow_variants.R <tissue> <tilesize> <binarize>
```

**Parameters:**
- `tissue`: One of {deepseq_488B, deepseq_489, lowseq_488B, lowseq_489, deepseq_combined, lowseq_combined}
- `tilesize`: 500 or 5000 (in base pairs)
- `binarize`: TRUE or FALSE (as string)

**Example:**
```bash
Rscript create_arrow_variants.R lowseq_489 5000 TRUE
```

**Outputs:**
- Arrow files: `Data/01_inputs/arrow/arrow_{binarize,not_binarize}/{tissue}_{tilesize}bp.arrow`
- Metrics: `analysis/binsize_comparison/{tissue}_{tilesize}bp_binarize{TRUE,FALSE}_metrics.txt`

**What it does:**
1. Loads input arrow file
2. Creates ArchR project with specified tile size
3. Applies QC filtering (TSS ≥ 3, nFrags ≥ 1000)
4. Computes tile matrix with specified binarization
5. Calculates sparsity, density, and coverage metrics
6. Saves metrics to file
7. Optionally saves modified arrow file to output directory

#### compare_arrow_sparsity.R

Compares sparsity metrics across all tile size and binarization combinations.

**Usage:**
```bash
Rscript compare_arrow_sparsity.R [tissue_name]
```

**Parameters:**
- `tissue_name` (optional): Filter to single tissue. If omitted, analyzes all tissues

**Example:**
```bash
# Analyze all tissues
Rscript compare_arrow_sparsity.R

# Analyze specific tissue
Rscript compare_arrow_sparsity.R deepseq_488B
```

**Outputs:**
- `sparsity_comparison_table.csv`: Full metrics table for all combinations
- `sparsity_comparison_plots.pdf`: Visualization comparing sparsity, density, and coverage
- `sparsity_recommendations.txt`: Summary with best parameters for each tissue

**Metrics included:**
- **Sparsity**: Proportion of zero-valued tiles (0.0-1.0)
- **Density**: Proportion of non-zero tiles (1-sparsity)
- **Cell Coverage**: Average number of tiles per cell (higher = more information)
- **Tile Coverage**: Average number of cells per tile (higher = more usage)

### Performance Interpretation

Lower sparsity = more information retained = better for downstream analysis

**Typical findings:**
- Smaller tile sizes (500bp) often have HIGHER sparsity (more zeros) but finer resolution
- Larger tile sizes (5000bp) often have LOWER sparsity (fewer zeros) but coarser resolution
- Binarization may affect sparsity depending on fragment distribution
- Tissue-specific effects: deepseq may differ from lowseq

### Troubleshooting

**Issue: "Arrow file not found"**
- Verify input arrow files exist: `ls -lh Data/01_inputs/arrow/`
- Check tissue names match: `deepseq_488B`, not `Deepseq_488B`

**Issue: "No matching cells found"**
- May indicate barcode format mismatch
- Check barcode files: `head Data/01_inputs/barcodes/tissue_barcodes/*/`

**Issue: Jobs stuck/failing**
- Check logs: `tail -f analysis/qsub_logs/create_arrow_*.err`
- Verify R modules load: `module load R && which Rscript`
- Test on compute node first: `qrsh -l h_rt=01:00:00 -pe omp 8 ...`

**Issue: Comparison script fails to find metrics**
- Ensure creation jobs completed successfully
- Check if metric files exist: `ls -lh analysis/binsize_comparison/`
- Manually run single job to test: `Rscript create_arrow_variants.R deepseq_488B 500 FALSE`

### Tips for Optimization

1. **Testing Strategy**: 
   - Start with one tissue/combination to verify setup works
   - Then submit full batch of 16 jobs
   - Run comparison only after all jobs complete

2. **Time Management**:
   - Each job takes ~20-40 minutes on 8 cores
   - Full batch (16 jobs parallel) takes ~40-60 minutes total
   - Comparison takes ~10-15 minutes

3. **Resource Management**:
   - Jobs request 8 cores and 32GB total (8G per core)
   - Monitor with `qstat` to ensure jobs aren't queued long
   - If queue is full, submit jobs in smaller batches (e.g., 4 at a time)

4. **Iterative Testing**:
   - If you want to test different parameters, modify the TISSUES, TILESIZES, or BINARIZE_OPTIONS lists in `submit_all_arrow_variants.sh`
   - Can also test combined tissues if desired

### Expected Output Example

```
====================================
Arrow Sparsity Comparison Summary
====================================

Overall Best Density Configuration:
  Tissue: deepseq_488B
  Tile Size: 5000 bp
  Binarize: TRUE
  Density: 0.087543
  Sparsity: 91.25%

Per-Tissue Recommendations:
  deepseq_488B: Best density=0.087543 with 5000 bp tiles, binarize=TRUE
  deepseq_489: Best density=0.105432 with 5000 bp tiles, binarize=FALSE
  lowseq_488B: Best density=0.042109 with 500 bp tiles, binarize=FALSE
  lowseq_489: Best density=0.063891 with 5000 bp tiles, binarize=TRUE
```

### File Organization After Completion

```
Data/01_inputs/arrow/
├── arrow_binarize/              # Binarized variants
│   ├── deepseq_488B_500bp.arrow
│   ├── deepseq_488B_5000bp.arrow
│   ├── deepseq_489_500bp.arrow
│   └── ...
├── arrow_not_binarize/          # Non-binarized variants
│   ├── deepseq_488B_500bp.arrow
│   ├── deepseq_488B_5000bp.arrow
│   └── ...
└── [original arrow files]       # Unchanged

analysis/binsize_comparison/
├── sparsity_comparison_table.csv
├── sparsity_comparison_plots.pdf
└── sparsity_recommendations.txt
```

### Phase 4: UMAP/Gene-Score Validation (Empirical Assessment)

After Phase 2 sparsity metrics, validate which tile-size/binarize parameters actually produce
reliable downstream ArchR results (quality cluster separation, stable gene-score patterns) by
running the full ArchR pipeline on each variant.

#### Step 4a: Submit UMAP/gene-score jobs (pilot on one tissue first)

```bash
cd /projectnb/paxlab/presh/projects/spatial_atac
bash analysis/src/build_tissue/binsize_test_archr/submit_archr_variant_pilot.sh
```

This submits:
- 4 parallel ArchR pipeline jobs (lowseq_489 × {500bp, 5000bp} × {binarize TRUE, FALSE})
  - Each runs: addIterativeLSI → addClusters → addUMAP → generate plots
  - Outputs UMAP+gene-score CSV and per-variant PDF with cluster + marker-gene plots
  - Uses house parameters (LSI: iterations=2, varFeatures=25000, dimsToUse=1:30;
    Clusters: method=Seurat, resolution=0.8; UMAP: nNeighbors=30, minDist=0.5, metric=cosine)
  - Logs to `analysis/qsub_logs/archr_variant_*.{log,err}`
  - ~30-60 minutes per job (iterative LSI is the bottleneck)
- 1 aggregation job (depends on all 4)
  - Generates faceted comparison plots: UMAP by cluster across all 4 variants, then per-marker-gene
  - Creates `archr_umap_genescore_comparison.pdf` with side-by-side variant assessment
  - Logs to `analysis/qsub_logs/archr_variant_aggregator_pilot.log`

#### Step 4b: Review comparison plots

```bash
# View the full comparison (all 4 variants side by side)
# analysis/binsize_comparison/archr_umap_genescore_comparison.pdf
#
# Key observations to make:
# 1. Cluster separation: Are clusters similarly separated across all variants?
# 2. Marker specificity: Do markers show expected patterns consistently?
# 3. Consistency: Which tile-size/binarize combo looks most robust?
# 4. Artifacts: Any noisy/unusual embedding patterns in particular variants?
```

#### Step 4c: Expand to all tissues (after pilot validation)

Once you're confident the pilot looks good, modify `submit_archr_variant_pilot.sh` line ~28
to include all tissues:

```bash
TISSUES=("deepseq_488B" "deepseq_489" "lowseq_488B" "lowseq_489")  # vs just ("lowseq_489")
```

Then resubmit — same script will generate 16 jobs instead of 4 and an expanded comparison PDF.

#### Step 4d: Document findings

Create a summary comparing Phase 2 (sparsity metrics alone) vs Phase 4 (empirical UMAP/gene-
scores). Do the sparsity-based rankings hold up when you look at actual cluster quality and
marker expression patterns? Often they do, but checking is critical before committing to a
parameter choice.

### Next Steps After Comparison

Once you've identified optimal parameters (combining Phase 2 sparsity metrics + Phase 4 UMAP validation):

1. **Update downstream scripts** to use best tile size and binarization
2. **Archive this directory** if you're confident with the results
3. **Document the findings** in your analysis notes (include both metrics + visual evidence)
4. **Apply to new tissues** if needed using the same workflow
5. **Consider replicating** the analysis if parameters unexpectedly differ between runs

### References

These citations support the parameter choices and methodology used in this workflow:

1. **Granja JM, Corces MR, Pierce SE, et al.** "ArchR is a scalable software package for integrative
   single-cell chromatin accessibility analysis." *Nat Genet.* 2021;53(3):403–411.
   - Describes ArchR's default 500bp TileMatrix and iterative LSI approach for dimensionality reduction

2. **Cusanovich DA, Hill AJ, Aghamirzaie D, et al.** "A Single-Cell Atlas of In Vivo Mammalian
   Chromatin Accessibility." *Cell.* 2018;174(5):1309–1324.
   - Genome-wide fixed-bin matrices; precedent for exploring larger bin sizes (e.g., 5kb) in
     scATAC dimensionality reduction without sacrificing clustering quality

3. **Fang R, Preissl S, Li Y, et al.** "Comprehensive analysis of single cell ATAC-seq data with
   SnapATAC." *Nat Commun.* 2021;12:1337.
   - Alternative pipeline that defaults to 5kb genome binning, demonstrates viability of larger
     tile sizes for sparse ATAC data

4. **Chen H, Lareau C, Andreani T, et al.** "Assessment of computational methods for the analysis
   of single-cell ATAC-seq data." *Genome Biol.* 2019;20:241.
   - Benchmark comparing scATAC analysis methods; discusses binning strategy, binarization, and
     their effects on clustering accuracy and downstream performance

### Additional Notes

- This workflow uses QC filtering (TSS ≥ 3, nFrags ≥ 1000) consistent with `0_create_archr_no_doublet.R`
- Metrics are computed on post-QC cells
- Arrow files contain unfiltered fragments; filtering happens during ArchR analysis
- Sparsity comparison focuses on tile matrix properties, not other ArchR metrics
- Phase 4 (UMAP/gene-score validation) goes beyond bulk sparsity to empirically assess
  downstream analysis quality — **this is critical** for choosing parameters, since two
  configurations with similar sparsity can produce different clustering/embedding quality
- For spatial ATAC with 1-10 cells per spot, aim for ≥100 tiles per cell on average to ensure
  per-spot CNV calling has enough resolution (though CNV pipelines build their own bins, the
  TileMatrix choice still impacts ArchR-based metadata and clustering used downstream)

---

*Last updated: 2026-06-22*
