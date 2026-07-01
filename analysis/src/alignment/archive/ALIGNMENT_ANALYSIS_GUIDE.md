# MOSAICField Alignment Analysis - Summary & Workflow Guide

## Date: June 8, 2026

## Overview

This document summarizes the MOSAICField alignment analysis and the updated R markdown workflow for mapping ATAC spots to Xenium cells.

### Key Files Generated

1. **Updated R Markdown Notebook**: `analysis/src/alignment/cell_to_spot_map.Rmd`
   - 588 lines with comprehensive sections for debugging and mapping
   - Includes diagnostic plots, validation checks, and QC metrics

2. **Debug R Script**: `analysis/src/alignment/debug_alignment.test.R`
   - Validates MOSAICField alignment quality
   - Creates diagnostic overlay plots
   - Checks spatial coordinate consistency

3. **Debug qsub Wrapper**: `analysis/src/alignment/debug_alignment.qsub.sh`
   - Submits debug script to compute cluster
   - Includes module initialization and output redirection
   - Job ID: 5985840 (status: running on scc-tc1)

### MOSAICField Output Files

Location: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/alignment/mosaicfield_outputs/`

- `atac_mosaicfield_outputs.h5ad` (12 GB) - Aligned ATAC spots
- `xenium_mosaicfield_outputs.h5ad` (87 MB) - Aligned Xenium cells
- MOSAICField visualization images:
  - `01_cropped_views.png` - Before/after alignment (zoomed)
  - `02_full_views.png` - Before/after alignment (full)
  - `03_overlay_before_full.png` - Overlay before alignment
  - `04_overlay_after_full.png` - Overlay after alignment

## Updated R Markdown Notebook Structure

The `cell_to_spot_map.Rmd` now includes:

### Section 1: Environment Setup
- Python/R configuration
- Library loading (Giotto, Scanpy, Arrow, SF, etc.)
- Path initialization

### Section 2: Debugging - Alignment Quality Validation
**Purpose**: Validate that MOSAICField correctly aligned the ATAC and Xenium data

Subsections:
1. **Load Objects** - Import aligned H5AD files via Giotto
2. **Validate Coordinates** - Check spatial dimension ranges
3. **Alignment Quality Assessment** - Compare bounding boxes, area ratios
4. **Diagnostic Overlay Plots** - Generate 3 comparison visualizations

Generated Plots:
- `01_debug_overlay_all.png` - ATAC spots (red) + Xenium cells (blue)
- `02_debug_atac_spots.png` - ATAC spot distribution
- `03_debug_xenium_cells.png` - Xenium cell distribution

Expected Outcomes:
- ✓ PASS: Objects have similar spatial extents (area ratio 0.5-2.0)
- ⚠ CAUTION: Different scales (ratio 0.1-10)
- ✗ FAIL: Very different scales (ratio <0.1 or >10)

### Section 3: ATAC-to-Xenium Mapping (Polygon-Based)
**Purpose**: Create geometric mapping between ATAC spots and Xenium cells

Subsections:
1. **Create 10µm Polygons** - Define square polygons around ATAC spot centroids
   - 5µm half-width = 10µm total edge length
   - 4 vertices per spot define bounding box

2. **Load Xenium Boundaries** - Read cell boundaries from parquet file
   - Source: `/projectnb/paxlab/DATA/DriesSpatial/Xenium/output-XETG00253__0048833__488B__20241217__182301/cell_boundaries.parquet`
   - Creates SF polygons for each Xenium cell

3. **Add Polygons to Objects** - Register both polygon sets with Giotto
   - ATAC bins: 10µm squares
   - Xenium cells: Original cell boundaries

4. **Calculate Overlap** - Use Giotto's `calculateOverlap()` function
   - Maps each Xenium cell to containing ATAC spot
   - Handles multiple cells per spot
   - Handles unmapped cells (outside ATAC region)

5. **Extract & Analyze Results** - Generate mapping statistics
   - Total mapped/unmapped cells
   - Occupied/empty ATAC spots
   - Cell density distribution

6. **Visualize Results** - Create publication-quality plots
   - `04_mapping_density_histogram.png` - Histogram of cells/spot
   - `05_mapping_spatial_distribution.png` - Spatial map with mapping status

7. **Save Results** - Export data to CSV
   - `xenium_to_atac_mapping.csv` - Detailed per-cell mapping
   - `atac_spot_density_stats.csv` - Per-spot statistics
   - `mapping_summary_report.csv` - Summary metrics

### Section 4: Validation Checks

Three automated QC checks:

1. **Mapping Coverage Check**
   - ✓ PASS: >90% of Xenium cells mapped
   - ⚠ CAUTION: 50-90% mapped
   - ✗ FAIL: <50% mapped (alignment issue)

2. **Spot Occupancy Check**
   - ✓ PASS: >70% of ATAC spots have ≥1 Xenium cell
   - ⚠ CAUTION: 30-70% occupied
   - ✗ FAIL: <30% occupied (misalignment)

3. **Cell Density Check**
   - ✓ PASS: 0.5-10 cells/spot (reasonable for 10µm)
   - ⚠ CAUTION: >10 cells/spot (possible over-mapping)
   - ✗ FAIL: 0 cells/spot (no mapping)

## Workflow Usage

### Option 1: Run Full Rmd Notebook in RStudio

```bash
# Set up environment
conda activate /projectnb/paxlab/presh/env/conda_env/giotto_env311
export LD_LIBRARY_PATH=/projectnb/paxlab/presh/env/conda_env/giotto_env311/lib:$LD_LIBRARY_PATH
export LD_PRELOAD=/projectnb/paxlab/presh/env/conda_env/giotto_env311/lib/libstdc++.so.6

# Start R
R --vanilla

# In R:
rmarkdown::render("/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/alignment/cell_to_spot_map.Rmd",
                  output_dir = "/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/alignment")
```

### Option 2: Create qsub Wrapper (Production)

```bash
# Create qsub script that runs Rmarkdown::render
# Submit to cluster for batch processing
qsub analysis/src/alignment/run_cell_to_spot_map.qsub.sh
```

### Option 3: Debug Individual Sections

```bash
# Run debug script standalone
Rscript analysis/src/alignment/debug_alignment.test.R

# Run as qsub job
qsub analysis/src/alignment/debug_alignment.qsub.sh
```

## Expected Output Files

### Diagnostic Plots
- `01_debug_overlay_all.png` - Full overlay comparison
- `02_debug_atac_spots.png` - ATAC spot locations
- `03_debug_xenium_cells.png` - Xenium cell locations
- `04_mapping_density_histogram.png` - Density distribution
- `05_mapping_spatial_distribution.png` - Spatial mapping status

### Data Files
- `xenium_to_atac_mapping.csv` - Detailed mapping (one row per Xenium cell)
- `atac_spot_density_stats.csv` - Density statistics (one row per ATAC spot)
- `mapping_summary_report.csv` - Summary metrics table

### Rmd Output
- `cell_to_spot_map.html` or `cell_to_spot_map.pdf` - Rendered report with all plots

## Key Improvements Over Original Workflow

1. **Comprehensive Debugging Section**
   - Validates alignment BEFORE mapping
   - Catches misalignment issues early
   - Provides diagnostic plots for visual assessment

2. **Better Error Handling**
   - Try/catch blocks around all critical operations
   - Graceful fallback to centroid-based mapping if boundaries unavailable
   - Detailed error messages for troubleshooting

3. **Enhanced Visualization**
   - All plots saved as high-resolution PNG (dpi=150)
   - Color coding for mapping status
   - Legends and descriptive titles
   - Summary statistics on plots

4. **Structured Output**
   - All results saved to CSV for downstream analysis
   - Consistent naming convention
   - Summary report for quick reference
   - Validation metrics logged to console

5. **Documentation**
   - Inline comments explaining each step
   - Progress indicators throughout
   - Clear section headers and subsections
   - Expected outcomes documented

## Troubleshooting Guide

### Issue: "No matching cell names between count_mat and df_allele"
- Cause: Different barcode sets used in different steps
- Solution: Review `xenium_to_atac_mapping.csv` - check if Xenium and ATAC cells align spatially

### Issue: Very few cells mapped (<50%)
- Cause: Possible misalignment or wrong cell boundary file
- Solution: Review `03_overlay_before_full.png` from MOSAICField - check if pre-alignment looked good

### Issue: Extremely high density (>50 cells/spot)
- Cause: Possible over-registration or wrong polygon size
- Solution: Reduce `dist` parameter from 5µm to 3µm to create smaller polygons

### Issue: Memory errors during object loading
- Cause: H5AD files are large (12GB ATAC)
- Solution: Ensure sufficient compute node memory (16GB+) or process in chunks

## Next Steps

1. **Review MOSAICField visualizations**
   - Check `01-04_*.png` in mosaicfield_outputs folder
   - Verify pre- and post-alignment look reasonable

2. **Run debug script**
   - Execute `debug_alignment.test.R` on compute node
   - Review diagnostic plots (01-03_debug_*.png)
   - Check console output for alignment status

3. **Run full cell_to_spot_map.Rmd**
   - Once debug looks good, run full notebook
   - Generate mapping plots and statistics

4. **Validate mapping results**
   - Review `04_mapping_density_histogram.png`
   - Check `mapping_summary_report.csv`
   - Confirm >70% ATAC spot occupancy

5. **Integration with downstream analysis**
   - Use `xenium_to_atac_mapping.csv` for joint analysis
   - Aggregate RNA counts by ATAC spot
   - Compare ATAC gene scores with Xenium RNA levels

## Reference

- MOSAICField GitHub: https://github.com/zhouyulab/MOSAICField
- Giotto Suite Documentation: https://ruv.github.io/GiottoClass/
- Xenium (10X) Documentation: https://www.10xgenomics.com/products/xenium

---

**Status**: Analysis notebooks and debug scripts complete and ready for execution.
**Last Updated**: 2026-06-08 15:13 UTC
**Maintained By**: Spatial ATAC Alignment Pipeline
