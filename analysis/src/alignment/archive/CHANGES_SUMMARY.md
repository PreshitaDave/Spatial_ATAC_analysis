# Summary of Changes to cell_to_spot_map.Rmd

## What's Been Added/Updated

### 1. **New Section 1: Environment Initialization**
- Proper library loading with error handling
- Python environment configuration
- Giotto backend initialization
- Path setup for inputs/outputs

### 2. **New Section 2: Debugging - Alignment Quality Validation**
**Subsections:**
- Load aligned ATAC and Xenium objects from MOSAICField outputs
- Validate spatial coordinates (x, y ranges)
- Assess alignment quality by comparing bounding box areas
- Generate 3 diagnostic overlay plots:
  - `01_debug_overlay_all.png` - Combined view
  - `02_debug_atac_spots.png` - ATAC only
  - `03_debug_xenium_cells.png` - Xenium only

**Purpose:** Catch alignment issues BEFORE proceeding to mapping

### 3. **Reorganized Section 3: ATAC-to-Xenium Mapping**
**Now includes 7 ordered subsections:**
1. Create 10µm square polygons from ATAC centroids
2. Load Xenium cell boundaries from parquet file
3. Add ATAC polygons to Xenium Giotto object
4. Calculate overlap using Giotto's `calculateOverlap()` function
5. Extract mapping results with statistics
6. Generate publication-quality visualization plots
7. Save all results to CSV files

**Key Improvements:**
- Explicit error handling with tryCatch
- Fallback to centroid-based mapping if boundaries unavailable
- Detailed console logging with timestamps
- Automatic plot generation and saving

### 4. **New Section 4: Validation Checks**
Automated QC metrics:
- Check 1: Mapping coverage (% cells mapped to ATAC)
- Check 2: Spot occupancy (% ATAC spots with cells)
- Check 3: Cell density reasonableness

Each check has PASS/CAUTION/FAIL verdicts

### 5. **Session Info Section**
Captures R/package versions for reproducibility

## Original Code Preserved

✅ All original mapping logic is intact:
- ATAC polygon creation (10µm squares)
- Xenium boundary loading
- Giotto polygon operations
- Overlap calculation
- Cell density statistics
- CSV export

The original code has been:
1. **Reorganized** into clear sections
2. **Enhanced** with error handling
3. **Extended** with debugging steps
4. **Documented** with progress indicators
5. **Visualized** with publication plots

## New Output Files Generated

### Diagnostic Plots
- `01_debug_overlay_all.png`
- `02_debug_atac_spots.png`
- `03_debug_xenium_cells.png`
- `04_mapping_density_histogram.png`
- `05_mapping_spatial_distribution.png`

### Data Tables
- `xenium_to_atac_mapping.csv`
- `atac_spot_density_stats.csv`
- `mapping_summary_report.csv`

## How to Run

### In RStudio (Interactive)
```R
rmarkdown::render(
  "/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/alignment/cell_to_spot_map.Rmd",
  output_dir = "/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/alignment"
)
```

### Via qsub (Batch)
Create `run_cell_to_spot_map.qsub.sh`:
```bash
#!/bin/bash
#$ -l h_rt=04:00:00
#$ -pe omp 4
#$ -P paxlab
#$ -l mem_per_core=8G

module load R

cd /projectnb/paxlab/presh/projects/spatial_atac

Rscript -e "rmarkdown::render('/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/alignment/cell_to_spot_map.Rmd')"
```

## Key Parameters

- **Polygon size**: 10µm squares (dist = 5µm half-width)
- **Cell boundaries file**: `/projectnb/paxlab/DATA/DriesSpatial/Xenium/output-XETG00253__0048833__488B__20241217__182301/cell_boundaries.parquet`
- **Input objects**: MOSAICField outputs (H5AD format)
- **Output directory**: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/alignment/`

## Expected Performance Metrics

**For well-aligned datasets:**
- ✓ Mapping coverage: >90% of Xenium cells
- ✓ Spot occupancy: >70% of ATAC spots have ≥1 cell
- ✓ Mean density: 0.5-10 cells/spot
- ✓ Area ratio (xenium/atac): 0.5-2.0

**Indicators of alignment issues:**
- ✗ <50% mapping coverage
- ✗ <30% spot occupancy
- ✗ Area ratio outside 0.1-10 range

## Debugging Workflow

1. Run `debug_alignment.test.R` to validate alignment
   ```bash
   qsub analysis/src/alignment/debug_alignment.qsub.sh
   ```

2. Review diagnostic plots:
   - Are ATAC and Xenium in similar spatial regions?
   - Is coordinate scale similar?
   - Do the overlay plots look reasonable?

3. If alignment looks good:
   - Proceed with full `cell_to_spot_map.Rmd`
   - Review mapping statistics
   - Check generated CSV files

4. If alignment looks poor:
   - Review MOSAICField parameters
   - Check input data quality
   - Re-run MOSAICField alignment

## Status

✅ **cell_to_spot_map.Rmd** - Updated and ready to use
✅ **debug_alignment.test.R** - Created for QC validation
✅ **debug_alignment.qsub.sh** - Created for cluster submission
✅ **ALIGNMENT_ANALYSIS_GUIDE.md** - Complete reference document

**Next Step**: Run debug job and validate alignment quality

---
Date: 2026-06-08
Job ID: 5985840 (debug_alignment.qsub.sh - running)
