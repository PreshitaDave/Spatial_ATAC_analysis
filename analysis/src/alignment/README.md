# MOSAICField Alignment Analysis — Workflow & Results

**Date**: June 2026  
**Tissue**: Breast cancer spatial ATAC-seq (Xenium + ArchR)  
**Goal**: Validate coordinate alignment quality between Xenium cells and ATAC spots

---

## Quick Summary

MOSAICField was used to spatially align Xenium (single-cell gene expression) coordinates to the ATAC spot grid. Three types of analyses were performed:

1. **Alignment quality**: affine vs nonlinear warping
2. **Cell-to-spot assignment**: nearest-neighbor vs Voronoi
3. **Gene-level validation**: Pearson correlation between Xenium RNA and ArchR GeneScore

**Bottom line**: Nonlinear alignment + Voronoi assignment is recommended; correlations are low (~0.017–0.095), suggesting good coarse-scale accuracy but fine-scale uncertainty.

---

## How to Read This Analysis

### Phase 1: Alignment visualization (`MOSAICField.ipynb`)
1. **Affine baseline**: `step0_affine_aligned.png` — before/after affine warp
2. **Input data**: `step1_crop.png`, `step1_pc_channels.png`, `step1_rasterized_shared.png`
3. **Nonlinear warp**: `step2_deformation_field.png`, `step2_nonlinear_warp.png`
4. **Warped coordinates**: `step3_before_after_arrows.png`, `step3_warped_coords.png`
5. **Cell assignment**: `step4_cutoff_comparison.png`, `step4b_method_comparison.png`, Voronoi zoom views
6. **Spot coverage**: `step5_*` (assignment overview + zoom), `step6_*` (cells per spot)
7. **Biology**: `step7_cell_types.png`, `step7_tumor_purity.png`, `step7_purity_spatial_bins.png`

### Phase 2: Gene-level validation (`gene_loss_evaluation.ipynb`)
8. **Start here**: `GENE_LOSS_EVALUATION_SUMMARY.txt` — text summary of correlations
9. **Per-gene Pearson**: `gene_loss_plot1_per_gene_pearson_boxplot_v2.png` — compare all 4 methods
10. **Resolution sweep**: `gene_loss_plot2_resolution_sweep_v2.png` — correlations improve at coarser scales
11. **Coverage vs precision**: `gene_loss_plot4_coverage_precision_v2.png` — NN vs Voronoi trade-off
12. **Per-spot agreement**: `gene_loss_plot3_per_spot_pearson_boxplot_v2.png`
13. **Method comparison**: `gene_loss_plot5_methodB_vs_A_paired_diff_v2.png`

### Phase 3: Interpret results
14. **Read**: Caveats section below (validation limitations)

---

## Current Analysis Plots

### Alignment Visualization (from `MOSAICField.ipynb`, all in `mosaicfield_outputs/`)

| Plot | What it shows | File |
|------|---------------|------|
| **Step 0** | ATAC spots before/after affine transformation | `step0_affine_aligned.png` |
| **Step 1: Crop** | Zoomed view of ATAC + Xenium before alignment | `step1_crop.png` |
| **Step 1: Channels** | PCA channels used for registration | `step1_pc_channels.png` |
| **Step 1: Rasterized** | ATAC accessibility heatmap (chromatin signal) | `step1_rasterized_shared.png` |
| **Step 2: Deformation** | Nonlinear warp vectors | `step2_deformation_field.png` |
| **Step 2: Warp** | ATAC spots after nonlinear alignment | `step2_nonlinear_warp.png` |
| **Step 2b: QC** | Jacobian severity map, distribution histogram, vortex hotspots, flow coherence, pass/fail verdict | `step2b_warp_quality_control.png` |
| **Step 3: Arrows** | Displacement vectors for nonlinear warping | `step3_before_after_arrows.png` |
| **Step 3: Coords** | Scatter: original vs warped positions | `step3_warped_coords.png` |
| **Step 4: Cutoff** | NN assignment vs soft Voronoi cutoff | `step4_cutoff_comparison.png` |
| **Step 4b: Comparison** | NN vs Voronoi head-to-head overlay | `step4b_method_comparison.png` |
| **Step 4b: Voronoi** | Voronoi cell assignments across full tissue | `step4b_voronoi_overview.png` |
| **Step 4b: Zoom wide** | Voronoi zoomed to ~1000 µm | `step4b_voronoi_zoom_wide.png` |
| **Step 4b: Zoom close** | Voronoi zoomed to ~200 µm | `step4b_voronoi_zoom_close.png` |
| **Step 5: Overview** | All Xenium cell assignments per spot | `step5_assignment_overview.png` |
| **Step 5: Single spot** | Example: 1 spot → multiple cell assignments | `step5_single_spot_zoom.png` |
| **Step 5: Zoom wide** | Assignment pattern ~1000 µm | `step5_zoom_wide.png` |
| **Step 5: Zoom close** | Assignment pattern ~200 µm | `step5_zoom_close.png` |
| **Step 6: Distribution** | Histogram of cell count per spot | `step6_cells_per_spot_dist.png` |
| **Step 6: Spatial** | Spatial map of cell density per spot | `step6_cells_per_spot_spatial.png` |
| **Step 6: Density bins** | Binned density map (25 µm bins) | `step6_spot_density_bins.png` |
| **Step 7: Cell types** | Xenium cell type composition per spot | `step7_cell_types.png` |
| **Step 7: Purity** | Fraction tumor cells per spot (spatial) | `step7_tumor_purity.png` |
| **Step 7: Purity bins** | Tumor purity in spatial bins | `step7_purity_spatial_bins.png` |

### Gene-Level Validation (from `gene_loss_evaluation.ipynb`)

| Plot | What it shows | File |
|------|---------------|------|
| **Gene loss plot 1** | Per-gene Pearson correlation boxplot (all 4 conditions) | `gene_loss_plot1_per_gene_pearson_boxplot_v2.png` |
| **Gene loss plot 2** | Median Pearson vs spatial resolution (0–400 µm) | `gene_loss_plot2_resolution_sweep_v2.png` |
| **Gene loss plot 3** | Per-spot Pearson correlation boxplot | `gene_loss_plot3_per_spot_pearson_boxplot_v2.png` |
| **Gene loss plot 4** | Cells kept vs correlation (coverage/precision trade-off) | `gene_loss_plot4_coverage_precision_v2.png` |
| **Gene loss plot 5** | Method A vs B paired difference (nonlinear fixed) | `gene_loss_plot5_methodB_vs_A_paired_diff_v2.png` |

---

## Results Summary

### Warp Quality Control (June 29, 2026)

**File**: `mosaicfield_outputs/step2b_warp_quality_control.png`  
**Notebook cell**: STEP 2b in `mosaic_run2-2.ipynb`  
**Method**: Displacement derived from `phi - identity_grid` (absolute coord field minus pixel positions); Jacobian computed via finite differences on forward displacement.

| Check | Result | Notes |
|---|---|---|
| Fold-free (det≥0) | FAIL* | 1,337 px (1.76%) with det<0 — all at canvas edges (boundary artifact) |
| Topology (<1% bad pixels) | WARN* | 1.82% total bad — same boundary artifact, not interior tissue |
| Flow (edge-bias >2) | **PASS** | 3.04× — large displacements concentrated at slide edges, consistent with global correction |
| Vortex (≤5 components) | WARN | 11 components, largest=67 px — all tiny and at canvas boundary, not interior singularities |

**\*False alarms**: The FAIL/WARN flags arise from the finite-difference Jacobian computation degrading at the canvas border (no pixels outside). The tissue interior is 97.8% acceptable (det≥0.5). No true tissue folding or interior vortices were detected. Warp is physically sound and safe to proceed.

**Displacement statistics**: mean=225 µm, max=545 µm. Large values are at the slide edges and represent real global rotation/stretch correction.

---

### Gene Loss Evaluation (June 26, 2026)

**File**: `GENE_LOSS_EVALUATION_SUMMARY.txt`

#### Correlations at Native Spot Resolution

| Condition | Median per-gene Pearson | Global Pearson | Cells kept |
|-----------|-------------------------|-----------------|-----------|
| Affine + NN (Method A) | 0.011 | 0.129 | 91.3% |
| Affine + Voronoi (Method B) | 0.011 | 0.129 | 99.7% |
| Nonlinear + NN (Method A) | 0.016 | 0.138 | 91.3% |
| **Nonlinear + Voronoi (Method B)** | **0.017** | **0.139** | **99.7%** |

**Best condition**: Nonlinear + Method B (Voronoi)

#### Correlation Improves at Coarser Spatial Scales

Nonlinear + Voronoi, median per-gene Pearson:
- **Native (0 µm)**: 0.017
- **25 µm bins**: 0.016
- **50 µm bins**: 0.022
- **100 µm bins**: 0.058
- **200 µm bins**: 0.082
- **400 µm bins**: 0.095

**Interpretation**: Agreement improves substantially at coarser scales, consistent with good alignment at tissue-region level but uncertainty at single-spot resolution.

#### Method A vs Method B (Nonlinear Fixed)

Voronoi (Method B) marginally better but retains more cells:
- **Cells kept**: 99.7% vs 91.3%
- **Mean assignment distance**: 13.0 µm vs 12.1 µm
- **Per-gene Pearson gain**: +0.0003 (56% of genes better, Wilcoxon p=0.049)
- **Per-spot Pearson gain**: +0.001 (18.5% better, p<0.0001)

---

## Important Caveats & Interpretation

### 1. Correlations Are Very Low
Median per-gene Pearson of **0.017 at native resolution** is close to noise. Even at 200 µm bins, it only reaches **~0.08–0.09**. These values do not strongly validate the alignment.

### 2. GeneScoreMatrix ≠ RNA Expression
The ArchR gene activity score is **model-inferred** from chromatin accessibility, not measured RNA. It's imputed twice:
- Once: ATAC → gene activity (ArchR's imputation)
- Twice: Xenium RNA assigned to ATAC spots (our spatial mapping)

Low correlations could reflect either bad alignment OR lossy imputation at each step.

### 3. What the Resolution Sweep Tells You
Correlation improving from **0.017 → 0.095** (native → 200 µm) suggests:
- ✓ Alignment is approximately correct at **tissue-region scale** (coarse)
- ✗ Cell-level assignment is **noisy** at single-spot resolution

The biology is real in aggregate; single-spot assignments are uncertain.

### 4. Nonlinear Advantage
Nonlinear alignment wins by ~**50% over affine** (0.017 vs 0.011), a meaningful relative improvement. This suggests the nonlinear warp captures real tissue deformation, but does **NOT** prove cell-level accuracy.

---

## Recommended Next Steps for Validation

### Option 1: scRNA-seq Label Transfer (Recommended First)
- Cluster scRNA-seq data, assign cell-type labels from RNA markers
- Project scRNA cell types onto spatial ATAC spots via nonlinear + Method B mapping
- Check if ArchR ATAC clusters correspond to spatially coherent scRNA cell types
- **Advantage**: Fully independent of GeneScoreMatrix; uses measured RNA

### Option 2: Direct Peak-to-Gene Correlation
- Use scRNA-seq gene expression as ground truth (instead of Xenium pseudobulk)
- Drop-in replacement for Xenium in the loss function
- Same notebook code is reusable

### Option 3: Co-Embedding / Integration
- WNN or Seurat v4 bridge integration of ATAC + RNA
- Joint embedding quality measures modality agreement
- Strongest validation but most effort

---

## Data Overview

### Input Files
- **ATAC**: ArchR project (11,319 cells after QC in v2)
- **Xenium**: Gene expression matrix (377 shared genes with ArchR)
- **Alignment**: MOSAICField nonlinear registration

### Generated Outputs

**Alignment & Mapping**:
- `atac_affine_aligned.h5ad` — ATAC spots after affine warp
- `atac_nonlinear_aligned.h5ad` — ATAC spots after nonlinear warp
- `xenium_to_atac_mapping.csv` — Nearest-neighbor assignments
- `xenium_to_atac_mapping_voronoi.csv` — Voronoi assignments
- `atac_xenium_mapping.csv` — Forward map (spot → cells)

**Gene Loss Analysis**:
- `gene_loss_inputs_v1/` — Original ArchR project (7,842 cells)
- `gene_loss_inputs_v2/` — Corrected ArchR project (11,319 cells)
- `per_spot_tumor_purity.csv` — Cell-type composition per spot

### What Changed: v1 vs v2

| Aspect | v1 | v2 | Reason |
|--------|----|----|--------|
| **ArchR cells** | 7,842 | 11,319 | Project re-built with corrected filtering |
| **Gene expression matches** | Lower | Higher | More cells → more reliable pseudobulk |
| **Correlation stability** | Noisier | Better | More cells reduce per-spot noise |
| **Used in notebook** | Early runs | Final analysis | Improved data quality |

The v2 inputs use a corrected ArchR project with more cells after proper QC filtering.

---

## How to Run the Full Analysis

### Prerequisites
```bash
# Install dependencies (if not already in environment)
pip install anndata scanpy mosaicfield
```

### Run Gene Loss Evaluation Notebook
```bash
cd /projectnb/paxlab/presh/projects/spatial_atac/analysis/src/alignment
jupyter notebook gene_loss_evaluation.ipynb
```

Key cells to run in order:
1. **Load libraries & setup** → define `gli_dir = gene_loss_inputs_v2` (v2 data)
2. **Load aligned h5ad files** → ATAC and Xenium coordinates
3. **Build pseudobulk** → Xenium RNA aggregated per spot
4. **Compute correlations** → per-gene and per-spot Pearson
5. **Resolution sweep** → binned correlations at 0/25/50/100/200/400 µm
6. **Plots** → all gene loss plots are generated

### Export Aligned Data
```bash
Rscript export_genescore_for_loss.R
```
This prepares the ArchR GeneScore matrix for correlation analysis.

---

## File Structure

```
analysis/src/alignment/
├── README.md                                    ← YOU ARE HERE
├── gene_loss_evaluation.ipynb                   ← Main analysis notebook
├── GENE_LOSS_EVALUATION_SUMMARY.txt             ← Results summary
├── MOSAICField/                                 ← Alignment software repo
├── export_genescore_for_loss.R                  ← Prepare ATAC data
├── create_h5ad_from_*.R / create_h5ad_from_*.py ← Data format conversion
└── mosaicfield_outputs/
    ├── atac_affine_aligned.h5ad                 ← Aligned ATAC (affine)
    ├── atac_nonlinear_aligned.h5ad              ← Aligned ATAC (nonlinear)
    ├── xenium_affine_aligned.h5ad               ← Aligned Xenium (affine)
    ├── xenium_nonlinear_aligned.h5ad            ← Aligned Xenium (nonlinear)
    ├── gene_loss_plot*_v2.png                   ← Gene loss evaluation plots (current)
    ├── gene_loss_inputs/                        ← v1: original ArchR cells (7,842)
    ├── gene_loss_inputs_v2/                     ← v2: corrected ArchR cells (11,319)
    ├── atac_xenium_mapping.csv                  ← Forward mapping (spot → cells)
    └── per_spot_tumor_purity.csv                ← Cell-type composition
```

---

## Questions / Troubleshooting

**Q: Why are correlations so low?**
A: See "Caveats" section above. GeneScore is imputed from ATAC; Xenium RNA is assigned to spots via spatial mapping. Both introduce noise. Coarse-scale correlations (~200 µm) are better (~0.08–0.09).

**Q: Which alignment should I use for downstream analysis?**
A: **Nonlinear + Voronoi (Method B)**. It has the best gene correlation (0.017 vs 0.011–0.016) and retains 99.7% of cells.

**Q: What if I have scRNA-seq data?**
A: Use scRNA-seq label transfer (Option 1 above) to validate independently. Don't rely solely on GeneScore correlations.

**Q: Can I trust single-spot assignments?**
A: Be cautious. Resolution sweep shows uncertainty increases at fine scale. Recommend binning to ≥100 µm or using aggregated cluster-level assignments.

---

**Last updated**: 2026-06-29  
**Notebook**: `gene_loss_evaluation.ipynb`  
**Contact**: preshita@bu.edu
