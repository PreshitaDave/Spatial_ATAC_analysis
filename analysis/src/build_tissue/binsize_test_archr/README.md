# Binsize Test ArchR

Compare tile sizes (500bp vs 5000bp) and binarization options for ArchR TileMatrix performance, with downstream cluster annotation and normal barcode extraction for CalicoST.

## Script Numbering (Execution Order)

| Script | Purpose |
|--------|---------|
| `1_create_arrow_variants.R` | Create arrow files for all tile/binarize combos (permissive QC, applies no-edge-effect whitelist) |
| `2_compare_arrow_sparsity.R` | Aggregate sparsity metrics → CSV + plots + recommendations |
| `3_build_archr_variant_project.R` | Build ArchR project per variant (LSI, clusters res=0.8, UMAP, gene scores) |
| `4_compare_archr_variant_umap.R` | Aggregate per-variant UMAPs into comparison PDF |
| `5_compare_spatial_smoothing_methods.R` | Build 3 smoothing variants (none / alpha-blend / iterative) per tissue |
| `6_add_spatial_smoothing_lsi.R` | Apply spatial kNN smoothing to LSI embedding |
| `7_annotate_clusters_by_markers.R` | Score clusters against 16-identity reference panel (fibroblast, T/B cell, macrophage, etc.) |
| `8_model_spot_count_distributions.R` | Per-tile count distributions (NB fitting, mean-variance diagnostics) |
| `9_extract_normal_barcodes_calicost.R` | Identify normal (immune/stromal) clusters via module scores; export barcodes for CalicoST `normalidx_file` |

### Deprecated scripts (kept for reference, do not use)
- `create_arrow_variants_v1_deprecated.R` — stricter QC, barcode whitelist not applied
- `build_archr_variant_project_v1_deprecated.R` — inherits arrow QC without post-load correction

## Quick Start

### Metrics recomputation (arrow files already exist)
```bash
cd /projectnb/paxlab/presh/projects/spatial_atac
bash analysis/src/build_tissue/binsize_test_archr/submit_metrics_recompute.sh
```
Submits 16 parallel metrics jobs + 1 aggregation job.

### UMAP/cluster validation (optional pilot)
```bash
bash analysis/src/build_tissue/binsize_test_archr/submit_archr_variant_pilot.sh
```

### Spatial smoothing (all 4 tissues)
```bash
bash analysis/src/build_tissue/binsize_test_archr/submit_spatial_smoothing.sh
```

### Extract normal barcodes for CalicoST (per tissue)
```bash
module load R
Rscript analysis/src/build_tissue/binsize_test_archr/9_extract_normal_barcodes_calicost.R lowseq_489
```
Output: `analysis/binsize_comparison/normal_barcodes/lowseq_489_normal_barcodes.csv`

## Chosen Parameters

For **lowseq_489**: `5000bp, binarize=FALSE` (selected based on sparsity analysis and CNV pipeline requirements).

## Output Locations

| Output | Path |
|--------|------|
| Sparsity metrics | `analysis/binsize_comparison/sparsity_comparison_table.csv` |
| Sparsity plots | `analysis/binsize_comparison/sparsity_comparison_plots.pdf` |
| UMAP CSVs | `analysis/binsize_comparison/archr_projects/<tissue>_<size>_<binarize>_umap.csv` |
| Smoothing PDFs | `analysis/binsize_comparison/spatial_smoothing/` |
| Cluster annotations | `analysis/binsize_comparison/cluster_annotations/` |
| Normal barcodes | `analysis/binsize_comparison/normal_barcodes/<tissue>_normal_barcodes.csv` |

See **WORKFLOW.md** for detailed technical documentation.
