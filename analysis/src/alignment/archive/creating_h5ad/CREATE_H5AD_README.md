# Creating h5ad Objects from ArchR and Giotto

This guide explains how to create h5ad objects in the same format as the MOSAICField input data.

## ATAC from ArchR

### Step 1: Export from ArchR (R)

```R
# Edit these parameters in create_h5ad_from_ArchR.R:
archR_project_path <- "/path/to/your/ArchRProject"
output_dir <- "./atac_export"

# Run:
Rscript create_h5ad_from_ArchR.R
```

This exports:
- `atac_peak_matrix.mtx` — peak accessibility matrix (peaks × cells, Matrix Market format)
- `atac_peak_names.csv` — peak identifiers (chr:start-end format, or your naming)
- `atac_cell_names.csv` — cell barcodes
- `atac_coords.csv` — spatial coordinates (x, y or xcor, ycor)
- `atac_obs.csv` — per-cell metadata (QC metrics, cluster assignments, etc.)
- `atac_var.csv` — per-peak statistics (mean, variance, dispersions)
- `atac_peak_coords.csv` — genomic coordinates parsed from peak names

### Step 2: Build h5ad (Python)

```bash
python3 create_h5ad_from_export.py atac_export atac_from_archR.h5ad
```

Output: `atac_from_archR.h5ad` in MOSAICField-compatible format

**Resulting h5ad structure:**
```
AnnData object with n_obs × n_var = (cells, peaks)
    obs: cell_id, [QC columns from metadata]
    var: peak, mean, variance, dispersions, [genomic info]
    obsm:
        'spatial': cell spatial coordinates (n_obs, 2)
    uns:
        'peak_coords': genomic coordinate dict
```

## Xenium from Giotto

### Step 1: Export from Giotto (R)

```R
# Create similar export structure from Giotto object
library(Giotto)
library(data.table)

giotto_obj <- loadGiottoObj("path/to/giotto.RDS")

# Get expression matrix (cells × genes)
expr_mat <- giotto_obj@expression$rna$raw

# Get spatial coordinates
spat_coords <- giotto_obj@spatial_locs@coordinates

# Save as CSVs
fwrite(data.table(gene = rownames(expr_mat)), 
       "xenium_gene_names.csv")
fwrite(data.table(cell_id = colnames(expr_mat)), 
       "xenium_cell_names.csv")

# Sparse matrix
Matrix::writeMM(t(expr_mat), "xenium_expression.mtx")

# Coordinates (should have colnames like c("x", "y") or c("sdimx", "sdimy"))
spat_df <- as.data.table(spat_coords)
setnames(spat_df, old = names(spat_df), new = c("cell_id", "x", "y"))
fwrite(spat_df, "xenium_coords.csv")

# Gene stats
gene_means <- Matrix::colMeans(expr_mat)
gene_vars <- apply(expr_mat, 2, function(x) var(as.numeric(x)))
gene_var_df <- data.table(
  gene = colnames(expr_mat),
  mean = gene_means,
  variance = gene_vars,
  dispersions = gene_vars / (gene_means + 1)
)
fwrite(gene_var_df, "xenium_var.csv")
```

### Step 2: Build h5ad (Python — same approach)

```bash
python3 create_h5ad_from_export.py xenium_export xenium_from_giotto.h5ad
```

Adapt the column names in the Python script if your coordinate columns differ.

## h5ad Format Specification

Both ATAC and Xenium h5ad objects follow this structure:

| Component | Content | Size |
|-----------|---------|------|
| **X** | Gene expression / peak accessibility (sparse CSR) | obs × var |
| **obs** | Per-cell/spot metadata (barcodes, QC, clusters) | obs × metadata_cols |
| **var** | Per-gene/peak statistics (mean, variance, dispersions) | var × stat_cols |
| **obsm['spatial']** | Spatial coordinates in micrometers (µm) | obs × 2 |
| **uns** | Unstructured metadata (optional: peak genomics, alignment params) | dict |

### Example: MOSAICField ATAC Format

```python
import anndata as ad

adata = ad.read_h5ad("atac_from_archR.h5ad")

print(adata)
# Output:
# AnnData object with n_obs × n_var = (11640, 2000)
#     obs: 'cell_id', 'n_fragment', 'frac_dup', 'frac_mito', 'cluster', ...
#     var: 'peak', 'mean', 'variance', 'dispersions', 'peak_idx', ...
#     obsm: 'spatial'
#     uns: 'peak_coords'

# Access data
coords = adata.obsm['spatial']  # (11640, 2) in micrometers
peaks = adata.var_names          # 2000 peak IDs
cells = adata.obs_names           # 11640 cell barcodes
```

## Tips

1. **Peak naming**: MOSAICField expects genomic bins in format `chr:start-end` (e.g., `chr1:10000-15000`). Adjust parsing in the Python script if your naming differs.

2. **Spatial coordinates**: Must be in **micrometers (µm)**, not pixels. ArchR/Giotto typically store original slide coordinates.

3. **Sparsity**: Keep matrices as sparse (CSR) for memory efficiency. The Python script preserves this.

4. **QC columns**: Include relevant QC metrics in `obs`:
   - ATAC: `n_fragment`, `TSSEnrichment`, `ReadsInPeaks`, `nMito`, `cluster`
   - Xenium: gene counts, UMI counts, etc.

5. **Testing**: After creating h5ad, load and inspect:
   ```python
   import anndata as ad
   adata = ad.read_h5ad("your_output.h5ad")
   print(adata)
   print(adata.obsm['spatial'].shape)
   print(adata.X.shape, adata.X.nnz / adata.X.shape[0] / adata.X.shape[1])
   ```

## Running MOSAICField

Once you have both `atac_from_archR.h5ad` and `xenium_from_giotto.h5ad`:

```python
import anndata as ad
slice_source = ad.read_h5ad("atac_from_archR.h5ad")
slice_target = ad.read_h5ad("xenium_from_giotto.h5ad")

# Feed into MOSAICField pipeline (see mosaic_run2-2.ipynb STEP 0)
```
