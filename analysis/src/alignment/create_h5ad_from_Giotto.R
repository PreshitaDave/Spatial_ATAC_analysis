# ==============================================================================
# Export Xenium data from Giotto object as CSVs for h5ad creation in Python
# ==============================================================================

library(Giotto)
library(data.table)
library(Matrix)

# --- Parameters ---
giotto_obj_path <- "path/to/your/giotto_object.RDS"  # CHANGE THIS
output_dir <- "./xenium_export"
dir.create(output_dir, showWarnings = FALSE)

set.seed(42)

# Load Giotto object
cat("Loading Giotto object...\n")
giotto_obj <- readRDS(giotto_obj_path)
cat(sprintf("Loaded: %d cells\n", ncol(giotto_obj)))

# ==============================================================================
# 1. Get gene expression matrix (cells × genes)
# ==============================================================================

cat("\nExtracting expression matrix...\n")

# Get raw expression (adjust slot if using different normalization)
# Common slots: 'raw', 'normalized', 'scaled'
expr_mat <- giotto_obj@expression$rna$raw

# If that doesn't work, try:
# expr_mat <- exprs(giotto_obj, assay_name = 'rna', values = 'raw')

if (is.null(expr_mat)) {
  stop("Could not extract expression matrix. Check object structure.")
}

cat(sprintf("Expression matrix: %d genes × %d cells\n", nrow(expr_mat), ncol(expr_mat)))

# Convert to sparse if dense
if (!inherits(expr_mat, "sparseMatrix")) {
  expr_mat <- Matrix::Matrix(expr_mat, sparse = TRUE)
}

# Get gene and cell IDs
gene_ids <- rownames(expr_mat)
cell_ids <- colnames(expr_mat)

cat(sprintf("Example genes: %s, %s, ...\n", gene_ids[1], gene_ids[2]))
cat(sprintf("Example cells: %s, %s, ...\n", cell_ids[1], cell_ids[2]))

# ==============================================================================
# 2. Get spatial coordinates
# ==============================================================================

cat("\nExtracting spatial coordinates...\n")

# Giotto stores spatial coords in spatialCoords slot
spat_locs <- giotto_obj@spatial_locs

# Extract coordinates (typically 2D: x, y)
coords_matrix <- spat_locs@coordinates

if (is.null(coords_matrix)) {
  stop("Could not find spatial coordinates. Check spatial_locs slot.")
}

cat(sprintf("Spatial coords: %d cells\n", nrow(coords_matrix)))
cat(sprintf("  Dims: x=[%.1f, %.1f], y=[%.1f, %.1f]\n",
            min(coords_matrix[,1]), max(coords_matrix[,1]),
            min(coords_matrix[,2]), max(coords_matrix[,2])))

# Convert to data.table with proper names
if (ncol(coords_matrix) >= 2) {
  coords_df <- data.table(
    cell_id = rownames(coords_matrix),
    x = coords_matrix[,1],
    y = coords_matrix[,2]
  )
  if (ncol(coords_matrix) > 2) {
    coords_df[, z := coords_matrix[,3]]
  }
} else {
  stop("Spatial coordinates should have at least 2 dimensions")
}

# ==============================================================================
# 3. Get cell metadata (optional but recommended)
# ==============================================================================

cat("\nExtracting cell metadata...\n")

# Giotto stores metadata in pDataDT
if (!is.null(giotto_obj@cell_metadata)) {
  obs_df <- as.data.table(giotto_obj@cell_metadata)
  obs_df[, cell_id := rownames(giotto_obj@cell_metadata)]
} else {
  obs_df <- data.table(cell_id = cell_ids)
}

setcolorder(obs_df, "cell_id")
cat(sprintf("Obs columns: %s\n", paste(colnames(obs_df), collapse = ", ")))

# ==============================================================================
# 4. Compute gene statistics
# ==============================================================================

cat("\nComputing gene statistics...\n")

gene_means <- Matrix::colMeans(expr_mat)
gene_vars <- apply(expr_mat, 2, function(x) var(as.numeric(x)))

var_df <- data.table(
  gene = gene_ids,
  mean = gene_means,
  variance = gene_vars,
  dispersions = gene_vars / (gene_means + 1),
  gene_idx = seq_along(gene_ids)
)

cat(sprintf("Gene stats: %d genes\n", nrow(var_df)))
cat(sprintf("  Mean range: [%.3f, %.3f]\n", min(var_df$mean), max(var_df$mean)))
cat(sprintf("  Variance range: [%.3f, %.3f]\n", min(var_df$variance), max(var_df$variance)))

# ==============================================================================
# 5. Export files
# ==============================================================================

cat("\nExporting files...\n")

# Expression matrix (cells × genes in MTX format)
# Note: mmwrite expects (rows × cols), so we transpose to (genes × cells)
writeMM(t(expr_mat), file = file.path(output_dir, "xenium_expression.mtx"))
cat(sprintf("Wrote xenium_expression.mtx (%d genes × %d cells)\n",
            nrow(expr_mat), ncol(expr_mat)))

# Gene and cell names
fwrite(data.table(gene = gene_ids),
       file.path(output_dir, "xenium_gene_names.csv"))
fwrite(data.table(cell_id = cell_ids),
       file.path(output_dir, "xenium_cell_names.csv"))

# Coordinates
fwrite(coords_df, file.path(output_dir, "xenium_coords.csv"))

# Obs metadata
fwrite(obs_df, file.path(output_dir, "xenium_obs.csv"))

# Var (gene statistics)
fwrite(var_df, file.path(output_dir, "xenium_var.csv"))

# ==============================================================================
# 6. Summary
# ==============================================================================

cat(sprintf("\n%s\n", strrep("=", 80)))
cat(sprintf("Export complete: %s\n", output_dir))
cat(sprintf("%s\n", strrep("=", 80)))
cat("Files:\n")
for (f in list.files(output_dir)) {
  fsize <- file.size(file.path(output_dir, f)) / 1024
  cat(sprintf("  %s (%.0f KB)\n", f, fsize))
}
cat("\nNext: run create_h5ad_from_export.py to build h5ad\n")
cat("Example:\n")
cat(sprintf("  python3 create_h5ad_from_export.py %s xenium_from_giotto.h5ad\n",
            output_dir))
