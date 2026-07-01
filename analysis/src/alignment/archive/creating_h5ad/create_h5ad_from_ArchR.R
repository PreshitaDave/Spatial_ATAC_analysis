# ==============================================================================
# Export ATAC data from ArchR project as CSVs for h5ad creation in Python
# ==============================================================================

library(ArchR)
library(data.table)
library(rhdf5)


# --- Parameters ---
archR_project_path <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/01_outputs/archR_objects/deepseq_488B/deepseq_488B_archR_project_final"  # CHANGE THIS
output_dir <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/alignment/atac_export"# CHANGE THIS
dir.create(output_dir, showWarnings = FALSE)

set.seed(42)

# Load ArchR project
proj <- loadArchRProject(archR_project_path)
cat(sprintf("Loaded ArchR project: %d cells\n", ncol(proj)))

# ==============================================================================
# 1. Get Tile matrix (cells × tiles)
# ==============================================================================

cat("\nExtracting Tile matrix...\n")

getMatrixFromArrow(
  ArrowFile = '/projectnb/paxlab/presh/projects/spatial_atac/Data/01_outputs/archR_objects/deepseq_488B/deepseq_488B_archR_project_final/ArrowFiles',
  useMatrix = "TileMatrix",
  useSeqnames = NULL,
  excludeChr = NULL,
  cellNames = NULL,
  ArchRProj = NULL,
  verbose = TRUE,
  binarize = FALSE,
  logFile = createLogFile("getMatrixFromArrow")
)

# Get tile matrix 
tile_mat <- getMatrixFromProject(proj, useMatrix = "TileMatrix", binarize = T)
cat(sprintf("tile matrix: %d cells × %d tiles\n", nrow(tile_mat), ncol(tile_mat)))

# Convert to sparse matrix if not already
if (!class(tile_mat)[1] %in% c("dgCMatrix", "dgTMatrix")) {
  tile_mat <- Matrix::Matrix(tile_mat, sparse = TRUE)
}

# Get tile names (seqnames:start-end)
tile_names <- rownames(tile_mat)
cat(sprintf("tile names format: %s\n", tile_names[1]))

# Get cell IDs
cell_ids <- colnames(tile_mat)
cat(sprintf("Example cell IDs: %s, %s\n", cell_ids[1], cell_ids[2]))

# ==============================================================================
# 2. Get spatial coordinates
# ==============================================================================

cat("\nExtracting spatial coordinates...\n")

# Get metadata with spatial coordinates
metadata <- getCellColData(proj)
head(metadata)

# ArchR uses '#ATAC#' prefix for barcode, with spatial coords in metadata
# Adjust these column names based on your metadata
if ("x" %in% colnames(metadata) && "y" %in% colnames(metadata)) {
  coords_df <- data.table(
    cell_id = rownames(metadata),
    x = metadata$x,
    y = metadata$y
  )
} else if ("xcor" %in% colnames(metadata) && "ycor" %in% colnames(metadata)) {
  coords_df <- data.table(
    cell_id = rownames(metadata),
    x = metadata$xcor,
    y = metadata$ycor
  )
} else {
  stop("Could not find spatial coordinates. Check metadata column names.")
}

cat(sprintf("Spatial coords: %d spots, x=[%.1f, %.1f], y=[%.1f, %.1f]\n",
            nrow(coords_df), min(coords_df$x), max(coords_df$x),
            min(coords_df$y), max(coords_df$y)))

# ==============================================================================
# 3. Get obs (per-spot metadata)
# ==============================================================================

cat("\nExtracting per-spot metadata...\n")

obs_df <- as.data.table(metadata)
obs_df[, cell_id := rownames(metadata)]
setcolorder(obs_df, "cell_id")

# Select QC columns (adjust to your metadata)
qc_cols <- intersect(colnames(obs_df),
                     c("nFrags", "TSSEnrichment", "ReadsIntiles",
                       "nMito", "log10_frags", "tsse", "cluster",
                       "sample", "condition"))
if (length(qc_cols) > 0) {
  obs_df <- obs_df[, c("cell_id", qc_cols), with = FALSE]
} else {
  obs_df <- obs_df[, .(cell_id)]
}

cat(sprintf("Obs columns: %s\n", paste(colnames(obs_df), collapse = ", ")))

# ==============================================================================
# 4. Get var (per-tile statistics)
# ==============================================================================

cat("\nExtracting per-tile statistics...\n")

# tile statistics: mean, variance, etc.
tile_means <- Matrix::rowMeans(tile_mat)
tile_vars <- Matrix::rowSums((tile_mat - tile_means)^2) / (ncol(tile_mat) - 1)

var_df <- data.table(
  tile = tile_names,
  mean = tile_means,
  variance = tile_vars,
  dispersions = tile_vars / (tile_means + 1),  # rough dispersion metric
  tile_idx = seq_len(nrow(tile_mat))
)

cat(sprintf("Var columns: %s\n", paste(colnames(var_df), collapse = ", ")))
cat(sprintf("tile stats: mean=[%.3f, %.3f], var=[%.3f, %.3f]\n",
            min(var_df$mean), max(var_df$mean),
            min(var_df$variance), max(var_df$variance)))

# ==============================================================================
# 5. Export to CSV/MTX format
# ==============================================================================

cat("\nExporting files...\n")

# tile accessibility matrix (sparse, Matrix Market format)
Matrix::writeMM(t(tile_mat),
                file = file.path(output_dir, "atac_tile_matrix.mtx"))
cat(sprintf("Wrote atac_tile_matrix.mtx (%d × %d)\n", ncol(tile_mat), nrow(tile_mat)))

# Row/col names for matrix
fwrite(data.table(tile = tile_names),
       file.path(output_dir, "atac_tile_names.csv"))
fwrite(data.table(cell_id = cell_ids),
       file.path(output_dir, "atac_cell_names.csv"))

# Coordinates
fwrite(coords_df, file.path(output_dir, "atac_coords.csv"))

# Obs
fwrite(obs_df, file.path(output_dir, "atac_obs.csv"))

# Var
fwrite(var_df, file.path(output_dir, "atac_var.csv"))

# tile genomic coordinates (parse from tile_names)
parse_tiles <- function(tile_str) {
  # Assuming format: chr1:10000-15000
  parts <- strsplit(tile_str, ":")[[1]]
  chr <- parts[1]
  coords <- strsplit(parts[2], "-")[[1]]
  start <- as.integer(coords[1])
  end <- as.integer(coords[2])
  data.table(tile = tile_str, chr = chr, start = start, end = end)
}

tile_coords <- rbindlist(lapply(tile_names, parse_tiles))
fwrite(tile_coords, file.path(output_dir, "atac_tile_coords.csv"))
cat(sprintf("tile coords: %d tiles parsed\n", nrow(tile_coords)))

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
