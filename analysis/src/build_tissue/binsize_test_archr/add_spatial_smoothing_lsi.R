#!/usr/bin/env Rscript
# ============================================================================
# add_spatial_smoothing_lsi.R
#
# Load an existing 5000bp ArchR project, attach spatial coordinates,
# build a spatial kNN graph, apply spatial smoothing to the LSI embedding,
# re-cluster and re-embed, and generate before/after comparison artifacts.
#
# Usage:
#   Rscript add_spatial_smoothing_lsi.R <tissue> [k=6] [alpha=0.5] [output_dir]
#   Example: Rscript add_spatial_smoothing_lsi.R lowseq_489
#
# ============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(FNN)
  library(Matrix)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

# Logging function
log_msg <- function(tag, msg) {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), tag, msg))
}

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  log_msg("error", "Usage: Rscript add_spatial_smoothing_lsi.R <tissue> [k=6] [alpha=0.5]")
  stop("Missing tissue argument")
}

tissue <- args[1]
k <- if (length(args) >= 2) as.integer(args[2]) else 6
alpha <- if (length(args) >= 3) as.numeric(args[3]) else 0.5

log_msg("start", sprintf("===== Spatial smoothing for: %s, k=%d, alpha=%.2f =====",
                         tissue, k, alpha))

# MUST SET THESE BEFORE CREATING/LOADING ArchRProject
set.seed(42)
addArchRGenome("hg38")
addArchRThreads(threads = 8)

# Setup paths
project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
proj_dir <- file.path(project_root, "analysis/binsize_comparison/archr_projects",
                      sprintf("%s_5000bp_binarizeFALSE", tissue))
spatial_coord_file <- file.path(project_root, "Data/01_inputs/spatial/tissue_positions_list.csv")
output_dir <- file.path(project_root, "analysis/binsize_comparison/spatial_smoothing")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Tissue-specific sample name prefix (exact casing from 0_create_archr_no_doublet.R)
tissue_name_map <- list(
  "lowseq_489" = "Lowseq_489",
  "lowseq_488B" = "Lowseq_488B",
  "deepseq_488B" = "Deepseq_488B",
  "deepseq_489" = "Deepseq_489"
)

sample_prefix <- tissue_name_map[[tissue]]
if (is.null(sample_prefix)) {
  log_msg("error", sprintf("Unknown tissue: %s", tissue))
  stop("Invalid tissue name")
}

# Verify project exists
if (!dir.exists(proj_dir)) {
  log_msg("error", sprintf("Project directory not found: %s", proj_dir))
  stop("Project missing")
}

log_msg("step", sprintf("Loading ArchR project from: %s", proj_dir))
tryCatch({
  proj <- loadArchRProject(path = proj_dir, force = TRUE)
  log_msg("step", sprintf("Loaded project with %d cells", nrow(proj@cellColData)))
}, error = function(e) {
  log_msg("error", sprintf("Failed to load project: %s", e$message))
  stop(e)
})

# ============================================================================
# Step 1: Attach spatial coordinates
# ============================================================================
log_msg("step", "Attaching spatial coordinates...")

tryCatch({
  tissue_locs <- read.csv(spatial_coord_file)

  # Filter to in_tissue == 1
  tissue_locs <- tissue_locs[tissue_locs$in_tissue == 1, ]

  # Create ArchR-compatible barcodes (format: "SampleName#barcode-1")
  tissue_locs$cellName <- paste0(sample_prefix, "#", tissue_locs$barcode, "-1")

  # Get current cellColData
  proj_meta <- data.frame(proj@cellColData)
  proj_meta$cellName <- rownames(proj_meta)

  # Merge on cellName
  proj_meta_spatial <- merge(proj_meta, tissue_locs,
                            by.x = "cellName", by.y = "cellName",
                            all.x = TRUE)

  # Set rownames and update cellColData
  rownames(proj_meta_spatial) <- proj_meta_spatial$cellName
  proj_meta_spatial$cellName <- NULL
  proj@cellColData <- as(proj_meta_spatial, "DFrame")

  # Check match rate
  n_matched <- sum(!is.na(proj@cellColData$x_spatial))
  match_rate <- n_matched / nrow(proj@cellColData)

  log_msg("step", sprintf("Spatial coordinate match rate: %d/%d cells (%.1f%%)",
                          n_matched, nrow(proj@cellColData), 100 * match_rate))

  if (match_rate < 0.8) {
    log_msg("warn", "Low match rate (<80%); check spatial coordinate join for this tissue")
  }
}, error = function(e) {
  log_msg("error", sprintf("Failed to attach spatial coordinates: %s", e$message))
  stop(e)
})

# Subset to cells with spatial coordinates
cells_with_spatial <- which(!is.na(proj@cellColData$x_spatial))
log_msg("step", sprintf("Subsetting to %d cells with spatial coordinates", length(cells_with_spatial)))

if (length(cells_with_spatial) == 0) {
  log_msg("error", "No cells have spatial coordinates")
  stop("Spatial join failed")
}

# ============================================================================
# Step 2: Build spatial kNN graph
# ============================================================================
log_msg("step", sprintf("Building spatial kNN graph with k=%d...", k))

tryCatch({
  coords <- proj@cellColData[cells_with_spatial, c("x_spatial", "y_spatial")]
  coords_matrix <- as.matrix(coords)

  # FNN::get.knn returns list with idx (cell indices) and dist
  knn_result <- get.knn(coords_matrix, k = k)
  knn_idx <- knn_result$nn.idx

  # Build sparse row-normalized weight matrix
  # Each row i has uniform weights 1/k for its k neighbors
  n_cells <- nrow(coords_matrix)
  i_idx <- rep(1:n_cells, each = k)
  j_idx <- as.vector(t(knn_idx))  # Transpose to vectorize column-wise

  # Create sparse matrix
  W <- sparseMatrix(i = i_idx, j = j_idx, x = rep(1/k, length(i_idx)),
                    dims = c(n_cells, n_cells))

  log_msg("step", sprintf("Built sparse weight matrix: %d x %d", nrow(W), ncol(W)))
}, error = function(e) {
  log_msg("error", sprintf("Failed to build kNN graph: %s", e$message))
  stop(e)
})

# ============================================================================
# Step 3: Extract and smooth LSI embedding
# ============================================================================
log_msg("step", "Extracting LSI embedding...")

tryCatch({
  # Get original LSI matrix from the full project
  # (not subsetted yet)
  lsi_obj <- proj@reducedDims$IterativeLSI
  if (is.null(lsi_obj)) {
    log_msg("error", "IterativeLSI not found in project")
    stop("LSI missing")
  }

  # Extract the cell x dims matrix
  # In ArchR, this is stored as matSVD (cells x dimsToUse)
  matSVD_full <- lsi_obj$matSVD
  log_msg("step", sprintf("Original LSI matrix: %d cells x %d dims",
                          nrow(matSVD_full), ncol(matSVD_full)))

  # Subset to cells with spatial coordinates
  # cells_with_spatial are indices into the full project
  matSVD_spatial <- matSVD_full[cells_with_spatial, , drop = FALSE]
  log_msg("step", sprintf("Subset LSI matrix: %d cells x %d dims",
                          nrow(matSVD_spatial), ncol(matSVD_spatial)))

  # Apply spatial smoothing
  # smoothed = alpha * original + (1 - alpha) * (neighbor average)
  matSVD_smoothed <- alpha * matSVD_spatial +
                     (1 - alpha) * (as.matrix(W %*% matSVD_spatial))

  log_msg("step", sprintf("Applied spatial smoothing: alpha=%.2f", alpha))

  # Sanity checks
  if (any(is.na(matSVD_smoothed))) {
    log_msg("warn", "NAs detected in smoothed LSI matrix; replacing with original")
    matSVD_smoothed[is.na(matSVD_smoothed)] <- matSVD_spatial[is.na(matSVD_smoothed)]
  }

}, error = function(e) {
  log_msg("error", sprintf("Failed to smooth LSI: %s", e$message))
  stop(e)
})

# ============================================================================
# Step 4: Re-subset project to spatial cells and register smoothed LSI
# ============================================================================
log_msg("step", "Re-subsetting project to cells with spatial coordinates...")

# Create a new ArchR project with only the spatial cells
# This is needed because the rest of the pipeline expects consistent dimensions
proj <- subsetArchRProject(proj, cells = rownames(proj@cellColData)[cells_with_spatial],
                          outputDirectory = proj_dir, force = TRUE)
log_msg("step", sprintf("Subset project now has %d cells", nrow(proj@cellColData)))

# Register smoothed LSI as a new reducedDims slot
log_msg("step", "Registering smoothed LSI...")

tryCatch({
  # Create a new reducedDims object mirroring IterativeLSI
  lsi_smooth <- list(
    matSVD = matSVD_smoothed,
    scaleDims = lsi_obj$scaleDims,
    outDir = lsi_obj$outDir,
    date = Sys.time(),
    useMatrix = lsi_obj$useMatrix
  )

  proj@reducedDims[["IterativeLSI_SpatialSmooth"]] <- lsi_smooth
  log_msg("step", "Registered IterativeLSI_SpatialSmooth")

}, error = function(e) {
  log_msg("error", sprintf("Failed to register smoothed LSI: %s", e$message))
  stop(e)
})

# ============================================================================
# Step 5: Re-cluster on smoothed LSI
# ============================================================================
log_msg("step", "Re-clustering on smoothed LSI...")

tryCatch({
  proj <- addClusters(
    input = proj,
    reducedDims = "IterativeLSI_SpatialSmooth",
    method = "Seurat",
    name = "Clusters_SpatialSmooth",
    resolution = 0.8,
    force = TRUE
  )
  log_msg("step", sprintf("Clustering complete: %d clusters",
                          length(unique(proj@cellColData$Clusters_SpatialSmooth))))
}, error = function(e) {
  log_msg("error", sprintf("Clustering failed: %s", e$message))
  stop(e)
})

# ============================================================================
# Step 6: Re-embed with UMAP on smoothed LSI
# ============================================================================
log_msg("step", "Computing UMAP on smoothed LSI...")

tryCatch({
  proj <- addUMAP(
    ArchRProj = proj,
    reducedDims = "IterativeLSI_SpatialSmooth",
    name = "UMAP_SpatialSmooth",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine",
    force = TRUE
  )
  log_msg("step", "UMAP computation complete")
}, error = function(e) {
  log_msg("error", sprintf("UMAP failed: %s", e$message))
  stop(e)
})

# ============================================================================
# Step 7: Compute spatial coherence metric (before/after)
# ============================================================================
log_msg("step", "Computing spatial coherence metric...")

tryCatch({
  # For each cell, compute fraction of spatial neighbors with same cluster label
  compute_coherence <- function(clusters, knn_idx) {
    coherence <- numeric(length(clusters))
    for (i in 1:length(clusters)) {
      neighbor_clusters <- clusters[knn_idx[i, ]]
      same_as_self <- sum(neighbor_clusters == clusters[i])
      coherence[i] <- same_as_self / length(neighbor_clusters)
    }
    mean(coherence)
  }

  coh_before <- compute_coherence(as.numeric(proj@cellColData$Clusters), knn_result$nn.idx)
  coh_after <- compute_coherence(as.numeric(proj@cellColData$Clusters_SpatialSmooth),
                                 knn_result$nn.idx)

  log_msg("step", sprintf("Coherence before: %.4f, after: %.4f (delta: +%.4f)",
                          coh_before, coh_after, coh_after - coh_before))

}, error = function(e) {
  log_msg("warn", sprintf("Failed to compute coherence: %s", e$message))
  coh_before <- NA
  coh_after <- NA
})

# ============================================================================
# Step 8: Export CSVs and PDFs
# ============================================================================
log_msg("step", "Exporting comparison data...")

tryCatch({
  # Extract embeddings
  umap_df <- data.frame(
    cellID = rownames(proj@cellColData),
    x_spatial = proj@cellColData$x_spatial,
    y_spatial = proj@cellColData$y_spatial,
    Clusters = proj@cellColData$Clusters,
    UMAP_1 = proj@embeddings$UMAP$df[, 1],
    UMAP_2 = proj@embeddings$UMAP$df[, 2],
    Clusters_SpatialSmooth = proj@cellColData$Clusters_SpatialSmooth,
    UMAP_1_smooth = proj@embeddings$UMAP_SpatialSmooth$df[, 1],
    UMAP_2_smooth = proj@embeddings$UMAP_SpatialSmooth$df[, 2],
    stringsAsFactors = FALSE
  )

  # Write CSV
  csv_file <- file.path(output_dir,
                       sprintf("%s_5000bp_spatial_smoothing_comparison.csv", tissue))
  write.csv(umap_df, csv_file, row.names = FALSE)
  log_msg("step", sprintf("Saved comparison CSV: %s", csv_file))

  # Generate PDF with 4 panels
  pdf_file <- file.path(output_dir,
                       sprintf("%s_5000bp_spatial_smoothing_plots.pdf", tissue))

  pdf(pdf_file, width = 14, height = 12)

  # Panel 1: UMAP before (original clusters)
  p1 <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "cellColData",
    name = "Clusters",
    embedding = "UMAP",
    size = 1.5
  )
  p1 <- p1 + ggtitle(sprintf("UMAP Before Spatial Smoothing (Clusters)"))
  print(p1)

  # Panel 2: UMAP after (smoothed clusters)
  p2 <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "cellColData",
    name = "Clusters_SpatialSmooth",
    embedding = "UMAP_SpatialSmooth",
    size = 1.5
  )
  p2 <- p2 + ggtitle(sprintf("UMAP After Spatial Smoothing (Clusters_SpatialSmooth)"))
  print(p2)

  # Panel 3: Spatial scatter before
  p3 <- ggplot(umap_df, aes(x = x_spatial, y = y_spatial, color = factor(Clusters))) +
    geom_point(size = 1.5) +
    labs(title = "Spatial Coordinates (Original Clusters)",
         x = "X spatial", y = "Y spatial", color = "Cluster") +
    theme_minimal()
  print(p3)

  # Panel 4: Spatial scatter after
  p4 <- ggplot(umap_df, aes(x = x_spatial, y = y_spatial, color = factor(Clusters_SpatialSmooth))) +
    geom_point(size = 1.5) +
    labs(title = "Spatial Coordinates (Spatial Smooth Clusters)",
         x = "X spatial", y = "Y spatial", color = "Cluster") +
    theme_minimal()
  print(p4)

  dev.off()
  log_msg("step", sprintf("Saved plots PDF: %s", pdf_file))

  # Write text summary
  summary_file <- file.path(output_dir,
                           sprintf("%s_5000bp_spatial_smoothing_summary.txt", tissue))

  sink(summary_file)
  cat(sprintf("===== Spatial Smoothing Summary: %s =====\n\n", tissue))
  cat(sprintf("Configuration:\n"))
  cat(sprintf("  k (spatial neighbors): %d\n", k))
  cat(sprintf("  alpha (smoothing weight): %.2f\n", alpha))
  cat(sprintf("  Input: 5000bp, non-binarized\n"))
  cat(sprintf("  Cells: %d\n\n", nrow(proj@cellColData)))

  cat(sprintf("Clustering Results:\n"))
  cat(sprintf("  Clusters (original): %d clusters\n",
              length(unique(proj@cellColData$Clusters))))
  cat(sprintf("  Clusters_SpatialSmooth: %d clusters\n",
              length(unique(proj@cellColData$Clusters_SpatialSmooth))))

  cat(sprintf("\nSpatial Coherence Metric (avg fraction of neighbors with same label):\n"))
  cat(sprintf("  Before spatial smoothing: %.4f\n", coh_before))
  cat(sprintf("  After spatial smoothing: %.4f\n", coh_after))
  cat(sprintf("  Improvement: +%.4f (%.2f%%)\n\n",
              coh_after - coh_before, 100 * (coh_after - coh_before) / coh_before))

  cat(sprintf("Outputs:\n"))
  cat(sprintf("  CSV: %s\n", csv_file))
  cat(sprintf("  PDF: %s\n", pdf_file))
  cat(sprintf("  This file: %s\n", summary_file))
  sink()

  log_msg("step", sprintf("Saved summary: %s", summary_file))

}, error = function(e) {
  log_msg("error", sprintf("Failed to export: %s", e$message))
  stop(e)
})

# ============================================================================
# Step 9: Save project
# ============================================================================
log_msg("step", "Saving project...")
tryCatch({
  proj <- saveArchRProject(proj)
  log_msg("step", "Project saved")
}, error = function(e) {
  log_msg("warn", sprintf("Failed to save project: %s", e$message))
})

log_msg("done", sprintf("Completed successfully: %s (k=%d, alpha=%.2f)", tissue, k, alpha))
