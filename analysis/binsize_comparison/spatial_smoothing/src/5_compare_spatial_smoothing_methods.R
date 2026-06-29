#!/usr/bin/env Rscript
# ============================================================================
# compare_spatial_smoothing_methods.R
#
# Compare 3 spatial-smoothing approaches on a single tissue:
# 1. Baseline (no smoothing) — original LSI, baseline clustering
# 2. Alpha-blend (current method) — k=6, 1-hop, alpha=0.5 blending
# 3. Iterative (Python-inspired) — k=8, 2-hop, self-inclusive, no blending
#
# For each method: cluster, UMAP, spatial plots, marker gene analysis.
# Output: single PDF with comparison (UMAPs, spatial scatter, marker heatmaps)
#         summary CSV with metrics table showing all 3 methods
#
# Usage:
#   Rscript compare_spatial_smoothing_methods.R <tissue>
#   Example: Rscript compare_spatial_smoothing_methods.R lowseq_489
#
# ============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(FNN)
  library(Matrix)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
})

log_msg <- function(tag, msg) {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), tag, msg))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  log_msg("error", "Usage: Rscript compare_spatial_smoothing_methods.R <tissue>")
  stop("Missing tissue argument")
}

tissue <- args[1]
log_msg("start", sprintf("===== Comparing 3 spatial smoothing methods: %s =====", tissue))

set.seed(42)
addArchRGenome("hg38")
addArchRThreads(threads = 8)

project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
proj_dir <- file.path(project_root, "analysis/binsize_comparison/archr_projects",
                      sprintf("%s_5000bp_binarizeFALSE", tissue))
spatial_coord_file <- file.path(project_root, "Data/01_inputs/spatial/tissue_positions_list.csv")
output_dir <- file.path(project_root, "analysis/binsize_comparison/spatial_smoothing")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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

if (!dir.exists(proj_dir)) {
  log_msg("error", sprintf("Project directory not found: %s", proj_dir))
  stop("Project missing")
}

log_msg("step", sprintf("Loading ArchR project from: %s", proj_dir))
proj <- loadArchRProject(path = proj_dir, force = TRUE)
log_msg("step", sprintf("Loaded project with %d cells", nrow(proj@cellColData)))

# ============================================================================
# Attach spatial coordinates
# ============================================================================
log_msg("step", "Attaching spatial coordinates...")

tissue_locs <- read.csv(spatial_coord_file)
tissue_locs <- tissue_locs[tissue_locs$in_tissue == 1, ]
tissue_locs$cellName <- paste0(sample_prefix, "#", tissue_locs$barcode, "-1")

# IMPORTANT: use match() (order-preserving) instead of merge() (which sorts
# rows by the join key by default). proj@embeddings / proj@reducedDims are
# never reordered, so cellColData must keep its original row order or every
# downstream positional index (cells_with_spatial, matSVD_full[idx,], etc.)
# silently points at the wrong cell.
match_idx <- match(rownames(proj@cellColData), tissue_locs$cellName)
proj@cellColData$x_spatial <- tissue_locs$x_spatial[match_idx]
proj@cellColData$y_spatial <- tissue_locs$y_spatial[match_idx]
proj@cellColData$array_row <- tissue_locs$array_row[match_idx]
proj@cellColData$array_col <- tissue_locs$array_col[match_idx]

n_matched <- sum(!is.na(proj@cellColData$x_spatial))
match_rate <- n_matched / nrow(proj@cellColData)
log_msg("step", sprintf("Spatial coordinate match rate: %d/%d cells (%.1f%%)",
                        n_matched, nrow(proj@cellColData), 100 * match_rate))

cells_with_spatial <- which(!is.na(proj@cellColData$x_spatial))
log_msg("step", sprintf("Found %d cells with spatial coordinates", length(cells_with_spatial)))

# ============================================================================
# Get spatial coordinates and baseline LSI (only for spatial cells, for kNN)
# ============================================================================
meta_df <- as.data.frame(proj@cellColData)
coords_spatial <- meta_df[cells_with_spatial, c("x_spatial", "y_spatial")]
coords_matrix <- as.matrix(coords_spatial)
n_spatial <- nrow(coords_matrix)

lsi_obj <- proj@reducedDims$IterativeLSI
matSVD_full <- lsi_obj$matSVD
matSVD_baseline <- matSVD_full  # Keep full for all cells

log_msg("step", sprintf("Full LSI: %d cells x %d dims; %d have spatial coordinates",
                        nrow(matSVD_full), ncol(matSVD_full), n_spatial))

# ============================================================================
# Distance-cutoff helper for spatial kNN weight matrices.
#
# get.knn() always returns exactly k matches regardless of true physical
# distance, so isolated/edge spots (no real neighbor within the tissue's
# normal grid pitch) get smoothed against distant, biologically irrelevant
# spots just to fill the k slots. Cutoff = 3x the median 1st-nearest-neighbor
# distance, computed per-tissue from the same coords_matrix (not a hardcoded
# constant) so it adapts to each tissue's actual spot pitch. Neighbors beyond
# the cutoff are dropped and remaining weights renormalized to sum to 1; a
# spot with zero neighbors within cutoff is left unsmoothed (self weight 1).
# ============================================================================
knn1 <- get.knn(coords_matrix, k = 1)
dist_cutoff <- 3 * median(knn1$nn.dist[, 1])
log_msg("step", sprintf("Spatial distance cutoff for smoothing: %.1f units (3x median 1-NN dist)", dist_cutoff))

build_distance_masked_weights <- function(nn_index, nn_dist, cutoff, n_spatial, include_self = FALSE) {
  k <- ncol(nn_index)
  i_list <- vector("list", n_spatial)
  j_list <- vector("list", n_spatial)
  x_list <- vector("list", n_spatial)
  for (i in 1:n_spatial) {
    keep <- nn_dist[i, ] <= cutoff
    js <- nn_index[i, keep]
    if (include_self) js <- c(i, js)
    if (length(js) == 0) js <- i  # fully isolated: leave unsmoothed
    w <- rep(1 / length(js), length(js))
    i_list[[i]] <- rep(i, length(js))
    j_list[[i]] <- js
    x_list[[i]] <- w
  }
  sparseMatrix(i = unlist(i_list), j = unlist(j_list), x = unlist(x_list),
               dims = c(n_spatial, n_spatial))
}

# ============================================================================
# Method B: Alpha-blend (k=6, 1-hop, alpha=0.5)
# Smooth only spatial cells, leave non-spatial unchanged
# ============================================================================
log_msg("step", "Method B: Alpha-blend (k=6, 1-hop, alpha=0.5)...")

knn_blend <- get.knn(coords_matrix, k = 6)
knn_idx_blend <- knn_blend$nn.index

# Build sparse matrix for spatial cells only, dropping neighbors beyond the
# distance cutoff and renormalizing remaining weights (see helper above).
W_blend <- build_distance_masked_weights(knn_idx_blend, knn_blend$nn.dist, dist_cutoff, n_spatial)

# Smooth only the spatial cells
matSVD_spatial_blend <- 0.5 * matSVD_full[cells_with_spatial, ] +
                       0.5 * as.matrix(W_blend %*% matSVD_full[cells_with_spatial, ])

# Create full matrix with smoothed spatial cells and original non-spatial
matSVD_blend <- matSVD_full
matSVD_blend[cells_with_spatial, ] <- matSVD_spatial_blend

log_msg("step", "Alpha-blend LSI computed")

# ============================================================================
# Method C: Iterative (k=8, 2-hop, self-inclusive, no blending)
# ============================================================================
log_msg("step", "Method C: Iterative (k=8, 2-hop, self-inclusive)...")

knn_iter <- get.knn(coords_matrix, k = 8)  # 8 distinct OTHER neighbors
# FNN::get.knn() always excludes self by construction (unlike Python's
# cKDTree.query against the same point set, which returns self at distance 0).
# Self is prepended inside build_distance_masked_weights() (include_self=TRUE)
# so the averaging window is self + up-to-8 distance-filtered neighbors,
# matching the reference Python implementation's k+1 query with col 0 = self.
W_iter <- build_distance_masked_weights(knn_iter$nn.index, knn_iter$nn.dist, dist_cutoff,
                                        n_spatial, include_self = TRUE)

# Smooth iteratively (2 rounds)
matSVD_spatial_iter <- matSVD_full[cells_with_spatial, ]
for (it in 1:2) {
  matSVD_spatial_iter <- as.matrix(W_iter %*% matSVD_spatial_iter)
}

# Create full matrix
matSVD_iterative <- matSVD_full
matSVD_iterative[cells_with_spatial, ] <- matSVD_spatial_iter

log_msg("step", "Iterative LSI computed (2 iterations)")

# ============================================================================
# Register all 3 LSI matrices as reducedDims
# ============================================================================
lsi_baseline <- lsi_obj
lsi_baseline$matSVD <- matSVD_baseline

lsi_blend <- lsi_obj
lsi_blend$matSVD <- matSVD_blend

lsi_iterative <- lsi_obj
lsi_iterative$matSVD <- matSVD_iterative

proj@reducedDims[["IterativeLSI_Blend"]] <- lsi_blend
proj@reducedDims[["IterativeLSI_Iterative"]] <- lsi_iterative
log_msg("step", "Registered all 3 LSI matrices")

# ============================================================================
# Cluster + UMAP for each method
# ============================================================================
log_msg("step", "Clustering and UMAP for all 3 methods...")

# Baseline: reuse existing on-disk Clusters/UMAP (matches archr_umap_cluster_comparison.pdf)
# instead of recomputing fresh, so baseline numbers are consistent across outputs.
# Only compute if missing (e.g. project never had clustering run).
if (!("Clusters" %in% colnames(proj@cellColData))) {
  log_msg("step", "No existing baseline Clusters found - computing fresh")
  proj <- addClusters(input = proj, reducedDims = "IterativeLSI",
                     method = "Seurat", name = "Clusters",
                     resolution = 0.8, force = TRUE)
} else {
  log_msg("step", "Reusing existing on-disk baseline Clusters")
}
if (!("UMAP" %in% names(proj@embeddings))) {
  log_msg("step", "No existing baseline UMAP found - computing fresh")
  proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI",
                 name = "UMAP", nNeighbors = 30, minDist = 0.5,
                 metric = "cosine", force = TRUE)
} else {
  log_msg("step", "Reusing existing on-disk baseline UMAP")
}

# Alpha-blend
proj <- addClusters(input = proj, reducedDims = "IterativeLSI_Blend",
                   method = "Seurat", name = "Clusters_Blend",
                   resolution = 0.8, force = TRUE)
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI_Blend",
               name = "UMAP_Blend", nNeighbors = 30, minDist = 0.5,
               metric = "cosine", force = TRUE)

# Iterative (use lower resolution to avoid too many clusters)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI_Iterative",
                   method = "Seurat", name = "Clusters_Iterative",
                   resolution = 0.4, force = TRUE)
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI_Iterative",
               name = "UMAP_Iterative", nNeighbors = 30, minDist = 0.5,
               metric = "cosine", force = TRUE)

log_msg("step", "Clustering and UMAP complete for all 3 methods")

# ============================================================================
# Compute spatial coherence for each method (only spatial cells)
# ============================================================================
log_msg("step", "Computing spatial coherence metrics...")

compute_coherence <- function(clusters, knn_idx) {
  if (is.null(clusters) || is.null(knn_idx) || length(clusters) != nrow(knn_idx)) {
    return(NA)
  }
  coherence <- numeric(length(clusters))
  for (i in 1:length(clusters)) {
    neighbor_clusters <- clusters[knn_idx[i, ]]
    same_as_self <- sum(neighbor_clusters == clusters[i], na.rm = TRUE)
    coherence[i] <- same_as_self / ncol(knn_idx)
  }
  result <- mean(coherence, na.rm = TRUE)
  if (is.na(result) || !is.finite(result)) return(NA)
  result
}

knn_coh <- get.knn(coords_matrix, k = 6)
knn_idx_coh <- knn_coh$nn.index

# Clusters only for spatial cells (extract properly from DFrame, preserve names)
meta_full <- as.data.frame(proj@cellColData)
all_clusters_baseline <- meta_full$Clusters
all_clusters_blend <- meta_full$Clusters_Blend
all_clusters_iterative <- meta_full$Clusters_Iterative

# Convert to character to preserve cluster names (C1, C2, etc)
clusters_baseline_spatial <- as.character(all_clusters_baseline[cells_with_spatial])
clusters_blend_spatial <- as.character(all_clusters_blend[cells_with_spatial])
clusters_iterative_spatial <- as.character(all_clusters_iterative[cells_with_spatial])

coh_baseline <- compute_coherence(clusters_baseline_spatial, knn_idx_coh)
coh_blend <- compute_coherence(clusters_blend_spatial, knn_idx_coh)
coh_iterative <- compute_coherence(clusters_iterative_spatial, knn_idx_coh)

log_msg("step", sprintf("Coherence — Baseline: %.4f, Blend: %.4f, Iterative: %.4f",
                        coh_baseline, coh_blend, coh_iterative))

# ============================================================================
# Marker gene analysis for each method
# ============================================================================
log_msg("step", "Running marker gene analysis for all 3 methods...")

run_marker_analysis <- function(proj, cluster_col, method_name) {
  tryCatch({
    n_clusters <- length(unique(proj@cellColData[[cluster_col]]))

    # Skip if too many clusters (likely fragments, hard to compute markers)
    if (n_clusters > 15) {
      log_msg("warn", sprintf("Marker analysis skipped for %s: %d clusters (too many, likely fragments)",
                             method_name, n_clusters))
      return(list(markers = NULL, marker_list = NULL, n_clusters = n_clusters,
                 total_markers = NA, mean_lfc = NA))
    }

    markers <- getMarkerFeatures(proj, useMatrix = "GeneScoreMatrix",
                                groupBy = cluster_col,
                                bias = c("TSSEnrichment", "log10(nFrags)"),
                                testMethod = "wilcoxon")

    if (is.null(markers)) {
      log_msg("warn", sprintf("Marker analysis: NULL result for %s", method_name))
      return(list(markers = NULL, marker_list = NULL, n_clusters = n_clusters,
                 total_markers = NA, mean_lfc = NA))
    }

    marker_list <- getMarkers(markers, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")

    total_markers <- sum(lengths(marker_list))
    lfc_values <- unlist(lapply(marker_list, function(x) {
      if (is.null(x) || nrow(x) == 0) return(numeric(0))
      abs(x$Log2FC)
    }))

    mean_lfc <- if (length(lfc_values) > 0) mean(lfc_values, na.rm = TRUE) else NA

    list(
      markers = markers,
      marker_list = marker_list,
      n_clusters = n_clusters,
      total_markers = total_markers,
      mean_lfc = mean_lfc
    )
  }, error = function(e) {
    log_msg("warn", sprintf("Marker analysis error for %s: %s", method_name, e$message))
    list(markers = NULL, marker_list = NULL, n_clusters = NA, total_markers = NA, mean_lfc = NA)
  })
}

markers_baseline <- run_marker_analysis(proj, "Clusters", "Baseline")
markers_blend <- run_marker_analysis(proj, "Clusters_Blend", "Alpha-blend")
markers_iterative <- run_marker_analysis(proj, "Clusters_Iterative", "Iterative")

log_msg("step", "Marker gene analysis complete")

# ============================================================================
# Generate comparison plots (spatial cells only)
# ============================================================================
log_msg("step", "Generating comparison plots...")

# Remove NA values for color palette computation
unique_baseline <- unique(na.omit(clusters_baseline_spatial))
unique_blend <- unique(na.omit(clusters_blend_spatial))
unique_iterative <- unique(na.omit(clusters_iterative_spatial))

n_colors <- max(length(unique_baseline), length(unique_blend), length(unique_iterative))
n_colors <- max(n_colors, 3)  # Ensure at least 3 colors

# Use discrete color palette (Set1 for up to 9, then polychrome for larger)
if (n_colors <= 9) {
  color_pal <- RColorBrewer::brewer.pal(n_colors, "Set1")
} else if (n_colors <= 12) {
  color_pal <- RColorBrewer::brewer.pal(12, "Set3")
} else {
  color_pal <- scales::hue_pal(l = 60)(n_colors)  # Better saturation for discrete
}

# Extract embeddings for spatial cells only
umap_baseline <- proj@embeddings$UMAP$df[cells_with_spatial, 1:2]
umap_blend <- proj@embeddings$UMAP_Blend$df[cells_with_spatial, 1:2]
umap_iterative <- proj@embeddings$UMAP_Iterative$df[cells_with_spatial, 1:2]

spatial_x <- meta_df$x_spatial[cells_with_spatial]
spatial_y <- meta_df$y_spatial[cells_with_spatial]

# Order cluster labels numerically (C1, C2, ..., C10) instead of the default
# factor() lexicographic sort (C1, C10, C11, C2, ...)
order_clusters_factor <- function(clusters) {
  lvls <- paste0("C", sort(unique(as.integer(gsub("[^0-9]", "", clusters)))))
  factor(clusters, levels = lvls)
}

# Build ggplot UMAP panels
plot_umap <- function(embedding, clusters, title, color_pal) {
  df <- data.frame(UMAP1 = embedding[, 1], UMAP2 = embedding[, 2], Cluster = order_clusters_factor(clusters))
  ggplot(df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
    geom_point(size = 0.5, alpha = 0.8) +
    scale_color_manual(values = color_pal[1:nlevels(df$Cluster)]) +
    labs(title = title) +
    theme_minimal() +
    coord_fixed()
}

p_umap_baseline <- plot_umap(umap_baseline, clusters_baseline_spatial, "Baseline (No Smoothing)", color_pal)
p_umap_blend <- plot_umap(umap_blend, clusters_blend_spatial, "Alpha-Blend (k=6)", color_pal)
p_umap_iterative <- plot_umap(umap_iterative, clusters_iterative_spatial, "Iterative (k=8, 2-hop)", color_pal)

# Build ggplot spatial panels
plot_spatial <- function(x, y, clusters, title, color_pal) {
  df <- data.frame(X = x, Y = y, Cluster = order_clusters_factor(clusters))
  ggplot(df, aes(x = X, y = Y, color = Cluster)) +
    geom_point(size = 0.5, alpha = 0.8) +
    scale_color_manual(values = color_pal[1:nlevels(df$Cluster)]) +
    labs(title = title, x = "X spatial", y = "Y spatial") +
    theme_minimal() +
    coord_fixed(ratio = 1)
}

p_spatial_baseline <- plot_spatial(spatial_x, spatial_y, clusters_baseline_spatial,
                                   "Baseline Spatial", color_pal)
p_spatial_blend <- plot_spatial(spatial_x, spatial_y, clusters_blend_spatial,
                               "Alpha-Blend Spatial", color_pal)
p_spatial_iterative <- plot_spatial(spatial_x, spatial_y, clusters_iterative_spatial,
                                   "Iterative Spatial", color_pal)

# ============================================================================
# Confusion heatmaps between the 3 methods (row-normalized %, numeric cluster order)
# ============================================================================
log_msg("step", "Building method-vs-method confusion heatmaps...")

plot_confusion_heatmap <- function(clusters_a, clusters_b, label_a, label_b) {
  fa <- order_clusters_factor(clusters_a)
  fb <- order_clusters_factor(clusters_b)
  tab <- table(fa, fb)
  prop <- prop.table(tab, margin = 1) * 100
  df <- as.data.frame(prop)
  colnames(df) <- c("A", "B", "Pct")
  df$A <- factor(df$A, levels = levels(fa))
  df$B <- factor(df$B, levels = levels(fb))
  ggplot(df, aes(x = B, y = A, fill = Pct)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(Pct >= 1, sprintf("%.0f", Pct), "")), size = 2.8, color = "black") +
    scale_fill_viridis_c(option = "viridis", limits = c(0, 100), name = "% of row") +
    labs(title = sprintf("%s vs %s", label_a, label_b), x = label_b, y = label_a) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, face = "bold"))
}

p_conf_base_blend <- plot_confusion_heatmap(clusters_baseline_spatial, clusters_blend_spatial, "Baseline", "Alpha-Blend")
p_conf_base_iter <- plot_confusion_heatmap(clusters_baseline_spatial, clusters_iterative_spatial, "Baseline", "Iterative")
p_conf_blend_iter <- plot_confusion_heatmap(clusters_blend_spatial, clusters_iterative_spatial, "Alpha-Blend", "Iterative")

confusion_pdf <- file.path(output_dir, sprintf("%s_5000bp_method_confusion_heatmaps.pdf", tissue))
ggsave(confusion_pdf,
      (p_conf_base_blend | p_conf_base_iter | p_conf_blend_iter) +
        patchwork::plot_annotation(title = sprintf("%s: cluster reassignment across smoothing methods (row-normalized %%)", tissue)),
      width = 16, height = 5.5)
log_msg("step", sprintf("Saved confusion heatmaps: %s", confusion_pdf))

best_match_pct <- function(clusters_a, clusters_b) {
  tab <- table(clusters_a, clusters_b)
  round(100 * sum(apply(tab, 1, max)) / sum(tab), 1)
}
agree_base_blend <- best_match_pct(clusters_baseline_spatial, clusters_blend_spatial)
agree_base_iter <- best_match_pct(clusters_baseline_spatial, clusters_iterative_spatial)
agree_blend_iter <- best_match_pct(clusters_blend_spatial, clusters_iterative_spatial)

# Per-cell cluster assignment CSV (for reproducible confusion-matrix rebuilding)
assign_df <- data.frame(
  cellName = rownames(proj@cellColData)[cells_with_spatial],
  Baseline = clusters_baseline_spatial,
  Alpha_Blend = clusters_blend_spatial,
  Iterative = clusters_iterative_spatial
)
assign_csv <- file.path(output_dir, sprintf("%s_5000bp_3method_cluster_assignments.csv", tissue))
write.csv(assign_df, assign_csv, row.names = FALSE)
log_msg("step", sprintf("Saved per-cell cluster assignments: %s", assign_csv))

# Save PDF
pdf_file <- file.path(output_dir, sprintf("%s_5000bp_3method_comparison.pdf", tissue))
pdf(pdf_file, width = 18, height = 6)

# Page 1: UMAP comparison
print(p_umap_baseline | p_umap_blend | p_umap_iterative)

# Page 2: Spatial comparison
print(p_spatial_baseline | p_spatial_blend | p_spatial_iterative)

# Page 3: Marker heatmaps (if available)
if (!is.null(markers_baseline$markers)) {
  tryCatch({
    h1 <- plotMarkerHeatmap(seMarker = markers_baseline$markers, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
    print(h1)
  }, error = function(e) {
    log_msg("warn", sprintf("Could not plot baseline heatmap: %s", e$message))
  })
}

if (!is.null(markers_blend$markers)) {
  tryCatch({
    h2 <- plotMarkerHeatmap(seMarker = markers_blend$markers, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
    print(h2)
  }, error = function(e) {
    log_msg("warn", sprintf("Could not plot blend heatmap: %s", e$message))
  })
}

if (!is.null(markers_iterative$markers)) {
  tryCatch({
    h3 <- plotMarkerHeatmap(seMarker = markers_iterative$markers, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
    print(h3)
  }, error = function(e) {
    log_msg("warn", sprintf("Could not plot iterative heatmap: %s", e$message))
  })
}

dev.off()
log_msg("step", sprintf("Saved comparison PDF: %s", pdf_file))

# ============================================================================
# Summary CSV
# ============================================================================
summary_df <- data.frame(
  Method = c("Baseline (No Smoothing)", "Alpha-Blend (k=6, 1-hop, alpha=0.5)", "Iterative (k=8, 2-hop, self-incl)"),
  N_Clusters = c(markers_baseline$n_clusters, markers_blend$n_clusters, markers_iterative$n_clusters),
  Total_Markers = c(markers_baseline$total_markers, markers_blend$total_markers, markers_iterative$total_markers),
  Mean_Log2FC = c(markers_baseline$mean_lfc, markers_blend$mean_lfc, markers_iterative$mean_lfc),
  Spatial_Coherence = c(coh_baseline, coh_blend, coh_iterative)
)

csv_file <- file.path(output_dir, sprintf("%s_5000bp_3method_comparison.csv", tissue))
write.csv(summary_df, csv_file, row.names = FALSE)
log_msg("step", sprintf("Saved comparison CSV: %s", csv_file))

# Summary text
summary_file <- file.path(output_dir, sprintf("%s_5000bp_3method_comparison_summary.txt", tissue))
sink(summary_file)
cat(sprintf("===== 3-Method Spatial Smoothing Comparison: %s =====\n\n", tissue))
cat(sprintf("Tissue: %s\nCells (total / with spatial coords): %d / %d\n\n", tissue, nrow(proj@cellColData), n_spatial))
cat("METHOD COMPARISON\n")
cat("=================\n\n")
print(summary_df)
cat("\n\nINTERPRETATION\n")
cat("==============\n")
cat("Spatial_Coherence: fraction of spatial neighbors with same cluster label (higher = more contiguous)\n")
cat("Total_Markers: count of significant markers (FDR<=0.05, Log2FC>=0.5) across all clusters\n")
cat("Mean_Log2FC: mean absolute Log2 fold-change of markers (higher = stronger separation)\n")
cat("N_Clusters: number of clusters detected at resolution 0.8\n")
cat("\n\nMETHOD CLUSTER AGREEMENT (best-match overlap %, cluster IDs are arbitrary per-run)\n")
cat("==============================================================================\n")
cat(sprintf("Baseline    vs Alpha-Blend: %.1f%%\n", agree_base_blend))
cat(sprintf("Baseline    vs Iterative:   %.1f%%\n", agree_base_iter))
cat(sprintf("Alpha-Blend vs Iterative:   %.1f%%\n", agree_blend_iter))
sink()

log_msg("step", sprintf("Saved summary: %s", summary_file))

log_msg("done", sprintf("Completed: %s (3-method comparison)", tissue))
