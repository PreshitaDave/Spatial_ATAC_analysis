#!/usr/bin/env Rscript
# 9_extract_normal_barcodes_calicost.R
#
# Identify normal (immune/stromal) spots in a spatial ATAC-seq tissue using
# ArchR cluster assignments + module gene scores, then export barcodes in the
# format expected by CalicoST normalidx_file.
#
# Normal cell identification follows the same approach as the tissue-specific
# Rmd (e.g., analysis/src/build_tissue/tissue_specific/lowseq_489.Rmd):
#   - 4 module scores: BScore (B cells), TScore (T cells),
#                      MScore (Myeloid), FScore (Fibroblast/CAF)
#   - Clusters enriched for immune/stromal markers are labelled "normal"
#   - Threshold: composite score > median + 1 SD across clusters
#
# Outputs:
#   analysis/binsize_comparison/normal_barcodes/<tissue>_normal_barcodes.csv
#   analysis/binsize_comparison/normal_barcodes/<tissue>_module_scores_umap.pdf
#   analysis/binsize_comparison/normal_barcodes/<tissue>_module_score_by_cluster.pdf
#   analysis/binsize_comparison/normal_barcodes/<tissue>_normal_tumor_spatial.pdf
#
# Usage:
#   Rscript 9_extract_normal_barcodes_calicost.R <tissue>
#   Example: Rscript 9_extract_normal_barcodes_calicost.R lowseq_489
#
# CalicoST normalidx_file format:
#   Single-column CSV, no header, barcode strings matching parsed_inputs/
#   barcodes (format: <SampleName>#<barcode>-1)

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# ── Args ──────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
tissue <- if (length(args) >= 1) args[1] else "lowseq_489"
cat(sprintf("[9_extract_normal_barcodes] Tissue: %s\n", tissue))

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_ROOT <- "/projectnb/paxlab/presh/projects/spatial_atac"

ARCHR_PATHS <- list(
  lowseq_489  = file.path(PROJECT_ROOT, "Data/01_outputs/archR_objects/lowseq_489/lowseq_489_archR_project_final"),
  lowseq_488B = file.path(PROJECT_ROOT, "Data/01_outputs/archR_objects/lowseq_488B/lowseq_488B_archR_project_final"),
  deepseq_488B = file.path(PROJECT_ROOT, "Data/01_outputs/archR_objects/deepseq_488B/deepseq_488B_archR_project_final"),
  deepseq_489 = file.path(PROJECT_ROOT, "Data/01_outputs/archR_objects/deepseq_489/deepseq_489_archR_project_final")
)

SPATIAL_PATHS <- list(
  lowseq_489  = file.path(PROJECT_ROOT, "Data/01_inputs/spatial/lowseq_489/spatial/tissue_positions_list.csv"),
  lowseq_488B = file.path(PROJECT_ROOT, "Data/01_inputs/spatial/lowseq_488B/spatial/tissue_positions_list.csv"),
  deepseq_488B = file.path(PROJECT_ROOT, "Data/01_inputs/spatial/deepseq_488B/spatial/tissue_positions_list.csv"),
  deepseq_489 = file.path(PROJECT_ROOT, "Data/01_inputs/spatial/deepseq_489/spatial/tissue_positions_list.csv")
)

archr_path   <- ARCHR_PATHS[[tissue]]
spatial_path <- SPATIAL_PATHS[[tissue]]
out_dir      <- file.path(PROJECT_ROOT, "analysis/binsize_comparison/normal_barcodes")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (is.null(archr_path) || !dir.exists(archr_path)) {
  stop(sprintf("ArchR project not found for tissue '%s': %s", tissue, archr_path))
}

addArchRThreads(threads = 4)
addArchRGenome("hg38")

# ── Load ArchR project ────────────────────────────────────────────────────────
cat("[9_extract_normal_barcodes] Loading ArchR project...\n")
proj <- loadArchRProject(path = archr_path, showLogo = FALSE)
cat(sprintf("  Loaded: %d cells, clusters: %s\n",
            nCells(proj),
            paste(sort(unique(proj$Clusters)), collapse = ", ")))

# ── Compute module scores manually from GeneScoreMatrix ──────────────────────
# Avoids addModuleScore() which has an Rle dimension mismatch bug in this ArchR version.
# Same gene sets as lowseq_489.Rmd
gene_groups <- list(
  BScore = c("MS4A1", "CD79A", "CD74", "CD19", "PAX5"),        # B cells
  TScore = c("CD3D", "CD8A", "GZMB", "CCR7", "LEF1"),          # T cells
  MScore = c("CD68", "LYZ", "CXCL10", "CD163"),                # Myeloid/Macrophage
  FScore = c("COL1A1", "COL1A2", "VIM", "DCN", "PDGFRA")      # Fibroblast/CAF
)

cat("[9_extract_normal_barcodes] Extracting GeneScoreMatrix...\n")
gsm_se <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix", verbose = FALSE)
gsm    <- assay(gsm_se, "GeneScoreMatrix")   # genes × cells (sparse)
gene_names <- rowData(gsm_se)$name

# L2-normalise per cell so scores are comparable across cells
col_norms <- Matrix::colSums(gsm)
col_norms[col_norms == 0] <- 1
gsm_norm  <- Matrix::t(Matrix::t(gsm) / col_norms) * 1e4  # analogous to CPM scaling

cat(sprintf("  GeneScoreMatrix: %d genes × %d cells\n", nrow(gsm_norm), ncol(gsm_norm)))

# Per-group mean score across member genes (mean of normalised scores)
module_scores <- lapply(names(gene_groups), function(grp) {
  genes <- gene_groups[[grp]]
  idx   <- which(gene_names %in% genes)
  found <- gene_names[idx]
  cat(sprintf("  %s: %d/%d genes found (%s)\n",
              grp, length(found), length(genes), paste(found, collapse = ", ")))
  if (length(idx) == 0) return(rep(0, ncol(gsm_norm)))
  if (length(idx) == 1) return(as.numeric(gsm_norm[idx, ]))
  Matrix::colMeans(gsm_norm[idx, , drop = FALSE])
})
names(module_scores) <- names(gene_groups)

# ── Extract per-cell data ─────────────────────────────────────────────────────
cell_data <- getCellColData(proj, select = "Clusters") %>% as.data.frame()
cell_data$barcode <- rownames(cell_data)
for (grp in names(module_scores)) {
  cell_data[[grp]] <- module_scores[[grp]]
}

# ── Score clusters ─────────────────────────────────────────────────────────────
cat("[9_extract_normal_barcodes] Scoring clusters by immune/stromal composite...\n")
cluster_scores <- cell_data %>%
  group_by(Clusters) %>%
  summarise(
    mean_BScore = mean(BScore, na.rm = TRUE),
    mean_TScore = mean(TScore, na.rm = TRUE),
    mean_MScore = mean(MScore, na.rm = TRUE),
    mean_FScore = mean(FScore, na.rm = TRUE),
    n_cells     = n(),
    .groups = "drop"
  ) %>%
  mutate(
    normal_composite = mean_BScore + mean_TScore + mean_MScore + mean_FScore
  ) %>%
  arrange(desc(normal_composite))

cat("\nCluster scores (descending normal composite):\n")
print(as.data.frame(cluster_scores), row.names = FALSE)

# Threshold: composite > median + 1 SD
threshold <- median(cluster_scores$normal_composite) + sd(cluster_scores$normal_composite)
normal_clusters <- cluster_scores$Clusters[cluster_scores$normal_composite > threshold]

cat(sprintf("\nThreshold: %.4f (median + 1 SD)\n", threshold))
cat(sprintf("Normal clusters: %s\n", paste(sort(normal_clusters), collapse = ", ")))
cat(sprintf("Tumor clusters:  %s\n",
            paste(sort(setdiff(cluster_scores$Clusters, normal_clusters)), collapse = ", ")))

# ── Export barcodes ───────────────────────────────────────────────────────────
normal_barcodes <- cell_data$barcode[cell_data$Clusters %in% normal_clusters]
tumor_barcodes  <- cell_data$barcode[!cell_data$Clusters %in% normal_clusters]

cat(sprintf("\nNormal spots: %d / %d total (%.1f%%)\n",
            length(normal_barcodes), nrow(cell_data),
            100 * length(normal_barcodes) / nrow(cell_data)))

out_csv <- file.path(out_dir, sprintf("%s_normal_barcodes.csv", tissue))
write.table(
  data.frame(barcode = normal_barcodes),
  file      = out_csv,
  sep       = ",",
  row.names = FALSE,
  col.names = FALSE,  # CalicoST expects no header
  quote     = FALSE
)
cat(sprintf("Saved: %s\n", out_csv))

# Also save a summary table
summary_csv <- file.path(out_dir, sprintf("%s_cluster_scores.csv", tissue))
write.csv(cluster_scores, summary_csv, row.names = FALSE, quote = FALSE)
cat(sprintf("Saved cluster scores: %s\n", summary_csv))

# ── Plot 1: Module score distributions per cluster ────────────────────────────
cat("[9_extract_normal_barcodes] Plotting module scores per cluster...\n")
plot_list <- lapply(names(gene_groups), function(sc) {
  ggplot(cell_data, aes(x = Clusters, y = .data[[sc]], fill = Clusters)) +
    geom_violin(scale = "width", alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.9) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = sc, x = NULL, y = "Module score") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
})
pdf(file.path(out_dir, sprintf("%s_module_score_by_cluster.pdf", tissue)),
    width = 10, height = 3 * length(gene_groups))
for (p in plot_list) print(p)
dev.off()

# Composite score bar plot
p_comp <- ggplot(cluster_scores,
                 aes(x = reorder(Clusters, -normal_composite),
                     y = normal_composite,
                     fill = Clusters %in% normal_clusters)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = threshold, linetype = "dashed", colour = "black") +
  annotate("text", x = Inf, y = threshold,
           label = sprintf("threshold = %.2f", threshold),
           hjust = 1.1, vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("TRUE" = "#4575b4", "FALSE" = "#d73027"),
                    labels = c("TRUE" = "Normal", "FALSE" = "Tumor")) +
  labs(title = sprintf("%s: Immune/stromal composite score per cluster", tissue),
       x = "Cluster", y = "Composite score (B+T+M+F)", fill = NULL) +
  theme_bw(base_size = 12)

pdf(file.path(out_dir, sprintf("%s_composite_score_bar.pdf", tissue)),
    width = 7, height = 4)
print(p_comp)
dev.off()
cat(sprintf("Saved: %s_module_score_by_cluster.pdf\n", tissue))
cat(sprintf("Saved: %s_composite_score_bar.pdf\n", tissue))

# ── Plot 2: UMAP colored by module scores ─────────────────────────────────────
cat("[9_extract_normal_barcodes] Plotting UMAPs...\n")

# Store normal/tumor label and module scores in cellColData for plotEmbedding
proj <- addCellColData(
  proj,
  data   = ifelse(proj$Clusters %in% normal_clusters, "normal", "tumor"),
  name   = "normal_tumor",
  cells  = getCellNames(proj),
  force  = TRUE
)
for (grp in names(gene_groups)) {
  proj <- addCellColData(
    proj,
    data   = module_scores[[grp]],
    name   = paste0("ModuleScore_", grp),
    cells  = getCellNames(proj),
    force  = TRUE
  )
}

# UMAP plots via ArchR
p_umap_clusters <- plotEmbedding(
  ArchRProj    = proj,
  embedding    = "UMAP",
  colorBy      = "cellColData",
  name         = "Clusters",
  title        = sprintf("%s: Clusters", tissue),
  plotAs       = "points",
  size         = 0.5
)
p_umap_normal <- plotEmbedding(
  ArchRProj    = proj,
  embedding    = "UMAP",
  colorBy      = "cellColData",
  name         = "normal_tumor",
  title        = sprintf("%s: Normal vs Tumor", tissue),
  plotAs       = "points",
  size         = 0.5
)
module_umaps <- lapply(names(gene_groups), function(sc) {
  plotEmbedding(
    ArchRProj = proj,
    embedding = "UMAP",
    colorBy   = "cellColData",
    name      = paste0("ModuleScore_", sc),
    title     = sprintf("%s: %s", tissue, sc),
    plotAs    = "points",
    size      = 0.3
  )
})

pdf(file.path(out_dir, sprintf("%s_module_scores_umap.pdf", tissue)),
    width = 10, height = 5)
print(p_umap_clusters)
print(p_umap_normal)
for (p in module_umaps) print(p)
dev.off()
cat(sprintf("Saved: %s_module_scores_umap.pdf\n", tissue))

# ── Plot 3: Spatial map normal vs tumor ───────────────────────────────────────
if (!is.null(spatial_path) && file.exists(spatial_path)) {
  cat("[9_extract_normal_barcodes] Plotting spatial normal/tumor map...\n")

  coords <- tryCatch({
    df <- read.csv(spatial_path, header = FALSE,
                   col.names = c("barcode", "in_tissue", "row", "col", "y", "x"))
    df <- df[df$in_tissue == 1, ]
    df$barcode_full <- paste0(
      sub("_", "#", tissue, fixed = FALSE), "#", df$barcode, "-1"
    )
    # Handle tissue name format (lowseq_489 → Lowseq_489)
    sample_prefix <- paste0(
      toupper(substr(tissue, 1, 1)), substr(tissue, 2, nchar(tissue))
    )
    df$barcode_full <- paste0(sample_prefix, "#", df$barcode, "-1")
    df
  }, error = function(e) NULL)

  if (!is.null(coords)) {
    coords$label <- ifelse(coords$barcode_full %in% normal_barcodes, "Normal", "Tumor")
    p_spatial <- ggplot(coords, aes(x = x, y = y, colour = label)) +
      geom_point(size = 0.8, alpha = 0.7) +
      scale_colour_manual(values = c("Normal" = "#4575b4", "Tumor" = "#d73027")) +
      scale_y_reverse() +
      coord_equal() +
      labs(title = sprintf("%s: Normal vs Tumor (spatial)", tissue),
           subtitle = sprintf("%d normal / %d total spots (threshold: composite > %.2f)",
                              length(normal_barcodes), nrow(cell_data), threshold),
           colour = NULL, x = "X", y = "Y") +
      theme_bw(base_size = 11)
    pdf(file.path(out_dir, sprintf("%s_normal_tumor_spatial.pdf", tissue)),
        width = 9, height = 4)
    print(p_spatial)
    dev.off()
    cat(sprintf("Saved: %s_normal_tumor_spatial.pdf\n", tissue))
  }
}

cat(sprintf(
  "\n[9_extract_normal_barcodes] Done.\n  Normal barcodes: %s\n  To use in CalicoST:\n    normalidx_file : %s\n",
  out_csv, out_csv
))
