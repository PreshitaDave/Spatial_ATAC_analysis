#!/usr/bin/env Rscript
# 00_prepare_scrna_489.R
# Finalize scRNA-seq QC for tissue 489 (WM_DR_SJ_487) and save processed Seurat object.
# Parameters follow tiss_488B_explore.Rmd: nFeature_RNA>200, percent.mt<20, dims=1:15.
# Outputs:
#   seurat_scrna_489_processed.rds  — full object (all cells)
#   seurat_scrna_489_balanced.rds   — downsampled to equal cells per cell type

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
})

set.seed(42)

# ── Paths ─────────────────────────────────────────────────────────────────────
# Note: sample 489 data lives in WM_DR_SJ_487 directory
DATA_10X <- "/projectnb/paxlab/DATA/DriesSpatial/scRNAseq/WM_DR_SJ_487/outs/filtered_feature_bc_matrix"
OUT_DIR  <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration/489"
PLOT_DIR <- file.path(OUT_DIR, "plots", "label_transfer")
OBJ_DIR  <- file.path(OUT_DIR, "objects")
TAB_DIR  <- file.path(OUT_DIR, "tables")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJ_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR,  recursive = TRUE, showWarnings = FALSE)

# ── Load & QC ─────────────────────────────────────────────────────────────────
cat("Loading 10x data from:", DATA_10X, "\n")
counts <- Read10X(DATA_10X)
obj    <- CreateSeuratObject(counts, project = "tiss_489", min.cells = 3, min.features = 200)
cat("Cells before QC:", ncol(obj), "\n")

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

p_qc <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3, pt.size = 0) &
        theme(plot.title = element_text(size = 10))
ggsave(file.path(PLOT_DIR, "qc_violin_prefilter.pdf"), p_qc, width = 12, height = 4)

# Rmd params: nFeature_RNA > 200, percent.mt < 20
obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 20)
cat("Cells after QC:", ncol(obj), "\n")

# ── Standard Seurat pipeline (Rmd params: dims=1:15, nfeatures=2000) ──────────
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = rownames(obj))
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

p_elbow <- ElbowPlot(obj, ndims = 30) + ggtitle("489 scRNA — PCA elbow plot")
ggsave(file.path(PLOT_DIR, "pca_elbow.pdf"), p_elbow, width = 6, height = 4)

obj <- FindNeighbors(obj, dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:15)

p_umap <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE) +
          ggtitle("489 scRNA-seq clusters")
ggsave(file.path(PLOT_DIR, "umap_clusters.pdf"), p_umap, width = 7, height = 6)

# ── Marker genes ──────────────────────────────────────────────────────────────
cat("Finding cluster markers...\n")
markers <- FindAllMarkers(obj, only.pos = TRUE, logfc.threshold = 0.5,
                          min.pct = 0.25, test.use = "wilcox")
write.csv(markers, file.path(TAB_DIR, "scrna_489_cluster_markers.csv"), row.names = FALSE)

# ── Cell type annotation ───────────────────────────────────────────────────────
canonical <- list(
  Tumor       = c("EPCAM", "KRT8", "KRT18", "ESR1", "PGR", "ERBB2", "KRT14", "KRT17"),
  T_cell      = c("CD3E", "CD3D", "CD3G", "CD8A", "CD4"),
  NK_cell     = c("GNLY", "NKG7", "KLRD1"),
  B_cell      = c("MS4A1", "CD79A", "CD79B"),
  Myeloid     = c("CD14", "LYZ", "FCGR3A", "CST3", "CD68"),
  Fibroblast  = c("COL1A1", "COL1A2", "DCN", "FAP", "VIM"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5")
)

p_markers <- FeaturePlot(obj, features = unlist(lapply(canonical, head, 2)),
                          ncol = 4, min.cutoff = "q05") &
             theme(axis.title = element_blank(), axis.text = element_blank())
ggsave(file.path(PLOT_DIR, "canonical_markers.pdf"), p_markers, width = 16, height = 14)

p_dot <- DotPlot(obj, features = lapply(canonical, head, 3),
                  cols = c("lightgrey", "red")) +
         RotatedAxis() + ggtitle("489: canonical markers by cluster")
ggsave(file.path(PLOT_DIR, "dotplot_canonical_markers.pdf"), p_dot, width = 14, height = 6)

cat("Scoring clusters by canonical markers...\n")
obj <- AddModuleScore(obj, features = canonical, name = "celltype_score_")

n_types    <- length(canonical)
type_names <- names(canonical)
score_cols <- paste0("celltype_score_", seq_len(n_types))

cluster_scores <- obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols), mean), .groups = "drop")

cluster_scores$predicted_cell_type <- type_names[
  apply(cluster_scores[, score_cols], 1, which.max)
]

cat("Cluster → cell type assignments:\n")
print(cluster_scores[, c("seurat_clusters", "predicted_cell_type")])

cluster_map   <- setNames(cluster_scores$predicted_cell_type, cluster_scores$seurat_clusters)
obj$cell_type <- unname(cluster_map[as.character(obj$seurat_clusters)])

p_celltypes <- DimPlot(obj, group.by = "cell_type", label = TRUE, repel = TRUE) +
               ggtitle("489 scRNA-seq — annotated cell types")
ggsave(file.path(PLOT_DIR, "umap_cell_types.pdf"), p_celltypes, width = 8, height = 6)

cat("\nCell type composition (489 — full):\n")
print(table(obj$cell_type))
write.csv(data.frame(table(obj$cell_type)),
          file.path(TAB_DIR, "scrna_489_celltype_counts.csv"), row.names = FALSE)
write.csv(obj@meta.data, file.path(TAB_DIR, "scrna_489_metadata.csv"))

# ── Save full object ───────────────────────────────────────────────────────────
out_full <- file.path(OBJ_DIR, "seurat_scrna_489_processed.rds")
saveRDS(obj, out_full)
cat("Saved full object:", out_full, "\n")

# ── Balanced downsampling ──────────────────────────────────────────────────────
# Downsample each cell type to the size of the smallest class to remove imbalance
# bias in label transfer.
cat("\nCreating balanced reference (equal cells per cell type)...\n")
ct_counts <- table(obj$cell_type)
min_n     <- min(ct_counts)
cat("Class sizes before balancing:\n"); print(ct_counts)
cat("Downsampling each class to:", min_n, "cells\n")

balanced_cells <- obj@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  group_by(cell_type) %>%
  slice_sample(n = min_n) %>%
  pull(barcode)

obj_balanced <- subset(obj, cells = balanced_cells)
cat("Balanced object:", ncol(obj_balanced), "cells\n")
cat("Balanced cell type distribution:\n")
print(table(obj_balanced$cell_type))

# Re-run PCA + UMAP on the balanced subset so the reference embedding is
# computed from balanced data (not carried over from the imbalanced object).
obj_balanced <- FindVariableFeatures(obj_balanced, selection.method = "vst", nfeatures = 2000)
obj_balanced <- ScaleData(obj_balanced, features = rownames(obj_balanced))
obj_balanced <- RunPCA(obj_balanced, npcs = 30, verbose = FALSE)
obj_balanced <- FindNeighbors(obj_balanced, dims = 1:15)
obj_balanced <- FindClusters(obj_balanced, resolution = 0.5)
obj_balanced <- RunUMAP(obj_balanced, dims = 1:15)

p_bal <- DimPlot(obj_balanced, group.by = "cell_type", label = TRUE, repel = TRUE) +
         ggtitle(paste0("489 scRNA — balanced reference (", min_n, " cells/type)"))
ggsave(file.path(PLOT_DIR, "umap_cell_types_balanced.pdf"), p_bal, width = 8, height = 6)

out_bal <- file.path(OBJ_DIR, "seurat_scrna_489_balanced.rds")
saveRDS(obj_balanced, out_bal)
cat("Saved balanced object:", out_bal, "\n")

write.csv(data.frame(table(obj_balanced$cell_type)),
          file.path(TAB_DIR, "scrna_489_celltype_counts_balanced.csv"), row.names = FALSE)

# ── Export matrices ────────────────────────────────────────────────────────────
cat("Exporting normalized expression matrices...\n")
for (tag in c("processed", "balanced")) {
  o  <- if (tag == "processed") obj else obj_balanced
  nm <- GetAssayData(o, layer = "data")
  writeMM(nm, file.path(OBJ_DIR, paste0("scrna_489_norm_matrix_", tag, ".mtx")))
  write.csv(rownames(nm), file.path(OBJ_DIR, paste0("scrna_489_genes_", tag, ".csv")),    row.names = FALSE, col.names = FALSE)
  write.csv(colnames(nm), file.path(OBJ_DIR, paste0("scrna_489_barcodes_", tag, ".csv")), row.names = FALSE, col.names = FALSE)
  write.csv(o@meta.data[, c("seurat_clusters", "cell_type", "nFeature_RNA", "nCount_RNA", "percent.mt")],
            file.path(OBJ_DIR, paste0("scrna_489_cell_metadata_", tag, ".csv")))
}

cat("\nDone. Output directory:", OUT_DIR, "\n")
