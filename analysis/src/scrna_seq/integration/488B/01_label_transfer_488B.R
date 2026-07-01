#!/usr/bin/env Rscript
# 01_label_transfer_488B.R
# CCA-based label transfer from scRNA → ATAC (tissue 488B / deepseq).
#
# Usage: Rscript 01_label_transfer_488B.R [full|balanced]
#   full     — use seurat_scrna_488B_processed.rds  (default, all cells)
#   balanced — use seurat_scrna_488B_balanced.rds   (equal cells per cell type)
#
# Prerequisites: 00_prepare_scrna_488B.R must have been run.

suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(ComplexHeatmap)
  library(circlize)
})

set.seed(42)
addArchRThreads(threads = 8)
addArchRGenome("hg38")

# ── Mode ──────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args) > 0 && args[1] == "balanced") "balanced" else "full"
cat("Running label transfer in MODE:", MODE, "\n")

# ── Paths ─────────────────────────────────────────────────────────────────────
ARCHR_DIR  <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/01_outputs/archR_objects/5000bp_non_binarize/deepseq_488B"
OBJ_BASE   <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration/488B/objects"
OUT_DIR    <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration/488B"

scrna_tag  <- if (MODE == "balanced") "balanced" else "processed"
SCRNA_RDS  <- file.path(OBJ_BASE, paste0("seurat_scrna_488B_", scrna_tag, ".rds"))
ARCHR_OUT  <- file.path(OBJ_BASE, paste0("archr_488B_with_integration_", MODE))
PLOT_DIR   <- file.path(OUT_DIR, "plots", paste0("label_transfer_", MODE))
TAB_DIR    <- file.path(OUT_DIR, "tables")

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR,  recursive = TRUE, showWarnings = FALSE)

if (!file.exists(SCRNA_RDS)) stop("Run 00_prepare_scrna_488B.R first: ", SCRNA_RDS)

cat("Loading ArchR project:", ARCHR_DIR, "\n")
proj  <- loadArchRProject(ARCHR_DIR, showLogo = FALSE)
cat("ATAC cells:", nCells(proj), "\n")

cat("Loading scRNA Seurat object:", SCRNA_RDS, "\n")
seRNA <- readRDS(SCRNA_RDS)
cat("scRNA cells:", ncol(seRNA), "  Cell types:", paste(sort(unique(seRNA$cell_type)), collapse=", "), "\n")
cat("Cell type distribution:\n"); print(table(seRNA$cell_type))

# ── Build Seurat ATAC object from GeneScoreMatrix ─────────────────────────────
cat("Extracting GeneScoreMatrix from ArchR...\n")
gs_se  <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
gs_mat <- assay(gs_se)
if (!is.null(rowData(gs_se)$name)) {
  rownames(gs_mat) <- make.unique(as.character(rowData(gs_se)$name))
}
cat("GeneScoreMatrix:", nrow(gs_mat), "genes ×", ncol(gs_mat), "cells\n")

atac_obj <- CreateSeuratObject(counts = gs_mat, assay = "ACTIVITY", project = "atac_488B")

avail_rd <- names(proj@reducedDims)
lsi_name <- avail_rd[grep("LSI|lsi", avail_rd, ignore.case = TRUE)][1]
lsi_mat  <- as.matrix(proj@reducedDims[[lsi_name]][["matSVD"]])
colnames(lsi_mat) <- paste0("LSI_", seq_len(ncol(lsi_mat)))

atac_bare <- sub(".*#", "", colnames(atac_obj))
lsi_bare  <- sub(".*#", "", rownames(lsi_mat))
lsi_mat   <- lsi_mat[match(atac_bare, lsi_bare), , drop = FALSE]
rownames(lsi_mat) <- colnames(atac_obj)

atac_obj[["lsi"]] <- CreateDimReducObject(embeddings = lsi_mat, key = "LSI_", assay = "ACTIVITY")

DefaultAssay(atac_obj) <- "ACTIVITY"
atac_obj <- NormalizeData(atac_obj)
atac_obj <- FindVariableFeatures(atac_obj, nfeatures = 2000)
atac_obj <- ScaleData(atac_obj, features = rownames(atac_obj))

# ── FindTransferAnchors ───────────────────────────────────────────────────────
shared_genes <- intersect(rownames(atac_obj), rownames(seRNA))
cat("Shared genes for anchor finding:", length(shared_genes), "\n")

if (!"pca" %in% names(seRNA@reductions)) {
  seRNA <- RunPCA(seRNA, npcs = 30, verbose = FALSE)
}

cat("Finding transfer anchors (CCA)...\n")
anchors <- FindTransferAnchors(
  reference = seRNA,
  query     = atac_obj,
  features  = shared_genes,
  reduction = "cca",
  dims      = 1:30,
  verbose   = TRUE
)
cat("Found", nrow(anchors@anchors), "anchors.\n")

# ── Transfer labels ───────────────────────────────────────────────────────────
cat("Transferring cell type labels...\n")
pred_labels <- TransferData(
  anchorset        = anchors,
  refdata          = seRNA$cell_type,
  weight.reduction = atac_obj[["lsi"]],
  dims             = 1:20
)
atac_obj <- AddMetaData(atac_obj, pred_labels)

cat("Cell type distribution (predicted):\n")
print(table(atac_obj$predicted.id))
cat("Mean score:", round(mean(atac_obj$prediction.score.max, na.rm = TRUE), 3), "\n")
cat("% score >0.5:", round(100 * mean(atac_obj$prediction.score.max > 0.5, na.rm = TRUE), 1), "%\n")

# ── Prediction score by cell type (task 2) ────────────────────────────────────
cat("\nPrediction score breakdown by predicted cell type:\n")
score_by_type <- atac_obj@meta.data %>%
  group_by(predicted.id) %>%
  summarise(
    n_cells    = n(),
    mean_score = round(mean(prediction.score.max, na.rm = TRUE), 3),
    pct_gt0.5  = round(100 * mean(prediction.score.max > 0.5, na.rm = TRUE), 1),
    pct_gt0.7  = round(100 * mean(prediction.score.max > 0.7, na.rm = TRUE), 1),
    .groups    = "drop"
  ) %>%
  arrange(desc(n_cells))
print(score_by_type)
write.csv(score_by_type,
          file.path(TAB_DIR, paste0("prediction_score_by_celltype_488B_", MODE, ".csv")),
          row.names = FALSE)

# ── Transfer imputed gene expression ──────────────────────────────────────────
cat("Transferring imputed gene expression...\n")
pred_expr   <- TransferData(
  anchorset        = anchors,
  refdata          = GetAssayData(seRNA, layer = "data"),
  weight.reduction = atac_obj[["lsi"]],
  dims             = 1:20
)
imputed_mat <- GetAssayData(pred_expr, layer = "data")
cat("Imputed expression matrix:", nrow(imputed_mat), "×", ncol(imputed_mat), "\n")

# ── Write predictions back to ArchR ───────────────────────────────────────────
cat("Writing predictions into ArchR project...\n")
pred_df <- data.frame(
  predictedGroup_co = atac_obj$predicted.id,
  predictedScore_co = atac_obj$prediction.score.max,
  row.names         = colnames(atac_obj)
)
proj_bare   <- sub(".*#", "", proj$cellNames)
pred_aligned <- pred_df[match(proj_bare, sub(".*#", "", rownames(pred_df))), ]
proj$predictedGroup_co <- pred_aligned$predictedGroup_co
proj$predictedScore_co <- pred_aligned$predictedScore_co

# ── UMAP + plots ──────────────────────────────────────────────────────────────
proj <- addUMAP(proj, reducedDims = "IterativeLSI", name = "UMAP_LSI", force = TRUE)

p1 <- plotEmbedding(proj, colorBy = "cellColData", name = "predictedGroup_co",
                    embedding = "UMAP_LSI", size = 0.5, plotAs = "points") +
      ggtitle(paste0("488B ATAC — predicted cell type (", MODE, " ref)"))
p2 <- plotEmbedding(proj, colorBy = "cellColData", name = "predictedScore_co",
                    embedding = "UMAP_LSI", size = 0.5, plotAs = "points") +
      ggtitle("488B ATAC — prediction confidence") + scale_color_viridis_c()
p3 <- plotEmbedding(proj, colorBy = "cellColData", name = "Clusters",
                    embedding = "UMAP_LSI", size = 0.5, plotAs = "points") +
      ggtitle("488B ATAC — original clusters")

pdf(file.path(PLOT_DIR, paste0("umap_label_transfer_488B_", MODE, ".pdf")), width = 18, height = 6)
print(p1 | p2 | p3)
dev.off()

# Confusion heatmap
meta      <- as.data.frame(getCellColData(proj, select = c("Clusters", "predictedGroup_co", "predictedScore_co")))
confusion     <- table(meta$Clusters, meta$predictedGroup_co)
confusion_pct <- sweep(confusion, 1, rowSums(confusion), "/")

pdf(file.path(PLOT_DIR, paste0("confusion_heatmap_488B_", MODE, ".pdf")), width = 10, height = 8)
Heatmap(confusion_pct,
        name = "Fraction",
        col  = colorRamp2(c(0, 0.5, 1), c("white", "orange", "darkred")),
        cluster_rows = TRUE, cluster_columns = TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (confusion_pct[i, j] > 0.15)
            grid.text(sprintf("%.0f%%", confusion_pct[i, j] * 100), x, y,
                      gp = gpar(fontsize = 8))
        },
        row_title = "ATAC cluster", column_title = "Predicted cell type",
        column_title_side = "bottom")
dev.off()

# Score distribution by cell type
p_score <- ggplot(meta, aes(x = predictedScore_co, fill = predictedGroup_co)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  facet_wrap(~predictedGroup_co, scales = "free_y") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0.7, linetype = "dashed", color = "blue") +
  labs(x = "Prediction score", y = "N cells",
       title = paste0("488B: label transfer confidence by cell type (", MODE, " ref)"),
       caption = "Red=0.5, Blue=0.7 thresholds") +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(PLOT_DIR, paste0("prediction_score_distribution_488B_", MODE, ".pdf")),
       p_score, width = 14, height = 8)

# ── Save ──────────────────────────────────────────────────────────────────────
saveArchRProject(proj, outputDirectory = ARCHR_OUT, overwrite = TRUE, load = FALSE)

meta_full <- as.data.frame(getCellColData(proj))
write.csv(meta_full, file.path(TAB_DIR, paste0("label_transfer_results_488B_", MODE, ".csv")))

writeMM(imputed_mat, file.path(OBJ_BASE, paste0("archr_488B_GeneIntegrationMatrix_", MODE, ".mtx")))
write.csv(data.frame(gene = rownames(imputed_mat)),
          file.path(OBJ_BASE, paste0("archr_488B_gene_names_", MODE, ".csv")), row.names = FALSE)
write.csv(data.frame(cell = colnames(imputed_mat)),
          file.path(OBJ_BASE, paste0("archr_488B_cell_names_", MODE, ".csv")), row.names = FALSE)

cat("\nLabel transfer (", MODE, ") complete for 488B.\n")
cat("ArchR project saved to:", ARCHR_OUT, "\n")
