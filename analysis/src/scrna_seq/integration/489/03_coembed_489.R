#!/usr/bin/env Rscript
# 03_coembed_489.R
# Option 3: Co-embedding for tissue 489 (lowseq_489).
#   3a: ArchR combined LSI+RNA dims → joint UMAP (uses predictedGroup_co from ArchR)
#   3b: Seurat bridge integration (GeneScoreMatrix → FindTransferAnchors)
#
# Usage: Rscript 03_coembed_489.R [full|balanced]
#   full     — uses archr_489_with_integration + seurat_scrna_489_processed.rds
#   balanced — uses archr_489_with_integration_balanced + seurat_scrna_489_balanced.rds
#
# Prerequisites: 01_label_transfer_489.R must have been run for the chosen mode.

suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

set.seed(42)
addArchRThreads(threads = 8)
addArchRGenome("hg38")

# ── Mode ──────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args) > 0 && args[1] == "balanced") "balanced" else "full"
cat("Running co-embedding in MODE:", MODE, "\n")

# ── Paths ─────────────────────────────────────────────────────────────────────
OBJ_DIR   <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration/489/objects"
OUT_DIR   <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration/489"

scrna_tag <- if (MODE == "balanced") "balanced" else "processed"
ARCHR_INT <- file.path(OBJ_DIR, paste0("archr_489_with_integration_", MODE))
SCRNA_RDS <- file.path(OBJ_DIR, paste0("seurat_scrna_489_", scrna_tag, ".rds"))
PLOT_DIR  <- file.path(OUT_DIR, "plots", paste0("coembed_", MODE))
TAB_DIR   <- file.path(OUT_DIR, "tables")

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR,  recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(ARCHR_INT)) stop("Run 01_label_transfer_489.R (", MODE, ") first: ", ARCHR_INT)

cat("Loading ArchR project:", ARCHR_INT, "\n")
proj  <- loadArchRProject(ARCHR_INT, showLogo = FALSE)
seRNA <- readRDS(SCRNA_RDS)
cat("ATAC cells:", nCells(proj), "| scRNA cells:", ncol(seRNA), "\n")

# ──────────────────────────────────────────────────────────────────────────────
# Approach 3a: ArchR combined LSI+RNA dims → joint UMAP
# ──────────────────────────────────────────────────────────────────────────────
cat("\n=== Approach 3a: ArchR combined LSI + RNA dimensions ===\n")

avail_dims <- getAvailableMatrices(proj)
cat("Available matrices:", paste(avail_dims, collapse=", "), "\n")
cat("Available reduced dims:", paste(names(proj@reducedDims), collapse=", "), "\n")

if ("LSI_ATAC_RNA" %in% names(proj@reducedDims) || "LSI_Combined" %in% names(proj@reducedDims)) {
  combined_dim_name <- if ("LSI_ATAC_RNA" %in% names(proj@reducedDims)) "LSI_ATAC_RNA" else "LSI_Combined"
  cat("Using existing combined dims:", combined_dim_name, "\n")
} else {
  cat("Creating combined LSI + GeneIntegration reduced dims...\n")
  proj <- addCombinedDims(
    ArchRProj   = proj,
    reducedDims = c("IterativeLSI"),
    name        = "LSI_ATAC_RNA"
  )
  combined_dim_name <- "LSI_ATAC_RNA"
}

proj <- addUMAP(proj, reducedDims = combined_dim_name,
                name = "UMAP_Combined", nNeighbors = 30, force = TRUE)
proj <- addClusters(proj, reducedDims = combined_dim_name,
                    name = "Clusters_Combined", resolution = 0.5, force = TRUE)

p_combined1 <- plotEmbedding(proj, colorBy = "cellColData", name = "predictedGroup_co",
                              embedding = "UMAP_Combined", size = 0.5, plotAs = "points") +
               ggtitle(paste0("489 joint UMAP — predicted cell type (", MODE, ")"))
p_combined2 <- plotEmbedding(proj, colorBy = "cellColData", name = "Clusters",
                              embedding = "UMAP_Combined", size = 0.5, plotAs = "points") +
               ggtitle("489 joint UMAP — original ATAC clusters")
p_combined3 <- plotEmbedding(proj, colorBy = "cellColData", name = "Clusters_Combined",
                              embedding = "UMAP_Combined", size = 0.5, plotAs = "points") +
               ggtitle("489 joint UMAP — combined clusters")

pdf(file.path(PLOT_DIR, paste0("3a_archr_joint_umap_489_", MODE, ".pdf")), width = 18, height = 6)
print(p_combined1 | p_combined2 | p_combined3)
dev.off()
cat("Saved: 3a_archr_joint_umap_489_", MODE, ".pdf\n")

# ──────────────────────────────────────────────────────────────────────────────
# Approach 3b: Seurat bridge integration via FindTransferAnchors
# ──────────────────────────────────────────────────────────────────────────────
cat("\n=== Approach 3b: Seurat bridge integration (gene activity → RNA) ===\n")
cat("Using scRNA reference:", SCRNA_RDS, "\n")

gs     <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
gs_mat <- assay(gs)
if (!is.null(rowData(gs)$name))
  rownames(gs_mat) <- make.unique(as.character(rowData(gs)$name))
cat("GeneScoreMatrix:", nrow(gs_mat), "genes ×", ncol(gs_mat), "cells\n")

atac_seurat <- CreateSeuratObject(counts = gs_mat, assay = "ACTIVITY", project = "atac_489")

avail_rd <- names(proj@reducedDims)
lsi_name <- avail_rd[grep("LSI|lsi", avail_rd, ignore.case=TRUE)][1]
lsi_emb  <- as.matrix(proj@reducedDims[[lsi_name]][["matSVD"]])
lsi_mat  <- lsi_emb
colnames(lsi_mat) <- paste0("LSI_", seq_len(ncol(lsi_mat)))

common_cells <- intersect(colnames(atac_seurat), rownames(lsi_mat))
atac_seurat  <- atac_seurat[, common_cells]
lsi_mat      <- lsi_mat[common_cells, ]

atac_seurat[["lsi"]] <- CreateDimReducObject(
  embeddings = lsi_mat, key = "LSI_", assay = "ACTIVITY"
)

DefaultAssay(atac_seurat) <- "ACTIVITY"
atac_seurat <- NormalizeData(atac_seurat)
atac_seurat <- FindVariableFeatures(atac_seurat)
atac_seurat <- ScaleData(atac_seurat, features = rownames(atac_seurat))

shared_genes <- intersect(rownames(atac_seurat), rownames(seRNA))
cat("Shared genes for anchor finding:", length(shared_genes), "\n")

if (!"pca" %in% names(seRNA@reductions))
  seRNA <- RunPCA(seRNA, npcs = 30, verbose = FALSE)

cat("Finding transfer anchors (CCA)...\n")
transfer_anchors <- FindTransferAnchors(
  reference           = seRNA,
  query               = atac_seurat,
  features            = shared_genes,
  reference.reduction = "pca",
  reduction           = "pcaproject",
  dims                = 1:30,
  verbose             = TRUE
)

n_anchors <- nrow(transfer_anchors@anchors)
cat("Anchors found:", n_anchors, "\n")
kw <- min(20, n_anchors - 1)
cat("Using k.weight =", kw, "\n")
predictions <- TransferData(
  anchorset        = transfer_anchors,
  refdata          = seRNA$cell_type,
  weight.reduction = atac_seurat[["lsi"]],
  dims             = 1:20,
  k.weight         = kw
)
atac_seurat <- AddMetaData(atac_seurat, metadata = predictions)

cat("Seurat transfer cell type distribution:\n")
print(table(atac_seurat$predicted.id))
cat("Mean prediction score:", round(mean(atac_seurat$prediction.score.max, na.rm=TRUE), 3), "\n")

atac_seurat <- RunUMAP(atac_seurat, reduction = "lsi", dims = 1:20,
                       reduction.name = "umap_lsi")

p_bridge <- DimPlot(atac_seurat, reduction = "umap_lsi",
                    group.by = "predicted.id", label = TRUE, repel = TRUE) +
            ggtitle(paste0("489 Seurat bridge: predicted cell type (", MODE, " ref)")) +
            theme(legend.position = "right")
ggsave(file.path(PLOT_DIR, paste0("3b_seurat_bridge_umap_489_", MODE, ".pdf")),
       p_bridge, width = 9, height = 7)

# ArchR vs Seurat comparison
meta_archr  <- as.data.frame(getCellColData(proj, select = c("predictedGroup_co", "predictedScore_co")))
meta_archr$cell <- rownames(meta_archr)
meta_seurat <- data.frame(
  cell         = colnames(atac_seurat),
  seurat_pred  = atac_seurat$predicted.id,
  seurat_score = atac_seurat$prediction.score.max
)
comparison <- merge(meta_archr, meta_seurat, by = "cell", all = FALSE)
comparison$methods_agree <- comparison$predictedGroup_co == comparison$seurat_pred
cat("ArchR vs Seurat agreement:", round(100*mean(comparison$methods_agree, na.rm=TRUE), 1), "%\n")

write.csv(comparison,
          file.path(TAB_DIR, paste0("coembed_method_comparison_489_", MODE, ".csv")),
          row.names = FALSE)

agree_table <- table(comparison$predictedGroup_co, comparison$seurat_pred)
p_agree <- ggplot(as.data.frame(agree_table), aes(x=Var2, y=Var1, fill=Freq)) +
  geom_tile() +
  geom_text(aes(label=Freq), size=3) +
  scale_fill_gradient(low="white", high="steelblue") +
  labs(x="Seurat prediction", y="ArchR prediction",
       title=paste0("489: ArchR vs Seurat agreement (", MODE, " ref)"),
       fill="N cells") +
  theme_bw() + theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(file.path(PLOT_DIR, paste0("3b_archr_vs_seurat_agreement_489_", MODE, ".pdf")),
       p_agree, width=9, height=8)

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(atac_seurat, file.path(OBJ_DIR, paste0("seurat_atac_bridge_489_", MODE, ".rds")))
saveArchRProject(proj, outputDirectory = ARCHR_INT, overwrite = TRUE, load = FALSE)
cat("Saved Seurat bridge object and updated ArchR project.\n")
cat("\nOption 3 co-embedding (", MODE, ") complete for 489.\n")
