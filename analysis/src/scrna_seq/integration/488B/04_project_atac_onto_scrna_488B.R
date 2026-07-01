#!/usr/bin/env Rscript
# 04_project_atac_onto_scrna_488B.R
# Project ATAC-seq cells (via GeneScoreMatrix) onto scRNA-seq UMAP using
# Seurat MapQuery. Produces overlay plots for both full and balanced references.
#
# Usage: Rscript 04_project_atac_onto_scrna_488B.R [full|balanced|both]
# Prerequisites: 00_prepare_scrna_488B.R and 01_label_transfer_488B.R must have run.

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

args  <- commandArgs(trailingOnly = TRUE)
modes <- if (length(args) > 0 && args[1] %in% c("full","balanced")) args[1] else c("full","balanced")

# ── Paths ─────────────────────────────────────────────────────────────────────
ARCHR_DIR <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/01_outputs/archR_objects/5000bp_non_binarize/deepseq_488B"
OBJ_BASE  <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration/488B/objects"
OUT_DIR   <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration/488B"
PLOT_DIR  <- file.path(OUT_DIR, "plots", "atac_projection")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Load ATAC GeneScoreMatrix (shared across modes) ───────────────────────────
cat("Loading ArchR project:", ARCHR_DIR, "\n")
proj <- loadArchRProject(ARCHR_DIR, showLogo = FALSE)
cat("ATAC cells:", nCells(proj), "\n")

cat("Extracting GeneScoreMatrix...\n")
gs_se  <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
gs_mat <- assay(gs_se)
if (!is.null(rowData(gs_se)$name))
  rownames(gs_mat) <- make.unique(as.character(rowData(gs_se)$name))

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

# ── Cell type color palette (consistent across plots) ─────────────────────────
ct_colors <- c(
  Tumor       = "#E41A1C",
  T_cell      = "#377EB8",
  NK_cell     = "#FF7F00",
  B_cell      = "#4DAF4A",
  Myeloid     = "#984EA3",
  Fibroblast  = "#A65628",
  Endothelial = "#F781BF"
)

# ── Run per mode ──────────────────────────────────────────────────────────────
for (MODE in modes) {
  cat("\n========== MODE:", MODE, "==========\n")

  scrna_tag <- if (MODE == "balanced") "balanced" else "processed"
  SCRNA_RDS <- file.path(OBJ_BASE, paste0("seurat_scrna_488B_", scrna_tag, ".rds"))
  if (!file.exists(SCRNA_RDS)) { cat("SKIP: missing", SCRNA_RDS, "\n"); next }

  seRNA <- readRDS(SCRNA_RDS)
  cat("scRNA cells:", ncol(seRNA), "\n")
  cat("Cell type distribution:\n"); print(table(seRNA$cell_type))

  if (!"pca" %in% names(seRNA@reductions))
    seRNA <- RunPCA(seRNA, npcs = 30, verbose = FALSE)
  # return.model = TRUE stores the uwot model so MapQuery can project into it
  seRNA <- RunUMAP(seRNA, dims = 1:15, verbose = FALSE, return.model = TRUE)

  # ── Find anchors ────────────────────────────────────────────────────────────
  shared_genes <- intersect(rownames(atac_obj), rownames(seRNA))
  cat("Shared genes:", length(shared_genes), "\n")

  cat("Finding transfer anchors...\n")
  anchors <- FindTransferAnchors(
    reference = seRNA,
    query     = atac_obj,
    features  = shared_genes,
    reduction = "cca",
    dims      = 1:30,
    verbose   = TRUE
  )
  cat("Anchors found:", nrow(anchors@anchors), "\n")

  # ── MapQuery: project ATAC cells into scRNA UMAP space ─────────────────────
  cat("Projecting ATAC onto scRNA UMAP via MapQuery...\n")
  atac_mapped <- MapQuery(
    anchorset         = anchors,
    query             = atac_obj,
    reference         = seRNA,
    refdata           = list(predicted_celltype = "cell_type"),
    reference.reduction = "pca",
    reduction.model   = "umap"
  )
  cat("Projected ATAC cells:", ncol(atac_mapped), "\n")
  cat("Predicted cell type distribution (projected):\n")
  print(table(atac_mapped$predicted.predicted_celltype))
  cat("Mean prediction score:", round(mean(atac_mapped$predicted.predicted_celltype.score, na.rm=TRUE), 3), "\n")

  # ── Score breakdown by cell type ────────────────────────────────────────────
  score_df <- data.frame(
    cell_type = atac_mapped$predicted.predicted_celltype,
    score     = atac_mapped$predicted.predicted_celltype.score
  ) %>%
    group_by(cell_type) %>%
    summarise(n=n(), mean_score=round(mean(score,na.rm=TRUE),3),
              pct_gt0.5=round(100*mean(score>0.5,na.rm=TRUE),1),
              pct_gt0.7=round(100*mean(score>0.7,na.rm=TRUE),1), .groups="drop") %>%
    arrange(desc(n))
  cat("\nScore breakdown by predicted type:\n"); print(score_df)
  write.csv(score_df,
    file.path(OUT_DIR, "tables", paste0("mapquery_score_by_celltype_488B_", MODE, ".csv")),
    row.names = FALSE)

  # ── Plots ───────────────────────────────────────────────────────────────────
  # 1. scRNA UMAP colored by cell type
  p_scrna <- DimPlot(seRNA, reduction = "umap", group.by = "cell_type",
                     cols = ct_colors[names(ct_colors) %in% unique(seRNA$cell_type)],
                     pt.size = 0.3, label = TRUE, repel = TRUE) +
             ggtitle(paste0("488B scRNA (", scrna_tag, ")\n", ncol(seRNA), " cells")) +
             theme(legend.position = "right") + NoAxes()

  # 2. ATAC projected onto scRNA UMAP, colored by predicted type
  p_atac_proj <- DimPlot(atac_mapped, reduction = "ref.umap",
                          group.by = "predicted.predicted_celltype",
                          cols = ct_colors[names(ct_colors) %in% unique(atac_mapped$predicted.predicted_celltype)],
                          pt.size = 0.3, label = TRUE, repel = TRUE) +
                 ggtitle(paste0("488B ATAC projected onto scRNA UMAP\n",
                                ncol(atac_mapped), " ATAC cells")) +
                 theme(legend.position = "right") + NoAxes()

  # 3. Overlay: scRNA (grey) + ATAC (colored)
  scrna_coords <- as.data.frame(Embeddings(seRNA, "umap"))
  colnames(scrna_coords) <- c("UMAP_1","UMAP_2")
  scrna_coords$cell_type <- seRNA$cell_type
  scrna_coords$modality  <- "scRNA"

  atac_coords  <- as.data.frame(Embeddings(atac_mapped, "ref.umap"))
  colnames(atac_coords) <- c("UMAP_1","UMAP_2")
  atac_coords$cell_type <- atac_mapped$predicted.predicted_celltype
  atac_coords$modality  <- "ATAC"

  p_overlay <- ggplot() +
    geom_point(data = scrna_coords, aes(x=UMAP_1, y=UMAP_2),
               color = "grey80", size = 0.4, alpha = 0.5) +
    geom_point(data = atac_coords, aes(x=UMAP_1, y=UMAP_2, color=cell_type),
               size = 0.5, alpha = 0.7) +
    scale_color_manual(values = ct_colors, name = "Predicted type") +
    ggtitle(paste0("488B: ATAC projected onto scRNA UMAP (", MODE, " ref)\ngrey=scRNA, colored=ATAC")) +
    theme_bw() + theme(axis.title=element_blank(), axis.text=element_blank(),
                       axis.ticks=element_blank(), panel.grid=element_blank()) +
    guides(color = guide_legend(override.aes = list(size=3, alpha=1)))

  # 4. Side-by-side score map: prediction confidence on projected UMAP
  atac_coords$score <- atac_mapped$predicted.predicted_celltype.score
  p_score <- ggplot(atac_coords, aes(x=UMAP_1, y=UMAP_2, color=score)) +
    geom_point(size=0.4, alpha=0.7) +
    scale_color_viridis_c(option="plasma", name="Pred. score") +
    ggtitle(paste0("488B ATAC: prediction confidence on scRNA UMAP (", MODE, ")")) +
    theme_bw() + theme(axis.title=element_blank(), axis.text=element_blank(),
                       axis.ticks=element_blank(), panel.grid=element_blank())

  # Save panels
  pdf(file.path(PLOT_DIR, paste0("scrna_umap_488B_", MODE, ".pdf")), width=8, height=6)
  print(p_scrna); dev.off()

  pdf(file.path(PLOT_DIR, paste0("atac_projected_488B_", MODE, ".pdf")), width=8, height=6)
  print(p_atac_proj); dev.off()

  pdf(file.path(PLOT_DIR, paste0("overlay_atac_on_scrna_488B_", MODE, ".pdf")), width=8, height=6)
  print(p_overlay); dev.off()

  pdf(file.path(PLOT_DIR, paste0("projection_score_map_488B_", MODE, ".pdf")), width=8, height=6)
  print(p_score); dev.off()

  # Combined 4-panel summary
  pdf(file.path(PLOT_DIR, paste0("summary_4panel_488B_", MODE, ".pdf")), width=18, height=12)
  print((p_scrna | p_atac_proj) / (p_overlay | p_score))
  dev.off()

  cat("Plots saved to:", PLOT_DIR, "\n")
}

cat("\nDone.\n")
