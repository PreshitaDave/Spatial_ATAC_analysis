#!/usr/bin/env Rscript
# ============================================================================
# annotate_clusters_by_markers.R
#
# Assign putative cell-type identities to clusters from all 3 spatial-
# smoothing methods (Baseline / Alpha-Blend / Iterative), using an expanded
# breast-cancer reference marker panel and a module-score approach (mean
# z-scored GeneScore across each identity's marker set, not a single noisy
# top gene) to be robust to GeneScoreMatrix sparsity.
#
# Usage:
#   Rscript annotate_clusters_by_markers.R <tissue>
#   Example: Rscript annotate_clusters_by_markers.R lowseq_489
# ============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

log_msg <- function(tag, msg) cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), tag, msg))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript annotate_clusters_by_markers.R <tissue>")
tissue <- args[1]
log_msg("start", sprintf("===== Annotating clusters by markers: %s =====", tissue))

set.seed(42)
addArchRGenome("hg38")
addArchRThreads(threads = 4)

project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
proj_dir <- file.path(project_root, "analysis/binsize_comparison/archr_projects",
                      sprintf("%s_5000bp_binarizeFALSE", tissue))
output_dir <- file.path(project_root, "analysis/binsize_comparison/spatial_smoothing")
assign_csv <- file.path(output_dir, sprintf("%s_5000bp_3method_cluster_assignments.csv", tissue))

if (!dir.exists(proj_dir)) stop(sprintf("Project directory not found: %s", proj_dir))
if (!file.exists(assign_csv)) stop(sprintf("Cluster assignment CSV not found (run compare_spatial_smoothing_methods.R first): %s", assign_csv))

# ============================================================================
# Reference marker panel — established breast-cancer single-cell/spatial
# atlas grouping (cf. Wu et al. 2021 Nat Genetics breast cancer atlas).
# Reference-anchored (not reference-free): this is what makes labels like
# "T cell" human-readable at all. The robustness fix vs the project's prior
# approach is the module score below (mean z-score across a gene set),
# not a single noisy top-DE-gene call.
# ============================================================================
marker_panel <- list(
  "Epithelial/Tumor (general)"       = c("EPCAM","KRT8","KRT18","KRT19","CDH1","ELF3"),
  "Luminal hormone-receptor+ tumor"  = c("ESR1","PGR","FOXA1","GATA3"),
  "HER2+ tumor"                      = c("ERBB2","GRB7"),
  "Basal/myoepithelial"              = c("KRT5","KRT14","TP63","ACTA2","OXTR"),
  "Fibroblast/CAF"                   = c("COL1A1","COL1A2","PDGFRB","FAP","DCN","LUM"),
  "Endothelial"                      = c("PECAM1","VWF","CDH5","CLDN5"),
  "T cell (general)"                 = c("CD3D","CD3E","CD2"),
  "Cytotoxic T cell"                 = c("CD8A"),
  "Helper/Regulatory T cell"         = c("CD4","FOXP3","CTLA4","IL2RA"),
  "B cell"                           = c("CD19","MS4A1","CD79A"),
  "Plasma cell"                      = c("MZB1","JCHAIN","IGHG1"),
  "Macrophage/Myeloid"               = c("CD68","CD14","ITGAM","LYZ","CSF1R"),
  "Dendritic cell"                   = c("CD1C","CLEC9A"),
  "Mast cell"                        = c("KIT","TPSAB1","CPA3"),
  "Proliferating"                    = c("MKI67","TOP2A","PCNA"),
  "DNA-repair-deficient tumor"       = c("BRCA1","BRCA2")
)
all_genes <- unique(unlist(marker_panel))
log_msg("step", sprintf("Reference panel: %d identities, %d unique genes", length(marker_panel), length(all_genes)))

log_msg("step", sprintf("Loading ArchR project from: %s", proj_dir))
proj <- loadArchRProject(path = proj_dir, force = TRUE)
log_msg("step", sprintf("Loaded project with %d cells", nrow(proj@cellColData)))

# ============================================================================
# Pull GeneScoreMatrix restricted to the reference panel
# ============================================================================
log_msg("step", "Pulling GeneScoreMatrix for reference panel genes...")
gsm <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
gene_names <- rowData(gsm)$name
present <- all_genes %in% gene_names
if (any(!present)) {
  log_msg("warn", sprintf("Genes not found in ArchR gene annotation, dropped: %s",
                          paste(all_genes[!present], collapse = ", ")))
}
genes_used <- all_genes[present]
gsm_mat <- assay(gsm)[match(genes_used, gene_names), , drop = FALSE]
rownames(gsm_mat) <- genes_used
colnames(gsm_mat) <- colnames(gsm)
log_msg("step", sprintf("GeneScoreMatrix subset: %d genes x %d cells", nrow(gsm_mat), ncol(gsm_mat)))

# ============================================================================
# Load per-cell cluster assignments (3 methods) and join by cellName
# ============================================================================
assign_df <- read.csv(assign_csv, stringsAsFactors = FALSE)
log_msg("step", sprintf("Loaded %d cells with 3-method cluster assignments", nrow(assign_df)))

match_idx <- match(assign_df$cellName, colnames(gsm_mat))
n_matched <- sum(!is.na(match_idx))
log_msg("step", sprintf("Matched %d/%d assignment cells to GeneScoreMatrix", n_matched, nrow(assign_df)))
assign_df <- assign_df[!is.na(match_idx), ]
match_idx <- match_idx[!is.na(match_idx)]
gsm_sub <- gsm_mat[, match_idx, drop = FALSE]  # columns now aligned to assign_df rows

order_clusters_factor <- function(clusters) {
  lvls <- paste0("C", sort(unique(as.integer(gsub("[^0-9]", "", clusters)))))
  factor(clusters, levels = lvls)
}

# ============================================================================
# Per method: mean GeneScore per (cluster x gene), z-score per gene across
# clusters, module score per identity, assign top identity + confidence gap.
# ============================================================================
annotate_method <- function(method_col, method_label) {
  clusters <- order_clusters_factor(assign_df[[method_col]])
  cluster_levels <- levels(clusters)

  mean_mat <- sapply(cluster_levels, function(cl) {
    idx <- which(clusters == cl)
    Matrix::rowMeans(gsm_sub[, idx, drop = FALSE])
  })
  rownames(mean_mat) <- genes_used
  colnames(mean_mat) <- cluster_levels

  z_mat <- t(scale(t(mean_mat)))  # z-score each gene across clusters
  z_mat[is.na(z_mat)] <- 0

  module_scores <- sapply(marker_panel, function(genes) {
    g <- intersect(genes, rownames(z_mat))
    if (length(g) == 0) return(rep(NA, ncol(z_mat)))
    colMeans(z_mat[g, , drop = FALSE])
  })
  rownames(module_scores) <- cluster_levels  # clusters x identities

  results <- lapply(cluster_levels, function(cl) {
    scores <- sort(module_scores[cl, ], decreasing = TRUE)
    n_cells <- sum(clusters == cl)
    data.frame(
      Method = method_label, Cluster = cl, N_cells = n_cells,
      Assigned_Identity = names(scores)[1], Module_Score = round(scores[1], 3),
      Second_Identity = names(scores)[2], Second_Module_Score = round(scores[2], 3),
      Confidence_Gap = round(scores[1] - scores[2], 3)
    )
  })
  identity_df <- do.call(rbind, results)

  list(identity_df = identity_df, z_mat = z_mat, cluster_levels = cluster_levels)
}

res_baseline <- annotate_method("Baseline", "Baseline")
res_blend <- annotate_method("Alpha_Blend", "Alpha-Blend")
res_iterative <- annotate_method("Iterative", "Iterative")

identity_all <- rbind(res_baseline$identity_df, res_blend$identity_df, res_iterative$identity_df)
identity_csv <- file.path(output_dir, sprintf("%s_cluster_identity.csv", tissue))
write.csv(identity_all, identity_csv, row.names = FALSE)
log_msg("step", sprintf("Saved cluster identity table: %s", identity_csv))

# ============================================================================
# Heatmap: cluster x gene z-scored GeneScore, one panel per method
# ============================================================================
plot_identity_heatmap <- function(z_mat, title) {
  df <- as.data.frame(z_mat)
  df$Gene <- rownames(df)
  df_long <- pivot_longer(df, -Gene, names_to = "Cluster", values_to = "Z")
  df_long$Cluster <- factor(df_long$Cluster, levels = colnames(z_mat))
  gene_identity <- sapply(df_long$Gene, function(g) {
    hits <- names(marker_panel)[sapply(marker_panel, function(x) g %in% x)]
    if (length(hits) == 0) return("Other")
    hits[1]
  })
  gene_order <- rownames(z_mat)[order(sapply(rownames(z_mat), function(g) {
    hits <- which(sapply(marker_panel, function(x) g %in% x))
    if (length(hits) == 0) length(marker_panel) + 1 else hits[1]
  }))]
  df_long$Gene <- factor(df_long$Gene, levels = gene_order)
  ggplot(df_long, aes(x = Cluster, y = Gene, fill = Z)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "z-score") +
    labs(title = title) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 6), plot.title = element_text(size = 10, face = "bold"))
}

p1 <- plot_identity_heatmap(res_baseline$z_mat, "Baseline")
p2 <- plot_identity_heatmap(res_blend$z_mat, "Alpha-Blend")
p3 <- plot_identity_heatmap(res_iterative$z_mat, "Iterative")

heatmap_pdf <- file.path(output_dir, sprintf("%s_cluster_identity_heatmap.pdf", tissue))
ggsave(heatmap_pdf, (p1 | p2 | p3) + patchwork::plot_annotation(title = sprintf("%s: cluster identity marker z-scores", tissue)),
      width = 20, height = 10)
log_msg("step", sprintf("Saved identity heatmap: %s", heatmap_pdf))

log_msg("done", sprintf("Completed: %s (cluster identity annotation)", tissue))
