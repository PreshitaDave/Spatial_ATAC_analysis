#!/usr/bin/env Rscript
# 05_compare_methods_489.R
# Compare all cell type assignment methods for lowseq_489:
#   - ArchR label transfer: full vs balanced reference
#   - Seurat bridge (03_coembed 3b): full vs balanced reference
# Outputs: confusion matrices, score distributions, spatial maps, summary table.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
})

set.seed(42)

TISSUE   <- "lowseq_489"
TAB_DIR  <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration/489/tables"
PLOT_DIR <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration/489/plots/method_comparison"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Consistent colour palette ──────────────────────────────────────────────────
ct_colors <- c(
  Tumor       = "#E41A1C",
  T_cell      = "#377EB8",
  NK_cell     = "#FF7F00",
  B_cell      = "#4DAF4A",
  Myeloid     = "#984EA3",
  Fibroblast  = "#A65628",
  Endothelial = "#F781BF",
  Unknown     = "#AAAAAA"
)

# ── Load data ──────────────────────────────────────────────────────────────────
cat("Loading tables...\n")

OBJ_DIR <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration/489/objects"

lt_full <- read.csv(file.path(TAB_DIR, "label_transfer_results_489.csv"),
                    row.names = 1, check.names = FALSE)
lt_bal  <- read.csv(file.path(TAB_DIR, "label_transfer_results_489_balanced.csv"),
                    row.names = 1, check.names = FALSE)

# Load Seurat bridge predictions directly from RDS objects (the CSVs had degenerate
# scores from an earlier run with k.weight=1; the RDS objects have corrected predictions)
cat("Loading Seurat bridge objects for correct predictions...\n")
bridge_full <- readRDS(file.path(OBJ_DIR, "seurat_atac_bridge_489.rds"))
bridge_bal  <- readRDS(file.path(OBJ_DIR, "seurat_atac_bridge_489_balanced.rds"))

cat("Seurat bridge (full) cell type table:\n")
print(table(bridge_full$predicted.id))
cat("Seurat bridge (balanced) cell type table:\n")
print(table(bridge_bal$predicted.id))

# Merge all 4 methods per cell using ArchR label transfer as the base (has spatial coords)
base <- lt_bal[, c("x_spatial","y_spatial","Clusters","predictedGroup_co","predictedScore_co")]
colnames(base) <- c("x","y","atac_cluster","archr_balanced","score_archr_balanced")

base$archr_full       <- lt_full[rownames(base), "predictedGroup_co"]
base$score_archr_full <- lt_full[rownames(base), "predictedScore_co"]

# Strip tissue prefix from Seurat bridge cell names to match ArchR rownames
# ArchR: "Lowseq_489#AAACGAAAGAACCGAA-1"  Seurat: "AAACGAAAGAACCGAA-1"
bridge_full_pred  <- setNames(bridge_full$predicted.id,          colnames(bridge_full))
bridge_full_score <- setNames(bridge_full$prediction.score.max,  colnames(bridge_full))
bridge_bal_pred   <- setNames(bridge_bal$predicted.id,           colnames(bridge_bal))
bridge_bal_score  <- setNames(bridge_bal$prediction.score.max,   colnames(bridge_bal))

# Try matching with and without tissue prefix
archr_cells <- rownames(base)
# strip prefix "Lowseq_489#" if present
seurat_cells_from_archr <- sub("^[^#]+#", "", archr_cells)

base$seurat_full          <- bridge_full_pred[seurat_cells_from_archr]
base$score_seurat_full    <- as.numeric(bridge_full_score[seurat_cells_from_archr])
base$seurat_balanced      <- bridge_bal_pred[seurat_cells_from_archr]
base$score_seurat_balanced <- as.numeric(bridge_bal_score[seurat_cells_from_archr])

# Fallback: try direct match if prefix stripping didn't work
if (all(is.na(base$seurat_full))) {
  base$seurat_full          <- bridge_full_pred[archr_cells]
  base$score_seurat_full    <- as.numeric(bridge_full_score[archr_cells])
  base$seurat_balanced      <- bridge_bal_pred[archr_cells]
  base$score_seurat_balanced <- as.numeric(bridge_bal_score[archr_cells])
}
cat("Seurat full matched:", sum(!is.na(base$seurat_full)), "of", nrow(base), "\n")
cat("Seurat balanced matched:", sum(!is.na(base$seurat_balanced)), "of", nrow(base), "\n")

cat("Cells with all 4 assignments:", sum(complete.cases(base[,c("archr_full","archr_balanced","seurat_full","seurat_balanced")])), "\n")

methods <- c("archr_full","archr_balanced","seurat_full","seurat_balanced")
method_labels <- c("ArchR (full)","ArchR (balanced)","Seurat bridge (full)","Seurat bridge (balanced)")
score_cols    <- c("score_archr_full","score_archr_balanced","score_seurat_full","score_seurat_balanced")

# ── 1. Cell type proportion bar chart — all 4 methods ─────────────────────────
cat("Plotting cell type proportions...\n")

prop_df <- lapply(seq_along(methods), function(i) {
  x <- base[[methods[i]]]
  x[is.na(x)] <- "Unknown"
  ct <- as.data.frame(table(x), stringsAsFactors = FALSE)
  colnames(ct) <- c("cell_type","n")
  ct$pct    <- 100 * ct$n / sum(ct$n)
  ct$method <- method_labels[i]
  ct
}) %>% bind_rows()

p_prop <- ggplot(prop_df, aes(x=method, y=pct, fill=cell_type)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=ct_colors, name="Cell type") +
  labs(x=NULL, y="% of ATAC spots",
       title=paste0(TISSUE,": cell type composition by method")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1, size=10))
ggsave(file.path(PLOT_DIR,"01_celltype_proportions.pdf"), p_prop, width=10, height=6)

# ── 2. Prediction score distributions — violin by method ──────────────────────
cat("Plotting score distributions...\n")

score_df <- lapply(seq_along(methods), function(i) {
  data.frame(
    method    = method_labels[i],
    cell_type = base[[methods[i]]],
    score     = base[[score_cols[i]]],
    stringsAsFactors = FALSE
  )
}) %>% bind_rows() %>% filter(!is.na(score), !is.na(cell_type))

# Overall score per method
p_score_overall <- ggplot(score_df, aes(x=method, y=score, fill=method)) +
  geom_violin(alpha=0.7, draw_quantiles=c(0.25,0.5,0.75)) +
  geom_hline(yintercept=0.5, linetype="dashed", color="red") +
  geom_hline(yintercept=0.7, linetype="dashed", color="blue") +
  scale_fill_brewer(palette="Set2") +
  labs(x=NULL, y="Prediction score",
       title=paste0(TISSUE,": prediction score by method"),
       caption="Red=0.5, Blue=0.7 thresholds") +
  theme_bw() + theme(axis.text.x=element_text(angle=35,hjust=1), legend.position="none")
ggsave(file.path(PLOT_DIR,"02_score_violin_by_method.pdf"), p_score_overall, width=9, height=5)

# Score by cell type for each method
p_score_ct <- ggplot(score_df %>% filter(!is.na(cell_type)),
                     aes(x=cell_type, y=score, fill=cell_type)) +
  geom_boxplot(outlier.size=0.3, alpha=0.8) +
  geom_hline(yintercept=0.5, linetype="dashed", color="red") +
  scale_fill_manual(values=ct_colors) +
  facet_wrap(~method, ncol=2) +
  labs(x=NULL, y="Prediction score",
       title=paste0(TISSUE,": score by cell type and method")) +
  theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none")
ggsave(file.path(PLOT_DIR,"03_score_by_celltype_method.pdf"), p_score_ct, width=12, height=8)

# ── 3. Summary stats table ─────────────────────────────────────────────────────
cat("Building summary table...\n")

summary_tbl <- lapply(seq_along(methods), function(i) {
  sc   <- base[[score_cols[i]]]
  pred <- base[[methods[i]]]
  top  <- names(sort(table(pred[!is.na(pred)]), decreasing=TRUE))[1]
  data.frame(
    method        = method_labels[i],
    n_assigned    = sum(!is.na(pred)),
    top_celltype  = top,
    pct_top       = round(100*sum(pred==top, na.rm=TRUE)/sum(!is.na(pred)),1),
    mean_score    = round(mean(sc, na.rm=TRUE),3),
    pct_gt0.5     = round(100*mean(sc>0.5, na.rm=TRUE),1),
    pct_gt0.7     = round(100*mean(sc>0.7, na.rm=TRUE),1),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

cat("\nMethod comparison summary:\n"); print(summary_tbl)
write.csv(summary_tbl, file.path(TAB_DIR,"method_comparison_summary_489.csv"), row.names=FALSE)

# ── 4. Confusion matrices between methods ─────────────────────────────────────
cat("Building confusion matrices...\n")

plot_confusion <- function(pred_a, pred_b, label_a, label_b, title, filename) {
  both <- complete.cases(data.frame(pred_a, pred_b))
  ct <- table(A=pred_a[both], B=pred_b[both])
  pct <- sweep(ct, 1, rowSums(ct), "/")
  all_types <- union(rownames(pct), colnames(pct))
  # Add % labels
  df <- as.data.frame(as.table(pct))
  colnames(df) <- c("MethodA","MethodB","Fraction")
  df$Count <- as.vector(ct)
  df$Label <- ifelse(df$Fraction > 0.05, paste0(round(df$Fraction*100), "%"), "")

  p <- ggplot(df, aes(x=MethodB, y=MethodA, fill=Fraction)) +
    geom_tile(color="white") +
    geom_text(aes(label=Label), size=3) +
    scale_fill_gradient2(low="white", mid="#FED976", high="#E31A1C",
                         midpoint=0.4, limits=c(0,1), name="Row %") +
    labs(x=label_b, y=label_a, title=title) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45,hjust=1),
          panel.grid=element_blank())
  ggsave(file.path(PLOT_DIR, filename), p, width=8, height=7)
  invisible(p)
}

# a) ArchR full vs ArchR balanced — shows impact of balancing
plot_confusion(base$archr_full, base$archr_balanced,
               "ArchR (full ref)", "ArchR (balanced ref)",
               paste0(TISSUE,": ArchR full vs balanced (row = full)"),
               "04_confusion_archr_full_vs_balanced.pdf")

# b) ArchR balanced vs Seurat bridge balanced — inter-method agreement
plot_confusion(base$archr_balanced, base$seurat_balanced,
               "ArchR (balanced)", "Seurat bridge (balanced)",
               paste0(TISSUE,": ArchR balanced vs Seurat bridge balanced"),
               "05_confusion_archr_vs_seurat_balanced.pdf")

# c) ArchR full vs Seurat bridge full
plot_confusion(base$archr_full, base$seurat_full,
               "ArchR (full)", "Seurat bridge (full)",
               paste0(TISSUE,": ArchR full vs Seurat bridge full"),
               "06_confusion_archr_vs_seurat_full.pdf")

# d) Seurat full vs Seurat balanced — impact of balancing on Seurat
plot_confusion(base$seurat_full, base$seurat_balanced,
               "Seurat bridge (full)", "Seurat bridge (balanced)",
               paste0(TISSUE,": Seurat bridge full vs balanced"),
               "07_confusion_seurat_full_vs_balanced.pdf")

# ── 5. Pairwise agreement rates table ─────────────────────────────────────────
pairs <- list(
  c("archr_full","archr_balanced"),
  c("archr_balanced","seurat_balanced"),
  c("archr_full","seurat_full"),
  c("seurat_full","seurat_balanced"),
  c("archr_full","seurat_balanced")
)
agree_tbl <- lapply(pairs, function(p) {
  both <- complete.cases(base[,p])
  agr  <- mean(base[[p[1]]][both] == base[[p[2]]][both], na.rm=TRUE)
  data.frame(method_A=p[1], method_B=p[2],
             n_cells=sum(both),
             pct_agree=round(100*agr,1))
}) %>% bind_rows()
cat("\nPairwise agreement rates:\n"); print(agree_tbl)
write.csv(agree_tbl, file.path(TAB_DIR,"method_pairwise_agreement_489.csv"), row.names=FALSE)

# ── 6. Spatial maps — predicted cell type per method ──────────────────────────
cat("Plotting spatial maps...\n")

spatial_plots <- lapply(seq_along(methods), function(i) {
  df <- base %>% filter(!is.na(.data[[methods[i]]]))
  ggplot(df, aes(x=x, y=y, color=.data[[methods[i]]])) +
    geom_point(size=0.6, alpha=0.8) +
    scale_color_manual(values=ct_colors, name="Cell type") +
    ggtitle(method_labels[i]) +
    coord_fixed() + theme_bw() +
    theme(axis.title=element_blank(), axis.text=element_blank(),
          axis.ticks=element_blank(), panel.grid=element_blank(),
          legend.position="right",
          plot.title=element_text(size=10, face="bold")) +
    guides(color=guide_legend(override.aes=list(size=2.5)))
})

pdf(file.path(PLOT_DIR,"08_spatial_maps_all_methods.pdf"), width=18, height=10)
print((spatial_plots[[1]] | spatial_plots[[2]]) / (spatial_plots[[3]] | spatial_plots[[4]]))
dev.off()

# ── 7. Score threshold analysis — what % high-confidence per method ────────────
cat("Plotting score threshold analysis...\n")

thresh_df <- lapply(seq_along(methods), function(i) {
  sc <- base[[score_cols[i]]]
  thresholds <- seq(0.1, 0.9, by=0.05)
  data.frame(
    method    = method_labels[i],
    threshold = thresholds,
    pct_above = sapply(thresholds, function(t) 100*mean(sc >= t, na.rm=TRUE))
  )
}) %>% bind_rows()

p_thresh <- ggplot(thresh_df, aes(x=threshold, y=pct_above, color=method)) +
  geom_line(size=1.2) + geom_point(size=2) +
  geom_vline(xintercept=c(0.5,0.7), linetype="dashed", color="grey50") +
  scale_color_brewer(palette="Set1", name="Method") +
  labs(x="Score threshold", y="% spots above threshold",
       title=paste0(TISSUE,": fraction of spots above score threshold by method")) +
  theme_bw()
ggsave(file.path(PLOT_DIR,"09_score_threshold_curve.pdf"), p_thresh, width=9, height=5)

# ── 8. Per-cluster cell type composition (ArchR balanced — best method) ────────
cat("Plotting cluster composition...\n")

clust_df <- base %>%
  filter(!is.na(archr_balanced), !is.na(atac_cluster)) %>%
  group_by(atac_cluster, archr_balanced) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(atac_cluster) %>%
  mutate(pct=100*n/sum(n)) %>%
  ungroup()

p_clust <- ggplot(clust_df, aes(x=atac_cluster, y=pct, fill=archr_balanced)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=ct_colors, name="Cell type") +
  labs(x="ATAC cluster", y="% of spots",
       title=paste0(TISSUE,": cell type composition per ATAC cluster (ArchR balanced)")) +
  theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(file.path(PLOT_DIR,"10_cluster_celltype_composition.pdf"), p_clust, width=10, height=6)

cat("\nAll plots saved to:", PLOT_DIR, "\n")
cat("Done.\n")
