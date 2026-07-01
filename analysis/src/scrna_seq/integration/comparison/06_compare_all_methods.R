#!/usr/bin/env Rscript
# 06_compare_all_methods.R
# Cross-tissue, cross-method summary of all three integration options.
# Requires Options 1-3 to have been run for both 488B and 489.
#
# Produces:
#   - Cell type purity per ATAC cluster (do clusters map to single RNA cell types?)
#   - Spatial coherence: are same predicted cell types spatially clustered?
#   - Prediction score comparison: ArchR label transfer vs Seurat bridge
#   - Cross-tissue cell type composition comparison

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
})

BASE_OUT <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration"
COMP_DIR <- file.path(BASE_OUT, "comparison")
dir.create(file.path(COMP_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)

TISSUES <- c("488B", "489")

# ── Load label transfer results ───────────────────────────────────────────────
load_tissue <- function(tissue) {
  lt_path  <- file.path(BASE_OUT, tissue, "tables", paste0("label_transfer_results_", tissue, ".csv"))
  coe_path <- file.path(BASE_OUT, tissue, "tables", paste0("coembed_method_comparison_", tissue, ".csv"))
  cor_path <- file.path(BASE_OUT, tissue, "tables", paste0("peak_gene_corr_resolution_sweep_", tissue, ".csv"))

  if (!file.exists(lt_path)) {
    message("SKIP ", tissue, ": label transfer CSV not found (", lt_path, "). Run scripts 01-03 first.")
    return(NULL)
  }

  lt  <- read.csv(lt_path)
  coe <- if (file.exists(coe_path)) read.csv(coe_path) else NULL
  cor <- if (file.exists(cor_path)) read.csv(cor_path) else NULL
  list(tissue=tissue, label_transfer=lt, coembed=coe, correlation=cor)
}

all_data <- lapply(TISSUES, load_tissue)
names(all_data) <- TISSUES
all_data <- Filter(Negate(is.null), all_data)

if (length(all_data) == 0) {
  stop("No results found. Run integration scripts 01-03 first.")
}

# ── 1. Cell type distribution per tissue ─────────────────────────────────────
ct_dist <- do.call(rbind, lapply(names(all_data), function(tissue) {
  lt <- all_data[[tissue]]$label_transfer
  if (!"predictedGroup_co" %in% colnames(lt)) return(NULL)
  n_total <- nrow(lt)
  as.data.frame(table(lt$predictedGroup_co)) %>%
    rename(cell_type=Var1, n=Freq) %>%
    mutate(tissue=tissue, pct=100*n/n_total)
}))

if (!is.null(ct_dist) && nrow(ct_dist) > 0) {
  p_dist <- ggplot(ct_dist, aes(x=reorder(cell_type, -pct), y=pct, fill=tissue)) +
    geom_col(position="dodge", alpha=0.85) +
    scale_fill_manual(values=c("488B"="steelblue","489"="tomato")) +
    labs(x="Predicted cell type", y="Fraction of ATAC spots (%)",
         title="Predicted cell type composition by tissue",
         subtitle="ArchR CCA label transfer (Option 1)") +
    theme_bw(base_size=12) +
    theme(axis.text.x=element_text(angle=35, hjust=1))
  ggsave(file.path(COMP_DIR, "plots", "celltype_composition_comparison.pdf"),
         p_dist, width=10, height=5)
  cat("Saved: celltype_composition_comparison.pdf\n")
}

# ── 2. Cell type purity per ATAC cluster ─────────────────────────────────────
purity_rows <- do.call(rbind, lapply(names(all_data), function(tissue) {
  lt <- all_data[[tissue]]$label_transfer
  if (!all(c("Clusters","predictedGroup_co") %in% colnames(lt))) return(NULL)
  lt %>%
    group_by(Clusters) %>%
    summarise(
      n = n(),
      top_ct = names(sort(table(predictedGroup_co), decreasing=TRUE))[1],
      purity = max(table(predictedGroup_co)) / n(),
      .groups = "drop"
    ) %>%
    mutate(tissue=tissue)
}))

if (!is.null(purity_rows) && nrow(purity_rows) > 0) {
  cat("\nCluster purity per tissue:\n")
  print(purity_rows %>% select(tissue, Clusters, n, top_ct, purity) %>%
          arrange(tissue, desc(purity)))

  p_purity <- ggplot(purity_rows, aes(x=purity, fill=tissue)) +
    geom_histogram(bins=20, alpha=0.7, position="identity") +
    geom_vline(data=purity_rows %>% group_by(tissue) %>% summarise(med=median(purity)),
               aes(xintercept=med, color=tissue), linetype="dashed", linewidth=1) +
    scale_fill_manual(values=c("488B"="steelblue","489"="tomato")) +
    scale_color_manual(values=c("488B"="steelblue","489"="tomato")) +
    labs(x="Cluster purity (max cell type fraction)", y="N clusters",
         title="Cell type purity per ATAC cluster") +
    theme_bw()
  ggsave(file.path(COMP_DIR, "plots", "cluster_purity.pdf"), p_purity, width=7, height=4)
  write.csv(purity_rows, file.path(COMP_DIR, "cluster_purity.csv"), row.names=FALSE)
  cat("Saved: cluster_purity.pdf + cluster_purity.csv\n")
}

# ── 3. Prediction score comparison across methods ────────────────────────────
score_rows <- do.call(rbind, lapply(names(all_data), function(tissue) {
  lt  <- all_data[[tissue]]$label_transfer
  coe <- all_data[[tissue]]$coembed

  rows <- NULL
  if ("predictedScore_co" %in% colnames(lt)) {
    rows <- rbind(rows, data.frame(
      tissue=tissue, method="ArchR CCA (Option 1)",
      score=lt$predictedScore_co
    ))
  }
  if (!is.null(coe) && "seurat_score" %in% colnames(coe)) {
    rows <- rbind(rows, data.frame(
      tissue=tissue, method="Seurat bridge (Option 3b)",
      score=coe$seurat_score
    ))
  }
  rows
}))

if (!is.null(score_rows) && nrow(score_rows) > 0) {
  p_score <- ggplot(score_rows, aes(x=score, fill=method)) +
    geom_density(alpha=0.5) +
    geom_vline(xintercept=0.5, linetype="dashed", color="red") +
    facet_wrap(~tissue) +
    scale_fill_brewer(palette="Set1") +
    labs(x="Prediction score", y="Density",
         title="Label transfer confidence: ArchR vs Seurat bridge") +
    theme_bw()
  ggsave(file.path(COMP_DIR, "plots", "prediction_score_methods.pdf"),
         p_score, width=10, height=5)
  cat("Saved: prediction_score_methods.pdf\n")

  cat("\nScore summary:\n")
  print(score_rows %>%
          group_by(tissue, method) %>%
          summarise(median=round(median(score,na.rm=TRUE),3),
                    pct_above_0.5=round(100*mean(score>0.5,na.rm=TRUE),1),
                    .groups="drop"))
}

# ── 4. Correlation sweep across tissues ──────────────────────────────────────
cor_rows <- do.call(rbind, lapply(names(all_data), function(tissue) {
  cor <- all_data[[tissue]]$correlation
  if (is.null(cor)) return(NULL)
  cor$tissue <- tissue
  cor
}))

xenium_baseline <- data.frame(
  bin_size_um=c(0,25,50,100,200,400),
  median_pearson=c(0.017,0.016,0.022,0.058,0.082,0.095),
  tissue="Xenium baseline"
)

if (!is.null(cor_rows) && nrow(cor_rows) > 0) {
  plot_cor <- bind_rows(cor_rows[,c("bin_size_um","median_pearson","tissue")], xenium_baseline)
  p_cor <- ggplot(plot_cor, aes(x=bin_size_um, y=median_pearson, color=tissue, linetype=tissue=="Xenium baseline")) +
    geom_line(linewidth=1.2) + geom_point(size=2.5) +
    scale_color_manual(values=c("488B"="steelblue","489"="tomato","Xenium baseline"="grey40")) +
    scale_linetype_manual(values=c("FALSE"="solid","TRUE"="dashed"), guide="none") +
    labs(x="Spatial bin size (µm)", y="Median per-bin Pearson",
         title="GeneScore vs imputed scRNA correlation: resolution sweep",
         subtitle="Compared to Xenium pseudobulk baseline (Option 2)",
         color="Tissue / Baseline") +
    theme_bw(base_size=12)
  ggsave(file.path(COMP_DIR, "plots", "correlation_sweep_both_tissues.pdf"),
         p_cor, width=8, height=5)
  write.csv(cor_rows, file.path(COMP_DIR, "correlation_sweep_all.csv"), row.names=FALSE)
  cat("Saved: correlation_sweep_both_tissues.pdf\n")
}

# ── 5. Overall summary table ──────────────────────────────────────────────────
summ <- do.call(rbind, lapply(names(all_data), function(tissue) {
  lt  <- all_data[[tissue]]$label_transfer
  coe <- all_data[[tissue]]$coembed
  cor <- all_data[[tissue]]$correlation

  data.frame(
    tissue=tissue,
    n_atac_cells=nrow(lt),
    n_predicted_types=if("predictedGroup_co" %in% colnames(lt)) length(unique(lt$predictedGroup_co)) else NA,
    median_archr_score=if("predictedScore_co" %in% colnames(lt)) round(median(lt$predictedScore_co,na.rm=TRUE),3) else NA,
    pct_score_above_0.5=if("predictedScore_co" %in% colnames(lt)) round(100*mean(lt$predictedScore_co>0.5,na.rm=TRUE),1) else NA,
    archr_seurat_agree=if(!is.null(coe)&&"methods_agree" %in% colnames(coe)) round(100*mean(coe$methods_agree,na.rm=TRUE),1) else NA,
    median_pearson_native=if(!is.null(cor)) round(cor$median_pearson[cor$bin_size_um==0],4) else NA,
    median_pearson_400um=if(!is.null(cor)) round(cor$median_pearson[cor$bin_size_um==400],4) else NA
  )
}))

cat("\n=== Integration Summary ===\n")
print(summ, row.names=FALSE)
write.csv(summ, file.path(COMP_DIR, "integration_summary.csv"), row.names=FALSE)
cat("\nSaved: integration_summary.csv\n")
cat("All comparison outputs in:", file.path(COMP_DIR, "plots"), "\n")
