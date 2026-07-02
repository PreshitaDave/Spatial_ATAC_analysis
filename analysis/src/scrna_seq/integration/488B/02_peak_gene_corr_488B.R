#!/usr/bin/env Rscript
# 02_peak_gene_corr_488B.R
# Option 2: Direct peak-to-gene correlation using scRNA as ground truth.
# Replaces Xenium pseudobulk (used in gene_loss_evaluation.ipynb, baseline Pearson=0.017)
# with scRNA expression imputed via GeneIntegrationMatrix from Option 1.
#
# Per ATAC spot: average GeneIntegrationMatrix across constituent ATAC cells,
# then correlate against GeneScoreMatrix. Resolution sweep: native → 400 µm.
#
# Prerequisites: 01_label_transfer_488B.R must have been run.

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(Matrix)
})

set.seed(42)
addArchRThreads(threads = 8)
addArchRGenome("hg38")

# ── Paths ────────────────────────────────────────────────────────────────────
OBJ_DIR  <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration/488B/objects"
ARCHR_INT <- file.path(OBJ_DIR, "archr_488B_with_integration")
OUT_DIR   <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration/488B"
PLOT_DIR  <- file.path(OUT_DIR, "plots", "peak_gene_corr")
TAB_DIR   <- file.path(OUT_DIR, "tables")

if (!dir.exists(ARCHR_INT)) {
  stop("Run 01_label_transfer_488B.R first: ", ARCHR_INT, " not found.")
}

cat("Loading integrated ArchR project:", ARCHR_INT, "\n")
proj <- loadArchRProject(ARCHR_INT, showLogo = FALSE)

# ── Extract matrices ──────────────────────────────────────────────────────────
cat("Extracting GeneScoreMatrix...\n")
gs  <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
gs_mat <- assay(gs)   # genes × cells
# ArchR stores gene names in rowData, not rownames
if (!is.null(rowData(gs)$name))
  rownames(gs_mat) <- make.unique(as.character(rowData(gs)$name))
gene_names_gs <- rownames(gs_mat)
cat("GeneScoreMatrix:", nrow(gs_mat), "genes ×", ncol(gs_mat), "cells\n")

cat("Loading GeneIntegrationMatrix from exported MTX (not in Arrow files)...\n")
gi_mat   <- readMM(file.path(OBJ_DIR, "archr_488B_GeneIntegrationMatrix.mtx"))
gi_genes <- read.csv(file.path(OBJ_DIR, "archr_488B_gene_names.csv"), header=TRUE)[,1]
gi_cells <- read.csv(file.path(OBJ_DIR, "archr_488B_cell_names.csv"), header=TRUE)[,1]
rownames(gi_mat) <- gi_genes
colnames(gi_mat) <- gi_cells
gi_mat <- as(gi_mat, "CsparseMatrix")
cat("GeneIntegrationMatrix:", nrow(gi_mat), "genes ×", ncol(gi_mat), "cells\n")

# Shared genes
shared_genes <- intersect(gene_names_gs, rownames(gi_mat))
cat("Shared genes between GeneScore and GeneIntegration:", length(shared_genes), "\n")
gs_mat <- gs_mat[shared_genes, ]
gi_mat <- gi_mat[shared_genes, ]

# Align cells: keep only cells present in both matrices
common_cells <- intersect(colnames(gs_mat), colnames(gi_mat))
cat("Common cells:", length(common_cells), "\n")
gs_mat <- gs_mat[, common_cells]
gi_mat <- gi_mat[, common_cells]

# ── Spatial coordinates per ATAC cell ────────────────────────────────────────
# ArchR stores spatial coords in cellColData if added; otherwise from nonlinear_aligned h5ad
meta <- as.data.frame(getCellColData(proj))
spatial_cols <- grep("spatial|x_um|y_um|X|Y", colnames(meta), value = TRUE, ignore.case = TRUE)
cat("Potential spatial columns:", paste(spatial_cols, collapse=", "), "\n")

# If spatial coords are not in cellColData, use the nonlinear-aligned h5ad
if (length(spatial_cols) < 2) {
  cat("Loading spatial coordinates from nonlinear aligned h5ad...\n")
  spatial_load_ok <- tryCatch({
    library(anndata)
    atac_h5ad <- read_h5ad("/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/alignment/mosaicfield_outputs/atac_nonlinear_aligned.h5ad")
    sp_coords <- atac_h5ad$obsm[["spatial"]]
    rownames(sp_coords) <- atac_h5ad$obs_names
    colnames(sp_coords) <- c("x_um", "y_um")
    # Match to ArchR cell names (strip sample prefix if needed)
    cell_bare <- sub(".*#", "", proj$cellNames)
    sp_coords_matched <- sp_coords[cell_bare, , drop = FALSE]
    meta$x_um <<- sp_coords_matched[, 1]
    meta$y_um <<- sp_coords_matched[, 2]
    TRUE
  }, error = function(e) {
    cat("  WARNING: failed to load spatial coordinates from h5ad (", conditionMessage(e), ")\n")
    cat("  Native-resolution correlation is unaffected; resolution sweep (bin_size > 0) will be skipped.\n")
    meta$x_um <<- NA_real_
    meta$y_um <<- NA_real_
    FALSE
  })
} else {
  meta$x_um <- meta[[spatial_cols[1]]]
  meta$y_um <- meta[[spatial_cols[2]]]
}

# ── Compute per-gene Pearson: GeneScore vs GeneIntegration, across units ───
# Correlate across CELLS (or spatial BINS), one r per gene — matches the
# gene_loss_evaluation.ipynb Xenium-baseline methodology (per-gene, across
# spots/cells). Correlating across GENES within one cell/bin (the previous
# approach) is dominated by absolute expression-magnitude differences between
# GeneScore and imputed-expression units and is not comparable to the
# Xenium baseline.
pergene_pearson <- function(gs_units, gi_units) {
  # gs_units, gi_units: genes x units matrices, same gene/unit order
  n_units <- ncol(gs_units)
  if (n_units < 3) return(rep(NA_real_, nrow(gs_units)))
  gs_m <- as.matrix(gs_units)
  gi_m <- as.matrix(gi_units)
  vapply(seq_len(nrow(gs_m)), function(i) {
    x <- gs_m[i, ]; y <- gi_m[i, ]
    if (var(x) > 0 && var(y) > 0) cor(x, y, method = "pearson") else NA_real_
  }, numeric(1))
}

cat("Computing per-gene Pearson correlation across cells (native resolution)...\n")
n_cells <- ncol(gs_mat)
r_gene_native <- pergene_pearson(gs_mat, gi_mat)
cat("Native resolution (per gene, across ATAC cells):\n")
cat("  Median Pearson:", round(median(r_gene_native, na.rm=TRUE), 4), "\n")
cat("  Mean Pearson:  ", round(mean(r_gene_native, na.rm=TRUE), 4), "\n")
cat("  Xenium baseline: 0.017\n")

# ── Resolution sweep ──────────────────────────────────────────────────────────
# Aggregate cells into spatial bins, then compute per-gene Pearson across bins.
bin_sizes <- c(0, 25, 50, 100, 200, 400)   # µm
sweep_results <- data.frame(bin_size_um = integer(), median_pearson = numeric(),
                             mean_pearson = numeric(), n_bins = integer())

for (bin_size in bin_sizes) {
  if (bin_size == 0) {
    # Native: per-gene across cells
    r_vals <- r_gene_native
    n_b <- n_cells
  } else {
    # Aggregate into spatial grid bins
    meta$bin_x <- floor(meta$x_um / bin_size)
    meta$bin_y <- floor(meta$y_um / bin_size)
    meta$bin_id <- paste(meta$bin_x, meta$bin_y, sep="_")
    bins <- unique(meta$bin_id[!is.na(meta$x_um)])
    gs_binned <- matrix(NA_real_, nrow = nrow(gs_mat), ncol = length(bins))
    gi_binned <- matrix(NA_real_, nrow = nrow(gi_mat), ncol = length(bins))
    for (b_idx in seq_along(bins)) {
      b <- bins[b_idx]
      cells_in_bin <- which(meta$bin_id == b & !is.na(meta$x_um))
      if (length(cells_in_bin) < 1) next
      gs_binned[, b_idx] <- rowMeans(as.matrix(gs_mat[, cells_in_bin, drop=FALSE]))
      gi_binned[, b_idx] <- rowMeans(as.matrix(gi_mat[, cells_in_bin, drop=FALSE]))
    }
    r_vals <- pergene_pearson(gs_binned, gi_binned)
    n_b <- length(bins)
  }
  sweep_results <- rbind(sweep_results, data.frame(
    bin_size_um    = bin_size,
    median_pearson = median(r_vals, na.rm=TRUE),
    mean_pearson   = mean(r_vals, na.rm=TRUE),
    n_bins         = n_b
  ))
  cat(sprintf("  %3d µm bins: median Pearson = %.4f  (n=%d bins)\n",
              bin_size, median(r_vals, na.rm=TRUE), n_b))
}

write.csv(sweep_results, file.path(TAB_DIR, "peak_gene_corr_resolution_sweep_488B.csv"), row.names=FALSE)

# Xenium baseline from gene_loss_evaluation
xenium_baseline <- data.frame(
  bin_size_um = c(0, 25, 50, 100, 200, 400),
  median_pearson = c(0.017, 0.016, 0.022, 0.058, 0.082, 0.095)
)

p_sweep <- ggplot() +
  geom_line(data = sweep_results,  aes(x=bin_size_um, y=median_pearson, color="scRNA (this analysis)"), linewidth=1.2) +
  geom_point(data = sweep_results, aes(x=bin_size_um, y=median_pearson, color="scRNA (this analysis)"), size=3) +
  geom_line(data = xenium_baseline, aes(x=bin_size_um, y=median_pearson, color="Xenium baseline"), linewidth=1.2, linetype="dashed") +
  geom_point(data = xenium_baseline, aes(x=bin_size_um, y=median_pearson, color="Xenium baseline"), size=3) +
  scale_color_manual(values = c("scRNA (this analysis)" = "steelblue", "Xenium baseline" = "tomato")) +
  labs(x = "Spatial bin size (µm)", y = "Median per-bin Pearson",
       title = "488B: GeneScore vs imputed scRNA — resolution sweep",
       subtitle = "Compared to Xenium-based baseline (nonlinear+Voronoi)",
       color = NULL) +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")
ggsave(file.path(PLOT_DIR, "correlation_resolution_sweep.pdf"), p_sweep, width=8, height=5)

# ── Top/bottom correlated genes ───────────────────────────────────────────────
gs_mean <- rowMeans(as.matrix(gs_mat))
gi_mean <- rowMeans(as.matrix(gi_mat))
gene_cors <- sapply(shared_genes, function(g) {
  cor(as.numeric(gs_mat[g, ]), as.numeric(gi_mat[g, ]), method = "pearson")
})
gene_df <- data.frame(gene=shared_genes, pearson=gene_cors,
                      mean_genescore=gs_mean, mean_imputed=gi_mean) %>%
  arrange(desc(pearson))
write.csv(gene_df, file.path(TAB_DIR, "per_gene_pearson_488B.csv"), row.names=FALSE)

top20 <- head(gene_df, 20)
bot20 <- tail(gene_df, 20)
plot_genes <- rbind(
  data.frame(top20, rank = "Top 20 (most correlated)"),
  data.frame(bot20, rank = "Bottom 20 (least correlated)")
)
p_genes <- ggplot(plot_genes, aes(x = reorder(gene, pearson), y = pearson, fill = rank)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~rank, scales="free_y") +
  scale_fill_manual(values = c("Top 20 (most correlated)"="steelblue",
                                "Bottom 20 (least correlated)"="tomato")) +
  labs(x = NULL, y = "Pearson r (GeneScore vs imputed scRNA)",
       title = "488B: top/bottom correlated genes") +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(PLOT_DIR, "top_bottom_genes.pdf"), p_genes, width=10, height=8)

cat("\nOption 2 (peak-gene correlation) complete for 488B.\n")
cat("Results in:", OUT_DIR, "\n")
