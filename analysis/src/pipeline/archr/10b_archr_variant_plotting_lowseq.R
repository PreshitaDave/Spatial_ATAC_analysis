#!/usr/bin/env Rscript
# 10b_archr_variant_plotting_lowseq.R
# Plot somatic variants on ArchR spatial + UMAP embeddings for LOWSEQ data
# Lowseq project must be rebuilt from arrow file since it was never saved

cat("=== LowSeq ArchR Somatic Variant Plotting ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# --- Setup ---
.libPaths(c('/projectnb/paxlab/presh/env/R_4.4/ArchR_libs', .libPaths()))

library(ArchR)
library(ggplot2)
library(parallel)
set.seed(1)
addArchRGenome("hg38")
addArchRThreads(threads = as.integer(Sys.getenv("NSLOTS", "6")))

out_dir <- "/projectnb/paxlab/presh/projects/spatial_atac/analysis/comparison/somatic"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================
# SECTION 1: Build lowseq ArchR project from arrow file
# =============================================================
lowseq_save_dir <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/lowseq_saveArchR"
save_rds <- file.path(lowseq_save_dir, "Save-ArchR-Project.rds")

if (file.exists(save_rds)) {
  cat("--- Loading existing lowseq ArchR project ---\n")
  lowseq.proj <- loadArchRProject(lowseq_save_dir)
  cat("Loaded:", nrow(getCellColData(lowseq.proj)), "cells\n")
  cat("Embeddings:", paste(names(lowseq.proj@embeddings), collapse = ", "), "\n")
  cat("CellColData cols:", paste(colnames(getCellColData(lowseq.proj)), collapse = ", "), "\n")
} else {
  cat("--- Building lowseq ArchR project ---\n")

  arrow_path <- "/projectnb/paxlab/yeting/SpatialATACseq/results/preprocessed/D01942_filtered_ArchR_Output/ArrowFiles"
  lowseq_arrow <- file.path(arrow_path, "Lowseq.arrow")

  dir.create(lowseq_save_dir, recursive = TRUE, showWarnings = FALSE)
  lowseq.proj <- ArchRProject(
    ArrowFiles = lowseq_arrow,
    outputDirectory = lowseq_save_dir,
    copyArrows = FALSE
  )
  cat("Lowseq project created:", nrow(getCellColData(lowseq.proj)), "cells\n")

  # --- Add spatial coordinates ---
  cat("Adding spatial coordinates...\n")
  spatial_locs <- read.csv("/projectnb/paxlab/yeting/SpatialATACseq/data/tissue_positions_list.csv",
                           header = FALSE)
  colnames(spatial_locs) <- c("barcode", "in_tissue", "array_row", "array_col", "x_spatial", "y_spatial")
  xy2 <- spatial_locs[spatial_locs$in_tissue == 1, ]
  xy2$tissue <- ifelse(xy2$y_spatial > 4000, "Tiss_488B", "Tiss_487")

  lowseq.meta <- data.frame(lowseq.proj@cellColData)
  lowseq.meta$barcode <- gsub("-1", "", gsub("Lowseq#", "", rownames(lowseq.meta)))
  lowseq.meta2 <- merge(lowseq.meta, xy2, by = "barcode", all.x = TRUE)
  match_indices <- match(lowseq.meta$barcode, lowseq.meta2$barcode)
  lowseq.meta2_ord <- lowseq.meta2[match_indices, ]

  lowseq.proj <- addCellColData(lowseq.proj, data = lowseq.meta2_ord$tissue,
                                 name = "tissue", cells = getCellNames(lowseq.proj))
  lowseq.proj <- addCellColData(lowseq.proj, data = lowseq.meta2_ord$x_spatial,
                                 name = "x_spatial", cells = getCellNames(lowseq.proj))
  lowseq.proj <- addCellColData(lowseq.proj, data = lowseq.meta2_ord$y_spatial,
                                 name = "y_spatial", cells = getCellNames(lowseq.proj))

  cat("Tissue distribution:\n")
  print(table(lowseq.meta2_ord$tissue, useNA = "ifany"))

  # --- Subset to Tiss_488B cells only ---
  tissue_vec <- getCellColData(lowseq.proj, select = "tissue", drop = TRUE)
  cells_488B_proj <- getCellNames(lowseq.proj)[!is.na(tissue_vec) & tissue_vec == "Tiss_488B"]
  cat("Subsetting lowseq.proj to Tiss_488B:", length(cells_488B_proj), "cells\n")
  lowseq.proj <- lowseq.proj[cells_488B_proj, ]

  # --- Add LSI, clustering, UMAP ---
  cat("\nAdding Iterative LSI (TileMatrix)...\n")
  lowseq.proj <- addIterativeLSI(
    ArchRProj = lowseq.proj,
    useMatrix = "TileMatrix",
    name = "IterativeTileLSI",
    iterations = 2,
    clusterParams = list(resolution = c(0.2), sampleCells = 10000, n.start = 10),
    varFeatures = 25000,
    dimsToUse = 1:30,
    force = TRUE
  )

  cat("Adding clusters...\n")
  lowseq.proj <- addClusters(
    input = lowseq.proj,
    reducedDims = "IterativeTileLSI",
    method = "Seurat",
    name = "Clusters_tile",
    resolution = 0.5,
    force = TRUE
  )

  cat("Adding UMAP...\n")
  lowseq.proj <- addUMAP(
    ArchRProj = lowseq.proj,
    reducedDims = "IterativeTileLSI",
    name = "UMAP_tile",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine",
    force = TRUE
  )

  cat("Lowseq project built.\n")
  cat("Embeddings:", paste(names(lowseq.proj@embeddings), collapse = ", "), "\n")
  cat("CellColData cols:", paste(colnames(getCellColData(lowseq.proj)), collapse = ", "), "\n")

  saveArchRProject(lowseq.proj, lowseq_save_dir)
  cat("Lowseq project saved to:", lowseq_save_dir, "\n")
}

# =============================================================
# SECTION 2: Load and filter lowseq somatic SNVs
# =============================================================
cat("\n--- Loading lowseq somatic SNV data ---\n")
base_dir <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling"
chromosomes <- paste0("chr", 1:22)

load_and_filter <- function(dataset, chr) {
  f <- file.path(base_dir, dataset, "somatic", paste0(chr, ".putativeSNVs.csv"))
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE)
  d <- d[d$Depth_ref > 5 & d$Depth_alt > 5, ]
  d <- d[d$BAF_alt < 0.5, ]
  d <- d[!is.na(d$LDrefine_merged_score) & d$LDrefine_merged_score > 0.25, ]
  d$snv_id <- paste0(d$chr, ":", d$pos, ":", d$Ref_allele, ":", d$Alt_allele)
  d
}

low_somatic <- do.call(rbind, lapply(chromosomes, function(chr) load_and_filter("lowseq", chr)))
cat("Lowseq filtered somatic SNVs:", nrow(low_somatic), "\n")

# =============================================================
# SECTION 3: Map SNVs to cells using SNV_mat.RDS
# =============================================================
cat("\n--- Mapping SNVs to cells ---\n")

load_snv_mat <- function(chr) {
  f <- file.path(base_dir, "lowseq", "somatic", paste0(chr, ".SNV_mat.RDS"))
  if (!file.exists(f)) return(NULL)
  tryCatch(readRDS(f), error = function(e) { cat("  Error loading", f, ":", e$message, "\n"); NULL })
}

# Get cell names from ArchR project
archr_cells <- getCellNames(lowseq.proj)
# Convert ArchR cell names to barcode format (Lowseq#BARCODE-1 -> BARCODE)
archr_barcodes <- gsub("-1$", "", gsub("^Lowseq#", "", archr_cells))
# Reverse lookup: 16bp barcode -> ArchR cell name
bc16_to_archr <- setNames(archr_cells, archr_barcodes)

# Track mutations per cell
cell_mutation_count <- setNames(rep(0, length(archr_cells)), archr_cells)

for (chr in chromosomes) {
  cat("  Processing", chr, "...\n")
  mat <- load_snv_mat(chr)
  if (is.null(mat)) { cat("    Skipped (no matrix)\n"); next }

  chr_snvs <- low_somatic[grepl(paste0("^", chr, ":"), low_somatic$snv_id), ]
  if (nrow(chr_snvs) == 0) next

  available_snvs <- intersect(chr_snvs$snv_id, rownames(mat))
  if (length(available_snvs) == 0) { cat("    No matching SNVs in matrix\n"); next }

  if (ncol(mat) < 19) { cat("    Matrix too few columns\n"); next }
  cell_cols <- mat[, 19:ncol(mat), drop = FALSE]

  # Lowseq SNV_mat has 16bp barcodes - directly map to ArchR names
  mat_barcodes <- colnames(cell_cols)

  for (snv in available_snvs) {
    snv_row <- cell_cols[snv, , drop = TRUE]
    alt_counts <- as.integer(sub("^[0-9]+\\|", "", snv_row))
    mut_cells_idx <- which(alt_counts > 0)
    if (length(mut_cells_idx) > 0) {
      mut_barcodes <- mat_barcodes[mut_cells_idx]
      archr_mut_names <- bc16_to_archr[mut_barcodes]
      archr_mut_names <- archr_mut_names[!is.na(archr_mut_names)]
      if (length(archr_mut_names) > 0) {
        cell_mutation_count[archr_mut_names] <- cell_mutation_count[archr_mut_names] + 1
      }
    }
  }
  cat("    Processed", length(available_snvs), "SNVs\n")
}

cat("\nMutation count distribution across all cells:\n")
print(summary(cell_mutation_count))
cat("Cells with >=1 mutation:", sum(cell_mutation_count > 0), "/", length(archr_cells), "\n")

# --- Filter to Tiss_488B and cap outliers ---
tissue_tmp <- tryCatch(getCellColData(lowseq.proj, select = "tissue", drop = TRUE),
                       error = function(e) NULL)
if (!is.null(tissue_tmp)) {
  cells_488B <- names(tissue_tmp)[!is.na(tissue_tmp) & tissue_tmp == "Tiss_488B"]
  cells_488B <- cells_488B[cells_488B %in% archr_cells]
} else {
  cells_488B <- archr_cells  # fallback: use all cells if tissue not yet annotated
}

# Determine outlier cap: 99th percentile or 60, whichever is lower
counts_488B <- cell_mutation_count[cells_488B]
if (length(counts_488B) > 0 && any(counts_488B > 0, na.rm = TRUE)) {
  mut_cap <- min(60, quantile(counts_488B, 0.99, na.rm = TRUE))
} else {
  mut_cap <- 60
}
if (!is.finite(mut_cap) || mut_cap < 1) mut_cap <- 60
cat("Mutation count cap:", mut_cap, "\n")
cats_outliers <- sum(counts_488B > mut_cap, na.rm = TRUE)
cat("Cells capped as outliers (>", mut_cap, "):", cats_outliers, "\n")

# Capped mutation count for colour scale only
cell_mut_capped <- pmin(cell_mutation_count, mut_cap)

# Barcode matching diagnostics
cat("\n--- Barcode matching diagnostics ---\n")
sample_mat <- load_snv_mat("chr1")
if (!is.null(sample_mat) && ncol(sample_mat) >= 19) {
  mat_bc <- colnames(sample_mat)[19:ncol(sample_mat)]
  cat("SNV matrix barcodes (sample):", paste(head(mat_bc, 3), collapse = ", "), "\n")
  cat("ArchR barcodes (sample):", paste(head(archr_barcodes, 3), collapse = ", "), "\n")
  cat("SNV matrix barcodes matched in ArchR:", sum(mat_bc %in% archr_barcodes), "/", length(mat_bc), "\n")
}

# =============================================================
# SECTION 4: Add mutation data and generate plots
# =============================================================
cat("\n--- Adding mutation data to ArchR project ---\n")
lowseq.proj <- addCellColData(lowseq.proj, data = cell_mutation_count[archr_cells],
                               name = "n_somatic_mutations", cells = archr_cells, force = TRUE)
has_mut <- ifelse(cell_mutation_count[archr_cells] > 0, "Mutated", "Reference")
lowseq.proj <- addCellColData(lowseq.proj, data = has_mut,
                               name = "has_mutation", cells = archr_cells, force = TRUE)

meta <- as.data.frame(getCellColData(lowseq.proj))
has_spatial <- "x_spatial" %in% colnames(meta) && "y_spatial" %in% colnames(meta)
embedding_names <- names(lowseq.proj@embeddings)
cat("Spatial coords:", has_spatial, "\n")
cat("Embeddings:", paste(embedding_names, collapse = ", "), "\n")

# --- Prepare Tiss_488B meta subsets (capped mutation count) ---
meta$n_mut_capped <- pmin(meta$n_somatic_mutations, mut_cap)
meta_488B <- meta[!is.na(meta$tissue) & meta$tissue == "Tiss_488B", ]
meta_488B <- meta_488B[order(meta_488B$n_mut_capped), ]

# Compute nice legend breaks (0, 1, a few steps up to cap)
legend_breaks <- unique(c(0, 1, round(mut_cap * c(0.25, 0.5, 0.75)), floor(mut_cap)))
legend_labels <- c(as.character(legend_breaks[-length(legend_breaks)]),
                   paste0("≥", floor(mut_cap)))

pdf(file.path(out_dir, "archr_variant_overlay_lowseq.pdf"), width = 14, height = 10)

# =============================================
# PLOT 0: Histogram - cells by # somatic mutations (488B, capped)
# =============================================
cat("  Plot 0: Histogram of somatic mutation counts per cell (488B)\n")
mut_488B <- cell_mutation_count[cells_488B]
df_hist <- data.frame(
  n_mut = mut_488B,
  group = ifelse(mut_488B == 0, "0", ifelse(mut_488B <= mut_cap, "1+", paste0(">", floor(mut_cap))))
)
# Cells with >=1 only histogram
p_hist <- ggplot(df_hist[df_hist$n_mut > 0 & df_hist$n_mut <= mut_cap, ],
                 aes(x = n_mut)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  scale_x_continuous(breaks = seq(0, mut_cap, by = max(1, floor(mut_cap / 10)))) +
  theme_classic() +
  labs(title = "LowSeq Tiss_488B: Cells by # Somatic Mutations (1 to cap)",
       x = "# Somatic Mutations", y = "# Cells") +
  annotate("text", x = mut_cap * 0.7, y = Inf, vjust = 2,
           label = paste0(cats_outliers, " cells capped (>", floor(mut_cap), ")"), size = 3.5)
print(p_hist)

# Bar chart: 0, 1, 2, 3, 4, 5+
df_cat <- data.frame(n_mut = mut_488B)
df_cat$cat <- cut(df_cat$n_mut, breaks = c(-1, 0, 1, 2, 3, 4, Inf),
                  labels = c("0", "1", "2", "3", "4", "5+"))
p_bar <- ggplot(df_cat, aes(x = cat)) +
  geom_bar(fill = c("gray70", "#FDE725", "#5DC863", "#21908C", "#3B528B", "#440154")) +
  theme_classic() +
  labs(title = "LowSeq Tiss_488B: Cell Count by Mutation Category",
       x = "# Somatic Mutations", y = "# Cells")
print(p_bar)

# =============================================
# PLOT 1-2: Spatial plots (488B only, capped scale)
# =============================================
if (has_spatial) {
  cat("  Plot 1: Spatial mutation count (488B, capped)\n")
  p <- ggplot(meta_488B, aes(x = x_spatial, y = y_spatial, color = n_mut_capped)) +
    geom_point(size = 0.7) +
    scale_color_viridis_c(option = "magma",
                          breaks = legend_breaks,
                          labels = legend_labels,
                          limits = c(0, mut_cap)) +
    theme_classic() +
    ggtitle(paste0("LowSeq Tiss_488B: Somatic Mutation Count (capped at ", floor(mut_cap), ")")) +
    labs(color = "# Mutations")
  print(p)

  cat("  Plot 2: Spatial mutation status (488B)\n")
  p <- ggplot(meta_488B, aes(x = x_spatial, y = y_spatial, color = has_mutation)) +
    geom_point(size = 0.7) +
    scale_color_manual(values = c("Reference" = "gray80", "Mutated" = "red"),
                       guide = guide_legend(override.aes = list(size = 3))) +
    theme_classic() +
    ggtitle("LowSeq Tiss_488B: Cells with Somatic Mutations")
  print(p)
}

# =============================================
# PLOT 3-4: UMAP plots (488B cells only, capped scale)
# =============================================
for (emb_name in embedding_names) {
  cat("  Plotting mutations on embedding:", emb_name, "\n")
  tryCatch({
    emb <- getEmbedding(lowseq.proj, embedding = emb_name)
    emb_df <- as.data.frame(emb)
    colnames(emb_df) <- c("dim1", "dim2")
    emb_df$cell <- rownames(emb_df)
    # Filter to 488B
    emb_df <- emb_df[emb_df$cell %in% cells_488B, ]
    emb_df$n_mutations <- cell_mutation_count[emb_df$cell]
    emb_df$n_mut_capped <- pmin(emb_df$n_mutations, mut_cap)
    emb_df$has_mutation <- ifelse(emb_df$n_mutations > 0, "Mutated", "Reference")

    if ("Clusters_tile" %in% colnames(meta)) {
      emb_df$cluster <- meta[emb_df$cell, "Clusters_tile"]
    }

    p <- ggplot(emb_df[order(emb_df$n_mut_capped), ],
                aes(x = dim1, y = dim2, color = n_mut_capped)) +
      geom_point(size = 0.7) +
      scale_color_viridis_c(option = "magma",
                            breaks = legend_breaks,
                            labels = legend_labels,
                            limits = c(0, mut_cap)) +
      theme_classic() +
      ggtitle(paste0("LowSeq Tiss_488B: Somatic Mutations on ", emb_name,
                     " (capped at ", floor(mut_cap), ")")) +
      labs(color = "# Mutations")
    print(p)

    p <- ggplot(emb_df[order(emb_df$n_mutations), ],
                aes(x = dim1, y = dim2, color = has_mutation)) +
      geom_point(size = 0.7) +
      scale_color_manual(values = c("Reference" = "gray80", "Mutated" = "red"),
                         guide = guide_legend(override.aes = list(size = 3))) +
      theme_classic() +
      ggtitle(paste0("LowSeq Tiss_488B: Mutated Cells on ", emb_name))
    print(p)
  }, error = function(e) cat("    Error:", e$message, "\n"))
}

# =============================================
# PLOT 5: Cluster composition with mutations (488B)
# =============================================
if ("Clusters_tile" %in% colnames(meta)) {
  cat("  Plot 5: Cluster mutation composition (488B)\n")
  meta_488B$has_mutation <- ifelse(cell_mutation_count[rownames(meta_488B)] > 0, "Mutated", "Reference")
  cluster_mut <- table(meta_488B$Clusters_tile, meta_488B$has_mutation)
  cluster_mut_pct <- prop.table(cluster_mut, margin = 1) * 100

  if (ncol(cluster_mut_pct) > 0) {
    par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
    barplot(t(cluster_mut), beside = TRUE, col = c("red", "gray80"),
            main = "LowSeq Tiss_488B: Mutation Status by Cluster",
            ylab = "Number of Cells", xlab = "Cluster", las = 1)
    legend("topright", c("Mutated", "Reference"), fill = c("red", "gray80"), cex = 0.8)

    if ("Mutated" %in% colnames(cluster_mut_pct)) {
      barplot(cluster_mut_pct[, "Mutated"], col = "coral",
              main = "LowSeq Tiss_488B: % Cells with Mutations per Cluster",
              ylab = "% Mutated", xlab = "Cluster", las = 1, ylim = c(0, 100))
    }
  }
}

# (Plot 6 removed - histogram now at Plot 0 above)

# =============================================
# PLOT 7: Top individual SNVs - spatial overlay
# =============================================
if (has_spatial) {
  cat("  Plot 7: Top individual SNV spatial plots\n")

  snv_cell_counts <- list()
  for (chr in chromosomes) {
    mat <- load_snv_mat(chr)
    if (is.null(mat) || ncol(mat) < 19) next

    chr_snvs <- low_somatic[grepl(paste0("^", chr, ":"), low_somatic$snv_id), ]
    available_snvs <- intersect(chr_snvs$snv_id, rownames(mat))
    cell_cols <- mat[, 19:ncol(mat), drop = FALSE]

    for (snv in available_snvs) {
      snv_row <- cell_cols[snv, , drop = TRUE]
      alt_vals <- as.integer(sub("^[0-9]+\\|", "", snv_row))
      ref_vals <- as.integer(sub("\\|[0-9]+$", "", snv_row))
      mut_count <- sum(alt_vals > 0, na.rm = TRUE)
      ref_count <- sum(ref_vals > 0 & alt_vals == 0, na.rm = TRUE)
      snv_cell_counts[[snv]] <- c(mut = mut_count, ref = ref_count)
    }
  }

  if (length(snv_cell_counts) > 0) {
    snv_df <- data.frame(
      snv_id = names(snv_cell_counts),
      n_mut = sapply(snv_cell_counts, `[`, "mut"),
      n_ref = sapply(snv_cell_counts, `[`, "ref"),
      stringsAsFactors = FALSE
    )
    snv_df <- snv_df[order(-snv_df$n_mut), ]

    top_snvs <- head(snv_df$snv_id, 6)
    cat("  Top SNVs:", paste(top_snvs, collapse = ", "), "\n")

    for (snv_id in top_snvs) {
      chr <- sub(":.*", "", snv_id)
      mat <- load_snv_mat(chr)
      if (is.null(mat) || ncol(mat) < 19) next
      if (!snv_id %in% rownames(mat)) next

      cell_cols <- mat[, 19:ncol(mat), drop = FALSE]
      snv_row <- cell_cols[snv_id, , drop = TRUE]
      mat_barcodes <- colnames(cell_cols)

      geno <- setNames(rep(0, length(archr_cells)), archr_cells)
      alt_vals <- as.integer(sub("^[0-9]+\\|", "", snv_row))
      ref_vals <- as.integer(sub("\\|[0-9]+$", "", snv_row))
      for (j in seq_along(mat_barcodes)) {
        an <- bc16_to_archr[mat_barcodes[j]]
        if (!is.na(an)) {
          if (!is.na(alt_vals[j]) && alt_vals[j] > 0) geno[an] <- 1
          else if (!is.na(ref_vals[j]) && ref_vals[j] > 0) geno[an] <- -1
        }
      }

      # Subset to 488B only
      meta_488B$snv_geno <- geno[rownames(meta_488B)]
      meta_488B$snv_status <- ifelse(meta_488B$snv_geno == 1, "Alt",
                                     ifelse(meta_488B$snv_geno == -1, "Ref", "No Data"))
      n_alt_488B <- sum(meta_488B$snv_status == "Alt")

      p <- ggplot(meta_488B[order(abs(meta_488B$snv_geno)), ],
                  aes(x = x_spatial, y = y_spatial, color = snv_status)) +
        geom_point(size = 0.6) +
        scale_color_manual(values = c("No Data" = "gray90", "Ref" = "steelblue", "Alt" = "red"),
                           guide = guide_legend(override.aes = list(size = 3))) +
        theme_classic() +
        ggtitle(paste0("LowSeq Tiss_488B: ", snv_id, "  (", n_alt_488B, " Alt cells)")) +
        theme(plot.title = element_text(size = 9))
      print(p)
    }
  }
}

dev.off()
cat("\nPlots saved to:", file.path(out_dir, "archr_variant_overlay_lowseq.pdf"), "\n")

cat("\nEnd time:", format(Sys.time()), "\n")
cat("=== Done ===\n")
