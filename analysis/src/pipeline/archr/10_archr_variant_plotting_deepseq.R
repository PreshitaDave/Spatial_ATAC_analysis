#!/usr/bin/env Rscript
# 10_archr_variant_plotting.R
# Plot somatic variants on ArchR spatial + UMAP embeddings
# Uses deepseq ArchR project (tiss_488B) and somatic SNV data from Monopogen

cat("=== ArchR Somatic Variant Plotting ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# --- Setup ---
.libPaths(c('/projectnb/paxlab/presh/env/R_4.4/ArchR_libs', .libPaths()))

library(ArchR)
library(ggplot2)
library(parallel)
library(Matrix)

set.seed(1)
addArchRGenome("hg38")
addArchRThreads(threads = as.integer(Sys.getenv("NSLOTS", "2")))

result_dir <- Sys.getenv(
  "RESULT_DIR",
  unset = "/projectnb/paxlab/presh/projects/spatial_atac/Data/05_results/variant_calling/somatic_comparison/tables"
)
plot_dir <- Sys.getenv(
  "PLOT_DIR",
  unset = "/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/comparison/somatic"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load ArchR project ----
cat("Loading saved ArchR project...\n")
proj <- loadArchRProject("/projectnb/paxlab/presh/projects/spatial_atac/Data/")
cat("Loaded", ncol(proj), "... no,", nrow(getCellColData(proj)), "cells\n")

# ---- Adding spatial locs ----

spatial_locs <- read.csv(file.path('/projectnb/paxlab/yeting/SpatialATACseq/data',
                                   "/tissue_positions_list.csv"),
                         header = FALSE)
colnames(spatial_locs) <- c("barcode", "in_tissue", "array_row", "array_col", "x_spatial", "y_spatial")

xy <- spatial_locs[, c("barcode", "x_spatial", "y_spatial")]
xy2 = spatial_locs[spatial_locs$in_tissue == 1, ]
rownames(xy2) = paste0("Deepseq#", xy2$barcode, "-1")
#subset to archR project cells in the same order 
xy3 = xy2[rownames(proj@cellColData), ]

new.meta = cbind(proj@cellColData, xy3)

#Update the archr project metdata
proj@cellColData = new.meta


# Check available embeddings and columns
cat("Available embeddings:", paste(names(proj@embeddings), collapse = ", "), "\n")
cat("Available cellColData:", paste(head(colnames(getCellColData(proj)), 50), collapse = ", "), "\n")

# ---- Load barcode mapping (8bp <-> 16bp) ----
cat("\nLoading barcode mapping...\n")
cell_data_16bp <- read.csv("/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq_cell_data.csv",
                           stringsAsFactors = FALSE)
cell_data_8bp  <- read.csv("/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq_cell_data_8bp.csv",
                           stringsAsFactors = FALSE)
# Map by id: 8bp barcode -> 16bp barcode
bc_map <- setNames(cell_data_16bp$cell, cell_data_8bp$cell)  # 8bp -> 16bp
cat("Barcode mapping loaded:", length(bc_map), "entries\n")

# Get cell names from ArchR project
archr_cells <- getCellNames(proj)
# Convert ArchR cell names to barcode format (Deepseq#BARCODE-1 -> BARCODE)
archr_barcodes <- gsub("-1$", "", gsub("^Deepseq#", "", archr_cells))
# Create reverse lookup: 16bp barcode -> ArchR cell name
bc16_to_archr <- setNames(archr_cells, archr_barcodes)



# ---- Load somatic SNV csv data ----
cat("\nLoading somatic SNV data...\n")
base_dir <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling"
chromosomes <- paste0("chr", 1:22)

# Load and filter following Monopogen pipeline
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

deep_somatic <- do.call(rbind, lapply(chromosomes, function(chr) load_and_filter("deepseq", chr)))
cat("Deepseq filtered somatic SNVs:", nrow(deep_somatic), "\n")

# ---- Load SNV-cell matrices ----
cat("\nLoading SNV-cell matrices for top chromosomes...\n")

# Determine which chromosomes have the most filtered variants
chr_counts <- table(deep_somatic$chr_name <- sub(":.*", "", deep_somatic$snv_id))

# Load SNV_mat.RDS for each chromosome
load_snv_mat <- function(chr) {
  f <- file.path(base_dir, "deepseq", "somatic", paste0(chr, ".SNV_mat.RDS"))
  if (!file.exists(f)) return(NULL)
  tryCatch(readRDS(f), error = function(e) { cat("  Error loading", f, ":", e$message, "\n"); NULL })
}

# Process each chromosome to get mutation profiles per cell
cat("Processing SNV matrices...\n")


all_mats_list <- list()
for (chr in chromosomes) {
  cat("  Processing", chr, "...\n")
  mat <- load_snv_mat(chr)
  if (is.null(mat)) { cat("    Skipped (no matrix)\n"); next }
  all_mats_list[[chr]] <- mat
  
}

big_snv_matrix <- do.call(rbind, unname(all_mats_list))

# subset to deep_somatic snvs 
big_snv_matrix_sub <- big_snv_matrix[deep_somatic$snv_id, ]

saveRDS(big_snv_matrix_sub, file.path(result_dir, "all_chr_snv_mat.rds"))


# ---- Encode SNV-cell matrices ----
encode_chrs <- function(mat) {
  vec <- as.character(as.matrix(mat))
  res <- rep(0, length(vec))
  res[vec == "-1"] <- NA
  res[grepl("^[1-9][0-9]*\\|0$", vec)] <- -1
  res[grepl("^0\\|[1-9][0-9]*$", vec)] <- 1
  
  out <- matrix(res, nrow = nrow(mat), ncol = ncol(mat))
  rownames(out) <- rownames(mat)
  colnames(out) <- colnames(mat)
  return(out)
}


cell_cols <- big_snv_matrix_sub[, 19:ncol(big_snv_matrix_sub), drop = FALSE]
big_snv_matrix_encode <- encode_chrs(cell_cols)

saveRDS(big_snv_matrix_encode, file.path(result_dir, "all_chr_snv_mat_encode.rds"))

# SANITY CHECK - 

# Wasn't sure about the values in the SNV_mat rds and the somatic csv files, what were the Ref and alt depths adding up to? 
# Turns out rds files add up to teh csv file values, but not within the mat fule before cell_counts was made

# process_mat <- as.matrix(cell_cols)
# process_mat[process_mat == "-1"] <- "0|0"
# 
# # 2. Extract Haplotype 1 (left of |) 
# # sub("\\|.*", "", x) removes everything from the pipe onwards
# hap1_counts <- matrix(as.numeric(sub("\\|.*", "", process_mat)), 
#                       nrow = nrow(process_mat))
# 
# # 3. Extract Haplotype 2 (right of |)
# # sub(".*\\|", "", x) removes everything up to and including the pipe
# hap2_counts <- matrix(as.numeric(sub(".*\\|", "", process_mat)), 
#                       nrow = nrow(process_mat))
# 
# # 4. Create the final data frame with row-wise sums
# pseudo_bulk_df <- data.frame(
#   Hap1_Total = rowSums(hap1_counts, na.rm = TRUE),
#   Hap2_Total = rowSums(hap2_counts, na.rm = TRUE),
#   row.names  = rownames(cell_cols)
# )
# 
# # View the result
# head(pseudo_bulk_df)

#---- Investigate the mutation distribution among cell types/clusters ----



library(dplyr)
library(stringr)
library(progress)

# Inputs you must have:
# big_snv_matrix_encode : numeric matrix (rows=SNVs, cols=cells), values -1,0,1,NA
# proj            : ArchRProject object
# celltype_colname      : name of the cell type column in proj@cellColData (default "CellType")
# out_dir               : output directory
# min_counts_per_celltype: minimum informative counts per cell type to test (default 3)

celltype_colname <- "Clusters_tile"   # change if your column is named differently
min_counts_per_celltype <- 3
p_adjust_method <- "BH"


big_snv_matrix_encode <- readRDS(file.path(result_dir, "all_chr_snv_mat_encode.rds"))

# Basic checks
if (!exists("big_snv_matrix_encode")) stop("big_snv_matrix_encode not found.")
if (!exists("proj")) stop("proj not found.")
if (!(celltype_colname %in% colnames(proj@cellColData))) {
  stop("Cell type column '", celltype_colname, "' not found in proj@cellColData.")
}

encoded_mat <- big_snv_matrix_encode
colnames(encoded_mat) <- paste0("Deepseq#", colnames(encoded_mat), '-1')
snv_ids <- rownames(encoded_mat)
matrix_barcodes <- colnames(encoded_mat)

# Get ArchR cell barcodes
archr_cells <- rownames(proj@cellColData)

# Align barcodes: prefer exact match; if not, try to find intersection
if (!all(matrix_barcodes %in% archr_cells)) {
  warning("Not all matrix barcodes are found in ArchR cellColData. Trying intersection.")
  common_barcodes <- intersect(matrix_barcodes, archr_cells)
  if (length(common_barcodes) == 0) stop("No overlapping barcodes between matrix and ArchR project.")
  # Subset matrix to common barcodes and reorder archr cellColData accordingly
  encoded_mat <- encoded_mat[, common_barcodes, drop = FALSE]
  message("Subsetting matrix to ", ncol(encoded_mat), " common barcodes.")
  # reorder ArchR metadata to match matrix columns
  archr_meta <- as.data.frame(proj@cellColData[common_barcodes, , drop = FALSE])
} else {
  # exact match, reorder metadata to matrix column order
  archr_meta <- as.data.frame(proj@cellColData[matrix_barcodes, , drop = FALSE])
}

# Extract cell types
cell_types <- as.character(archr_meta[[celltype_colname]])
unique_celltypes <- sort(unique(cell_types))
message("Found ", length(unique_celltypes), " cell types: ", paste(unique_celltypes, collapse = ", "))

# Precompute masks for -1 and +1
is_neg1 <- (encoded_mat == -1)
is_pos1 <- (encoded_mat == 1)

total_neg1_per_snv <- rowSums(is_neg1, na.rm = TRUE)
total_pos1_per_snv <- rowSums(is_pos1, na.rm = TRUE)

# Fisher helper
fisher_for_celltype <- function(obs_pos, obs_neg, tot_pos, tot_neg) {
  mat <- matrix(c(obs_pos, obs_neg, tot_pos - obs_pos, tot_neg - obs_neg),
                nrow = 2, byrow = FALSE)
  if (any(mat < 0) || sum(mat) == 0) return(NA_real_)
  ft <- try(fisher.test(mat), silent = TRUE)
  if (inherits(ft, "try-error")) return(NA_real_)
  return(ft$p.value)
}

# Loop through SNVs
results_list <- vector("list", nrow(encoded_mat))
names(results_list) <- snv_ids
pb <- progress_bar$new(total = nrow(encoded_mat), format = "[:bar] :current/:total :percent eta: :eta")

for (i in seq_len(nrow(encoded_mat))) {
  pb$tick()
  snv <- snv_ids[i]
  tot_neg <- total_neg1_per_snv[i]
  tot_pos <- total_pos1_per_snv[i]
  tot_info <- tot_neg + tot_pos
  if (is.na(tot_info) || tot_info < 3) {
    results_list[[i]] <- NULL
    next
  }
  
  per_ct_results <- data.frame(celltype=character(0),
                               obs_neg=integer(0),
                               obs_pos=integer(0),
                               pval=numeric(0),
                               stringsAsFactors = FALSE)
  for (ct in unique_celltypes) {
    cols_ct <- which(cell_types == ct)
    if (length(cols_ct) == 0) next
    obs_neg_ct <- sum(is_neg1[i, cols_ct], na.rm = TRUE)
    obs_pos_ct <- sum(is_pos1[i, cols_ct], na.rm = TRUE)
    if ((obs_neg_ct + obs_pos_ct) < min_counts_per_celltype) next
    pval <- fisher_for_celltype(obs_pos = obs_pos_ct, obs_neg = obs_neg_ct, tot_pos = tot_pos, tot_neg = tot_neg)
    per_ct_results <- rbind(per_ct_results, data.frame(celltype = ct,
                                                       obs_neg = obs_neg_ct,
                                                       obs_pos = obs_pos_ct,
                                                       pval = pval,
                                                       stringsAsFactors = FALSE))
  }
  if (nrow(per_ct_results) == 0) {
    results_list[[i]] <- NULL
    next
  }
  per_ct_results <- per_ct_results %>% mutate(pval_adj = p.adjust(pval, method = p_adjust_method))
  best_idx <- which.min(per_ct_results$pval_adj)
  best_row <- per_ct_results[best_idx, , drop = FALSE]
  results_list[[i]] <- data.frame(SNV = snv,
                                  celltype = best_row$celltype,
                                  obs_neg = best_row$obs_neg,
                                  obs_pos = best_row$obs_pos,
                                  pval = best_row$pval,
                                  pval_adj = best_row$pval_adj,
                                  tot_neg = tot_neg,
                                  tot_pos = tot_pos,
                                  stringsAsFactors = FALSE)
}

res_df <- do.call(rbind, results_list)
if (is.null(res_df) || nrow(res_df) == 0) {
  message("No SNVs passed testing criteria.")
  res_df <- data.frame()
} else {
  res_df <- res_df %>% arrange(pval_adj)
  rownames(res_df) <- NULL
}

message("Tested SNVs: ", length(results_list))
message("Significant hits (pval_adj < 0.05): ", sum(res_df$pval_adj < 0.05, na.rm = TRUE))

saveRDS(res_df, file.path(result_dir, "deepseq_snv_fisher_test.rds"))

# Result columns
# 
# SNV: SNV identifier (rownames of encoded matrix)
# celltype: the cell type with the smallest adjusted p-value (best candidate enrichment)
# obs_neg: number of −1 calls in that cell type for this SNV
# obs_pos: number of +1 calls in that cell type
# pval: raw Fisher p-value for that cell type vs rest
# pval_adj: p-value adjusted across the cell types tested for this SNV (BH)
# tot_neg: total −1 across all cells for this SNV
# tot_pos: total +1 across all cells for this SNV



# ---- Generate plots ----
cat("\n--- Generating variant overlay plots ---\n")

# Get cell metadata with spatial coords
meta <- as.data.frame(getCellColData(proj))

# Check if spatial coords exist
has_spatial <- "x_spatial" %in% colnames(meta) && "y_spatial" %in% colnames(meta)
if (has_spatial) cat("Spatial coordinates available.\n") else cat("WARNING: No spatial coordinates found.\n")

# Available embeddings
embedding_names <- names(proj@embeddings)
cat("Available embeddings:", paste(embedding_names, collapse = ", "), "\n")

pdf(file.path(plot_dir, "archr_variant_overlay.pdf"), width = 14, height = 10)




# -------------------------
# Parameters (adjust if needed)
# -------------------------
cell_col_start     <- 19       # column index where cell columns start in raw per-chr matrices
top_n_snvs         <- 6        # number of SNVs to show spatially
min_cells_show     <- 1        # minimum number of ALT-carrying cells to include SNV
cell_barcode_type  <- "archr"  # "archr" or "short"

# -------------------------
# Basic checks
# -------------------------
if (!exists("proj")) stop("proj not found. Load your ArchRProject as `proj`.")
if (!"Clusters_tile" %in% colnames(proj@cellColData)) {
  stop("proj@cellColData must contain a 'Clusters_tile' column. Rename or create before running.")
}
archr_cells <- rownames(proj@cellColData)
if (is.null(archr_cells) || length(archr_cells) == 0) stop("No ArchR cell barcodes found in proj@cellColData rownames.")

# -------------------------
# Build per-cell mutation counts and has_mutation
# -------------------------

# At this point big_snv_matrix_encode should exist
big_snv_matrix_encode <- readRDS(file.path(result_dir, "all_chr_snv_mat_encode.rds"))

if (!exists("big_snv_matrix_encode")) stop("big_snv_matrix_encode not available after attempted build.")

# Compute per-cell alt counts and has_mutation depending on barcode type
if (cell_barcode_type == "archr") {
  # columns are ArchR barcodes
  if (!all(colnames(big_snv_matrix_encode) %in% archr_cells)) {
    warning("Not all matrix columns match ArchR barcodes. Will intersect and align.")
    colnames(big_snv_matrix_encode) = paste0("Deepseq#",colnames(big_snv_matrix_encode),"-1",sep = '')
    common <- intersect(colnames(big_snv_matrix_encode), archr_cells)
    if (length(common) == 0) stop("No overlapping barcodes between encoded matrix and ArchR.")
    mat_sub <- big_snv_matrix_encode[, common, drop = FALSE]
  } else {
    mat_sub <- big_snv_matrix_encode[, archr_cells, drop = FALSE]  # reorder to archr order
  }
  # per-cell mutation counts: number of ALT calls (1)
  cell_alt_count <- colSums(mat_sub == 1, na.rm = TRUE)
  # per-cell informative counts (non-zero)
  cell_inform_count <- colSums(!is.na(mat_sub) & mat_sub != 0)
  # has_mutation flag
  has_mutation_vec <- cell_alt_count > 0
  # attach to ArchR metadata (safe assignment)
  proj@cellColData$has_mutation <- FALSE
  proj@cellColData$cell_mutation_count <- 0L
  proj@cellColData[names(cell_alt_count), "has_mutation"] <- has_mutation_vec
  proj@cellColData[names(cell_alt_count), "cell_mutation_count"] <- cell_alt_count
  # Also make R variables for plotting code
  cell_mutation_count <- cell_alt_count[rownames(proj@cellColData)]
} 


cat("Attached has_mutation and cell_mutation_count to ArchR metadata.\n")

# Create a working copy of metadata for plotting (do not mutate proj in place further)
meta_plot <- as.data.frame(proj@cellColData)
meta_plot$CellType <- as.character(meta_plot$CellType)  # ensure char

# -------------------------
# PLOT 5: Cluster mutation composition
# -------------------------
if ("Clusters_tile" %in% colnames(meta_plot)) {
  cat("Plot 5: Cluster mutation composition\n")
  if (!"has_mutation" %in% colnames(meta_plot)) {
    warning("meta_plot does not contain 'has_mutation'; skipping Plot 5")
  } else {
    cluster_mut <- table(meta_plot$Clusters_tile, meta_plot$has_mutation)
    cluster_mut_pct <- prop.table(cluster_mut, margin = 1) * 100
    
    if (ncol(cluster_mut_pct) > 0) {
      op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
      par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
      
      # Absolute counts
      barplot(t(cluster_mut), beside = TRUE, col = c("red", "gray80"),
              main = "Mutation Status by Cluster",
              ylab = "Number of Cells", xlab = "Cluster", las = 1)
      legend("topright", c("Mutated", "Reference"), fill = c("red", "gray80"), cex = 0.8)
      
      # Percentage mutated (try to find the TRUE/Mutated column)
      pct_col <- if ("TRUE" %in% colnames(cluster_mut_pct)) "TRUE" else if ("Mutated" %in% colnames(cluster_mut_pct)) "Mutated" else colnames(cluster_mut_pct)[1]
      barplot(cluster_mut_pct[, pct_col], col = "coral",
              main = "% Cells with Mutations per Cluster",
              ylab = "% Mutated", xlab = "Cluster", las = 1, ylim = c(0, 100))
    }
  }
} else {
  message("Clusters_tile not found in metadata; skipping Plot 5.")
}

# -------------------------
# PLOT 6: Mutation frequency distribution
# -------------------------
cat("Plot 6: Mutation frequency distribution\n")
op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
par(mfrow = c(1, 1), mar = c(5, 5, 3, 1))

if (!exists("cell_mutation_count")) {
  warning("cell_mutation_count not found — skipping Plot 6")
} else {
  mut_cells_only <- cell_mutation_count[cell_mutation_count > 0]
  if (length(mut_cells_only) > 0) {
    max_val <- max(mut_cells_only, na.rm = TRUE)
    breaks_n <- min(max(50, max_val), 200)
    hist(mut_cells_only, breaks = breaks_n,
         col = "steelblue",
         main = "Distribution of Mutations per Cell\n(Cells with >=1 mutation)",
         xlab = "Number of Somatic Mutations", ylab = "Number of Cells")
  } else {
    message("No cells with mutations to plot for histogram.")
  }
}

# -------------------------
# PLOT 7: Top mutated SNVs spatial overlay
# -------------------------
# ---- pick top SNVs ----
# choose ranking field; here use obs_pos (cells with ALT) — change to n_mut or pval if you prefer
results_table <- res_df[order(res_df$pval_adj),]
results_table <- results_table[results_table$celltype == 'C2',]
n_top = 6

if (!"obs_pos" %in% colnames(results_table)) stop("results_table must contain obs_pos column")
top_snvs <- head(results_table %>% arrange(desc(obs_pos)) %>% pull(SNV), n_top)
message("Top SNVs: ", paste(top_snvs, collapse = ", "))

# ---- helper to get per-cell encoded status for a given SNV ----
get_snv_genotype <- function(snv_id) {
  # returns a named vector of length archr_cells with values: "Alt","Ref","No Data"
  geno_vec <- setNames(rep("No Data", length(archr_cells)), archr_cells)
  
  # prefer encoded matrix if present and contains this SNV
  if (exists("big_snv_matrix_encode") && snv_id %in% rownames(big_snv_matrix_encode)) {
    row_vals <- big_snv_matrix_encode[snv_id, , drop = TRUE]  # numeric -1,0,1,NA
    if (cell_barcode_type == "archr") {
      # columns should be archr barcodes; align by names
      common <- intersect(names(row_vals), archr_cells)
      if (length(common) == 0) stop("No overlapping barcodes between encoded matrix and ArchR")
      geno_vec[common][row_vals[common] == 1] <- "Alt"
      geno_vec[common][row_vals[common] == -1] <- "Ref"
    } else {
      # short barcodes -> need mapping
      if (!exists("bc_map") || !exists("bc16_to_archr")) stop("bc_map and bc16_to_archr required for short barcode mode")
      short_bcs <- names(row_vals)
      bc16 <- bc_map[short_bcs]
      archr_ids <- bc16_to_archr[bc16]
      mapped_idx <- which(!is.na(archr_ids) & archr_ids %in% archr_cells)
      if (length(mapped_idx) > 0) {
        for (j in mapped_idx) {
          a <- archr_ids[j]
          v <- row_vals[j]
          if (is.na(v)) next
          if (v == 1) geno_vec[a] <- "Alt"
          else if (v == -1) geno_vec[a] <- "Ref"
        }
      }
    }
    return(geno_vec)
  }
  
  # fallback: use per-chromosome raw matrix parsing
  chr <- sub(":.*", "", snv_id)
  if (!exists("load_snv_mat")) stop("load_snv_mat function not found for fallback")
  mat <- load_snv_mat(chr)
  if (is.null(mat) || !snv_id %in% rownames(mat)) {
    warning("SNV ", snv_id, " not found in matrices; returning No Data.")
    return(geno_vec)
  }
  cell_cols <- mat[, 19:ncol(mat), drop = FALSE]
  mat_barcodes <- colnames(cell_cols)
  # parse alt/ref counts
  alt_vals <- suppressWarnings(as.integer(sub("^[^\\|]*\\|", "", cell_cols[snv_id, ])))
  ref_vals <- suppressWarnings(as.integer(sub("\\|[^\\|]*$", "", cell_cols[snv_id, ])))
  if (cell_barcode_type == "archr") {
    # assume mat_barcodes are archr barcodes
    common <- intersect(mat_barcodes, archr_cells)
    idxs <- match(common, mat_barcodes)
    for (k in seq_along(common)) {
      b <- common[k]; j <- idxs[k]
      a_val <- alt_vals[j]; r_val <- ref_vals[j]
      if (!is.na(a_val) && a_val > 0) geno_vec[b] <- "Alt"
      else if (!is.na(r_val) && r_val > 0) geno_vec[b] <- "Ref"
    }
  } else {
    # short barcode mapping
    if (!exists("bc_map") || !exists("bc16_to_archr")) stop("bc_map and bc16_to_archr required for short barcode mode")
    bc16 <- bc_map[mat_barcodes]
    archr_ids <- bc16_to_archr[bc16]
    for (j in seq_along(mat_barcodes)) {
      archr_id <- archr_ids[j]
      if (is.na(archr_id) || !(archr_id %in% archr_cells)) next
      a_val <- alt_vals[j]; r_val <- ref_vals[j]
      if (!is.na(a_val) && a_val > 0) geno_vec[archr_id] <- "Alt"
      else if (!is.na(r_val) && r_val > 0) geno_vec[archr_id] <- "Ref"
    }
  }
  return(geno_vec)
}

# ---- plotting loop: build one ggplot per SNV and arrange ----
plots_list <- vector("list", length(top_snvs))
for (i in seq_along(top_snvs)) {
  snv <- top_snvs[i]
  geno_vec <- get_snv_genotype(snv)
  plot_df <- archr_meta %>%
    mutate(SNV_status = factor(geno_vec[rownames(archr_meta)], levels = c("No Data", "Ref", "Alt")))
  p <- ggplot(plot_df, aes(x = x_spatial, y = y_spatial, color = SNV_status)) +
    geom_point(size = 0.6) +
    scale_color_manual(values = c("No Data" = "gray90", "Ref" = "blue", "Alt" = "red")) +
    theme_classic() +
    ggtitle(snv) +
    theme(plot.title = element_text(size = 10),
          legend.position = "right",
          legend.title = element_blank())
  plots_list[[i]] <- p
}

# Arrange plots in grid (2x3)
# do.call(grid.arrange, c(plots_list, ncol = 2))



out_pdf <- file.path(plot_dir, "top_snvs_spatial.pdf")
pdf(out_pdf, width = 12, height = 6)   # open PDF device

for (p in plots_list) {
  print(p)    # each print() writes a new page
}

dev.off()     # close device
message("Saved PDF: ", out_pdf)


# -------------------------
# PLOT 8: Top mutated SNVs HEATMAP 
# -------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)  # optional; we use ggplot below

# CONFIG
results_table <- res_df         # your results table (SNV, celltype, obs_pos, ...)
clusters <- c("C1","C2","C3","C4")
top_n <- 25
out_pdf <- "snv_cluster_heatmap_top25_per_cluster.pdf"

# sanity checks
stopifnot(exists("big_snv_matrix_encode"))
stopifnot(is.matrix(big_snv_matrix_encode) || is.data.frame(big_snv_matrix_encode))
stopifnot(all(colnames(big_snv_matrix_encode) %in% rownames(proj@cellColData)))

# 1) pick top N SNVs per cluster by obs_pos
top_snvs_by_cluster <- res_df %>%
  filter(celltype %in% clusters) %>%
  group_by(celltype) %>%
  arrange(desc(obs_pos)) %>%
  slice_head(n = top_n) %>%
  ungroup() %>%
  pull(SNV) %>%
  unique()

length(top_snvs_by_cluster)  # up to 4*25 = 100 SNVs (fewer if overlap)

# 2) build per-cluster percent ALT matrix for those SNVs
# matrix rows = SNV, cols = clusters, values = fraction_of_cells_with_ALT (0-1)
mat <- big_snv_matrix_encode  # numeric matrix: rows=SNV, cols=ArchR cell barcodes

# ensure order of columns matches archr metadata
archr_meta <- as.data.frame(proj@cellColData)
mat <- mat[, rownames(archr_meta), drop = FALSE]  # reorder columns to ArchR order

# function to compute fraction ALT per cluster for a vector of SNVs
compute_cluster_frac <- function(snv_ids, mat, meta, clusters) {
  res_list <- list()
  for (snv in snv_ids) {
    if (!snv %in% rownames(mat)) {
      # if missing, return zeros
      res_list[[snv]] <- setNames(rep(0, length(clusters)), clusters)
      next
    }
    row_vals <- mat[snv, , drop = TRUE]  # per-cell values: -1,0,1,NA
    # ALT presence defined as ==1
    per_cluster_frac <- sapply(clusters, function(cl) {
      cells <- which(as.character(meta$Clusters_tile) == cl)
      if (length(cells) == 0) return(NA_real_)
      alt_count <- sum(row_vals[cells] == 1, na.rm = TRUE)
      # use denominator = number of cells in that cluster (or number of covered cells if you prefer)
      denom <- length(cells)
      frac <- alt_count / denom
      return(frac)
    })
    res_list[[snv]] <- per_cluster_frac
  }
  mat_out <- do.call(rbind, res_list)
  rownames(mat_out) <- snv_ids
  colnames(mat_out) <- clusters
  return(mat_out)
}

frac_mat <- compute_cluster_frac(top_snvs_by_cluster, mat, archr_meta, clusters)

# 3) reorder rows by which cluster has highest value (for triangular diagonal as example)
# We'll compute cluster of maximum fraction per SNV and sort by cluster then by fraction descending
# Debug: check lengths and NAs
max_cluster <- apply(frac_mat, 1, function(x) {
  idx <- which.max(x)
  if (length(idx) == 0 || is.na(idx)) NA else clusters[idx]
})
max_frac <- apply(frac_mat, 1, function(x) max(x, na.rm = TRUE))

# Check for problematic rows
cat("max_cluster length:", length(max_cluster), "\n")
cat("max_frac length:", length(max_frac), "\n")
cat("NAs in max_cluster:", sum(is.na(max_cluster)), "\n")
cat("NAs in max_frac:", sum(is.na(max_frac)), "\n")

# If there are rows with all NAs, remove them
valid_rows <- !is.na(max_cluster) & !is.na(max_frac)
if (!all(valid_rows)) {
  message("Removing ", sum(!valid_rows), " rows with all NAs")
  frac_mat <- frac_mat[valid_rows, , drop = FALSE]
  max_cluster <- max_cluster[valid_rows]
  max_frac <- max_frac[valid_rows]
}

# Now order safely (convert clusters to numeric for sorting)
cluster_order <- as.numeric(factor(max_cluster, levels = clusters))
order_idx <- order(cluster_order, -max_frac)
frac_mat_ord <- frac_mat[order_idx, , drop = FALSE]

cat("Final matrix dimensions:", dim(frac_mat_ord), "\n")

# Optionally scale rows (e.g., by max) or not; example uses raw fractions (0-1)
# Make annotation for rows: show SNV names on right side — we'll use pheatmap and adjust

# 4) Plot heatmap using pheatmap: rows = SNVs, columns = clusters
library(pheatmap)

# create a color palette (white -> pale red -> dark red)
library(RColorBrewer)
cols <- colorRampPalette(c("white", "#f7b0b0", "#d7301f"))(50)

pdf(file.path(plot_dir, out_pdf), width = 6, height = max(4, nrow(frac_mat_ord) * 0.12 + 1))
pheatmap(frac_mat_ord,
         color = cols,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 10,
         border_color = "black",
         cellwidth = 40,
         cellheight = 8,
         legend = TRUE,
         main = paste0("Top SNVs per cluster (top ", top_n, " by obs_pos)"),
         display_numbers = FALSE,
         angle_col = 45,
         na_col = "gray95",
         breaks = seq(0, 1, length.out = length(cols)+1)
)

dev.off()
message("Saved heatmap: ", out_pdf)






library(pheatmap)
library(RColorBrewer)

cols <- colorRampPalette(c("white", "#f7b0b0", "#d7301f"))(50)

pdf(file.path(plot_dir, out_pdf), width = 6, height = max(4, nrow(frac_mat_ord) * 0.12 + 1))

pheatmap(frac_mat_ord,
         color = cols,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 10,
         # border_color = "black",    # REMOVE THIS LINE
         cellwidth = 40,
         cellheight = 8,
         legend = TRUE,
         main = paste0("Top SNVs per cluster (top ", top_n, " by obs_pos)")
)

dev.off()
message("Saved heatmap: ", out_pdf)


cat("\nEnd time:", format(Sys.time()), "\n")
cat("=== Done ===\n")



