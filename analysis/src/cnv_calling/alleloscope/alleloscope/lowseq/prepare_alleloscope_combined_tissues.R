#!/usr/bin/env Rscript
# =============================================================================
# Prepare combined lowseq Alleloscope inputs for tissues 488B + 489
# 
# Combines:
#   1. SNP matrices (alt_all.mtx, ref_all.mtx) — union variants, cbind cells
#   2. Barcodes (concatenate tissue-specific barcodes)
#   3. Variant VCF (union of variants)
#   4. Fragment counts (cbind across tissues)
#   5. Segmentation table (use from 488B - no need to regenerate)
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
output_dir <- file.path(project_root, "Data/04_analysis/cnv/alleloscope/lowseq_tissue_from_existing/combined_488B_489")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║  Alleloscope Combined Lowseq Prep (488B + 489)                ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
message(sprintf("[%s] START: Preparing combined Alleloscope inputs", format(Sys.time(), "%H:%M:%S")))
message(sprintf("[%s] Output: %s", format(Sys.time(), "%H:%M:%S"), output_dir))
cat("\n")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Load tissue-specific data
# ─────────────────────────────────────────────────────────────────────────────
cat("█ STEP 1: Loading tissue-specific data\n")
message(sprintf("[%s] ├─ Loading 488B matrices and metadata...", format(Sys.time(), "%H:%M:%S")))

tissues <- c("488B", "489")
tissue_data <- list()

for (tissue in tissues) {
  tissue_dir <- file.path(project_root, "Data/04_analysis/cnv/alleloscope/lowseq_tissue_from_existing", tissue)
  
  # Load matrices
  alt_mat <- readMM(file.path(tissue_dir, "alt_all.mtx"))
  ref_mat <- readMM(file.path(tissue_dir, "ref_all.mtx"))
  
  # Load barcodes
  barcodes <- readLines(file.path(tissue_dir, "barcodes.tsv"))
  barcodes <- barcodes[nzchar(barcodes)]
  colnames(alt_mat) <- barcodes
  colnames(ref_mat) <- barcodes
  
  # Load VCF
  vcf_df <- read.table(
    file.path(tissue_dir, "var_all.vcf"),
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "#",
    skip = 1
  )
  
  # Load raw counts
  raw_counts <- read.table(
    file.path(tissue_dir, "chr1000k_fragments.tsv"),
    sep = "\t",
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
  )
  colnames(raw_counts) <- gsub("[.]", "-", colnames(raw_counts))
  
  # Load segmentation table (just from 488B - no need to regenerate)
  seg_table <- readRDS(file.path(tissue_dir, "seg_table.rds"))
  
  tissue_data[[tissue]] <- list(
    alt_mat = alt_mat,
    ref_mat = ref_mat,
    barcodes = barcodes,
    vcf_df = vcf_df,
    raw_counts = raw_counts,
    seg_table = seg_table
  )
  
  message(sprintf(
    "[%s] ├─ %s: ✓ %d variants × %d cells | Fragments: %d bins",
    format(Sys.time(), "%H:%M:%S"), tissue,
    nrow(alt_mat), ncol(alt_mat),
    nrow(raw_counts)
  ))
}

message(sprintf("[%s] └─ Data loading: COMPLETE\n", format(Sys.time(), "%H:%M:%S")))

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Combine SNP matrices (union variants, cbind cells)
# ─────────────────────────────────────────────────────────────────────────────
cat("█ STEP 2: Combining SNP matrices\n")
message(sprintf("[%s] ├─ Creating variant union...", format(Sys.time(), "%H:%M:%S")))

vcf_488B <- tissue_data[["488B"]]$vcf_df
vcf_489 <- tissue_data[["489"]]$vcf_df

vcf_488B$variant_key <- with(vcf_488B, paste0(V1, ":", V2, ":", V4, ":", V5))
vcf_489$variant_key <- with(vcf_489, paste0(V1, ":", V2, ":", V4, ":", V5))

all_variant_keys <- unique(c(vcf_488B$variant_key, vcf_489$variant_key))
n_variants_union <- length(all_variant_keys)

idx_488B <- match(all_variant_keys, vcf_488B$variant_key)
idx_489 <- match(all_variant_keys, vcf_489$variant_key)

message(sprintf("[%s] ├─ Union: 488B=%d variants, 489=%d variants → %d total",
  format(Sys.time(), "%H:%M:%S"),
  nrow(vcf_488B), nrow(vcf_489), n_variants_union))

# Get matrices
alt_488B_mat <- tissue_data[["488B"]]$alt_mat
ref_488B_mat <- tissue_data[["488B"]]$ref_mat
alt_489_mat <- tissue_data[["489"]]$alt_mat
ref_489_mat <- tissue_data[["489"]]$ref_mat

n_cells_488B <- ncol(alt_488B_mat)
n_cells_489 <- ncol(alt_489_mat)
n_cells_total <- n_cells_488B + n_cells_489

message(sprintf("[%s] ├─ Building sparse matrices (triplet format)...", format(Sys.time(), "%H:%M:%S")))

# Extract non-zero entries in triplet format
alt_488B_triplet <- as(alt_488B_mat, "TsparseMatrix")
ref_488B_triplet <- as(ref_488B_mat, "TsparseMatrix")
alt_489_triplet <- as(alt_489_mat, "TsparseMatrix")
ref_489_triplet <- as(ref_489_mat, "TsparseMatrix")

# Filter NAs and build combined triplet coordinates
alt_488B_idx <- idx_488B[alt_488B_triplet@i + 1]
alt_488B_keep <- !is.na(alt_488B_idx)
alt_i_488B <- alt_488B_idx[alt_488B_keep]
alt_j_488B <- alt_488B_triplet@j[alt_488B_keep] + 1
alt_x_488B <- alt_488B_triplet@x[alt_488B_keep]

alt_489_idx <- idx_489[alt_489_triplet@i + 1]
alt_489_keep <- !is.na(alt_489_idx)
alt_i_489 <- alt_489_idx[alt_489_keep]
alt_j_489 <- alt_489_triplet@j[alt_489_keep] + 1 + n_cells_488B
alt_x_489 <- alt_489_triplet@x[alt_489_keep]

alt_i <- c(alt_i_488B, alt_i_489)
alt_j <- c(alt_j_488B, alt_j_489)
alt_x <- c(alt_x_488B, alt_x_489)

# Same for ref
ref_488B_idx <- idx_488B[ref_488B_triplet@i + 1]
ref_488B_keep <- !is.na(ref_488B_idx)
ref_i_488B <- ref_488B_idx[ref_488B_keep]
ref_j_488B <- ref_488B_triplet@j[ref_488B_keep] + 1
ref_x_488B <- ref_488B_triplet@x[ref_488B_keep]

ref_489_idx <- idx_489[ref_489_triplet@i + 1]
ref_489_keep <- !is.na(ref_489_idx)
ref_i_489 <- ref_489_idx[ref_489_keep]
ref_j_489 <- ref_489_triplet@j[ref_489_keep] + 1 + n_cells_488B
ref_x_489 <- ref_489_triplet@x[ref_489_keep]

ref_i <- c(ref_i_488B, ref_i_489)
ref_j <- c(ref_j_488B, ref_j_489)
ref_x <- c(ref_x_488B, ref_x_489)

# Create combined matrices
alt_combined <- sparseMatrix(i = alt_i, j = alt_j, x = alt_x, 
                             dims = c(n_variants_union, n_cells_total), 
                             dimnames = list(NULL, c(colnames(alt_488B_mat), colnames(alt_489_mat))))
ref_combined <- sparseMatrix(i = ref_i, j = ref_j, x = ref_x,
                             dims = c(n_variants_union, n_cells_total),
                             dimnames = list(NULL, c(colnames(ref_488B_mat), colnames(ref_489_mat))))

message(sprintf(
  "[%s] └─ Matrices combined: %d variants × %d cells (488B:%d + 489:%d)\n",
  format(Sys.time(), "%H:%M:%S"),
  nrow(alt_combined), ncol(alt_combined),
  n_cells_488B, n_cells_489
))

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: Combine VCF (vectorized approach)
# ─────────────────────────────────────────────────────────────────────────────
cat("█ STEP 3: Combining VCF files\n")
message(sprintf("[%s] ├─ Vectorized VCF indexing for %d variants...", format(Sys.time(), "%H:%M:%S"), n_variants_union))

vcf_cols <- colnames(vcf_488B)[!colnames(vcf_488B) %in% "variant_key"]
vcf_combined <- data.frame(matrix(NA_character_, nrow = n_variants_union, ncol = length(vcf_cols)),
                          stringsAsFactors = FALSE)
colnames(vcf_combined) <- vcf_cols

# Vectorized indexing (much faster than loops)
use_488B <- !is.na(idx_488B)
use_489 <- !use_488B

safe_idx_488B <- ifelse(is.na(idx_488B), 1, idx_488B)
safe_idx_489 <- ifelse(is.na(idx_489), 1, idx_489)

vcf_combined[use_488B, ] <- vcf_488B[safe_idx_488B[use_488B], vcf_cols]
vcf_combined[use_489, ] <- vcf_489[safe_idx_489[use_489], vcf_cols]

message(sprintf(
  "[%s] └─ VCF combined: %d variants (488B=%d, 489=%d)\n",
  format(Sys.time(), "%H:%M:%S"),
  nrow(vcf_combined), sum(use_488B), sum(use_489)
))

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4: Combine barcodes
# ─────────────────────────────────────────────────────────────────────────────
cat("█ STEP 4: Combining barcodes\n")

barcodes_488B <- tissue_data[["488B"]]$barcodes
barcodes_489 <- tissue_data[["489"]]$barcodes
barcodes_combined <- c(barcodes_488B, barcodes_489)

colnames(alt_combined) <- barcodes_combined
colnames(ref_combined) <- barcodes_combined

message(sprintf(
  "[%s] └─ Barcodes combined: %d total (488B=%d + 489=%d)\n",
  format(Sys.time(), "%H:%M:%S"),
  length(barcodes_combined), length(barcodes_488B), length(barcodes_489)
))

# ─────────────────────────────────────────────────────────────────────────────
# STEP 5: Combine fragment counts
# ─────────────────────────────────────────────────────────────────────────────
cat("█ STEP 5: Combining fragment counts\n")

raw_counts_488B <- tissue_data[["488B"]]$raw_counts
raw_counts_489 <- tissue_data[["489"]]$raw_counts

if (!identical(rownames(raw_counts_488B), rownames(raw_counts_489))) {
  stop("ERROR: Fragment bin names differ between tissues", call. = FALSE)
}

raw_counts_combined <- cbind(raw_counts_488B, raw_counts_489)
colnames(raw_counts_combined) <- barcodes_combined

message(sprintf(
  "[%s] └─ Fragments combined: %d bins × %d cells (488B=%d + 489=%d)\n",
  format(Sys.time(), "%H:%M:%S"),
  nrow(raw_counts_combined), ncol(raw_counts_combined),
  ncol(raw_counts_488B), ncol(raw_counts_489)
))

# ─────────────────────────────────────────────────────────────────────────────
# STEP 6: Use segmentation table from 488B (no regeneration needed)
# ─────────────────────────────────────────────────────────────────────────────
cat("█ STEP 6: Segmentation table\n")

seg_combined <- tissue_data[["488B"]]$seg_table

message(sprintf(
  "[%s] └─ Using 488B segmentation: %d segments\n",
  format(Sys.time(), "%H:%M:%S"),
  nrow(seg_combined)
))

# ─────────────────────────────────────────────────────────────────────────────
# STEP 7: Write outputs
# ─────────────────────────────────────────────────────────────────────────────
cat("█ STEP 7: Writing output files\n")

message(sprintf("[%s] ├─ Writing alt_all.mtx...", format(Sys.time(), "%H:%M:%S")))
writeMM(alt_combined, file.path(output_dir, "alt_all.mtx"))

message(sprintf("[%s] ├─ Writing ref_all.mtx...", format(Sys.time(), "%H:%M:%S")))
writeMM(ref_combined, file.path(output_dir, "ref_all.mtx"))

message(sprintf("[%s] ├─ Writing barcodes.tsv (%d barcodes)...", format(Sys.time(), "%H:%M:%S"), length(barcodes_combined)))
writeLines(barcodes_combined, file.path(output_dir, "barcodes.tsv"))

message(sprintf("[%s] ├─ Writing var_all.vcf (%d variants)...", format(Sys.time(), "%H:%M:%S"), nrow(vcf_combined)))
write.table(
  vcf_combined,
  file.path(output_dir, "var_all.vcf"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

message(sprintf("[%s] ├─ Writing chr1000k_fragments.tsv...", format(Sys.time(), "%H:%M:%S")))
df <- as.data.frame(as.matrix(raw_counts_combined))
df <- cbind(row.names(df), df)
colnames(df)[1] <- "bin"
write.table(df, file.path(output_dir, "chr1000k_fragments.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message(sprintf("[%s] ├─ Writing seg_table.rds...", format(Sys.time(), "%H:%M:%S")))
saveRDS(seg_combined, file.path(output_dir, "seg_table.rds"))

message(sprintf("[%s] └─ Output files written\n", format(Sys.time(), "%H:%M:%S")))

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
cat("╔════════════════════════════════════════════════════════════════╗\n")
message(sprintf("[%s] ✓ COMPLETE: All inputs combined successfully", format(Sys.time(), "%H:%M:%S")))
cat("╠════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("│ Output directory: %s\n", output_dir))
cat(sprintf("│ Total variants:   %d (488B:%d + 489:%d)\n", n_variants_union, nrow(vcf_488B), nrow(vcf_489)))
cat(sprintf("│ Total cells:      %d (488B:%d + 489:%d)\n", n_cells_total, n_cells_488B, n_cells_489))
cat(sprintf("│ Fragment bins:    %d\n", nrow(raw_counts_combined)))
cat(sprintf("│ Segments:         %d\n", nrow(seg_combined)))
cat("╚════════════════════════════════════════════════════════════════╝\n")
message(sprintf("[%s] Ready for Alleloscope analysis!", format(Sys.time(), "%H:%M:%S")))

invisible(
  list(
    alt_all = file.path(output_dir, "alt_all.mtx"),
    ref_all = file.path(output_dir, "ref_all.mtx"),
    barcodes = file.path(output_dir, "barcodes.tsv"),
    var_all = file.path(output_dir, "var_all.vcf"),
    raw_counts = file.path(output_dir, "chr1000k_fragments.tsv"),
    seg_table = file.path(output_dir, "seg_table.rds")
  )
)

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Load tissue-specific data
# ─────────────────────────────────────────────────────────────────────────────
message(sprintf("[%s] Step 1: Loading tissue-specific data", Sys.time()))

tissues <- c("488B", "489")
tissue_data <- list()

for (tissue in tissues) {
  tissue_dir <- file.path(project_root, "Data/04_analysis/cnv/alleloscope/lowseq_tissue_from_existing", tissue)
  
  message(sprintf("[%s] Loading %s from %s", Sys.time(), tissue, tissue_dir))
  
  # Load matrices
  alt_mat <- readMM(file.path(tissue_dir, "alt_all.mtx"))
  ref_mat <- readMM(file.path(tissue_dir, "ref_all.mtx"))
  
  # Load barcodes
  barcodes <- readLines(file.path(tissue_dir, "barcodes.tsv"))
  barcodes <- barcodes[nzchar(barcodes)]
  colnames(alt_mat) <- barcodes
  colnames(ref_mat) <- barcodes
  
  # Load VCF
  vcf_df <- read.table(
    file.path(tissue_dir, "var_all.vcf"),
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "#",
    skip = 1  # Skip the column header line (#CHROM  POS  ...)
  )
  
  # Load raw counts
  raw_counts <- read.table(
    file.path(tissue_dir, "chr1000k_fragments.tsv"),
    sep = "\t",
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
  )
  colnames(raw_counts) <- gsub("[.]", "-", colnames(raw_counts))
  
  # Load segmentation table
  seg_table <- readRDS(file.path(tissue_dir, "seg_table.rds"))
  
  tissue_data[[tissue]] <- list(
    alt_mat = alt_mat,
    ref_mat = ref_mat,
    barcodes = barcodes,
    vcf_df = vcf_df,
    raw_counts = raw_counts,
    seg_table = seg_table
  )
  
  message(sprintf(
    "[%s]   %s: alt=%d variants x %d cells, ref=%d variants x %d cells",
    Sys.time(), tissue,
    nrow(alt_mat), ncol(alt_mat),
    nrow(ref_mat), ncol(ref_mat)
  ))
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Combine SNP matrices efficiently (union variants, cbind cells)
# ─────────────────────────────────────────────────────────────────────────────
message(sprintf("[%s] Step 2: Combining SNP matrices (union of variants, cbind cells)", Sys.time()))

vcf_488B <- tissue_data[["488B"]]$vcf_df
vcf_489 <- tissue_data[["489"]]$vcf_df

message(sprintf(
  "[%s] Tissue VCF sizes: 488B=%d variants, 489=%d variants",
  Sys.time(), nrow(vcf_488B), nrow(vcf_489)
))

# Create variant key (chr:pos:ref:alt) for matching
vcf_488B$variant_key <- with(vcf_488B, paste0(V1, ":", V2, ":", V4, ":", V5))
vcf_489$variant_key <- with(vcf_489, paste0(V1, ":", V2, ":", V4, ":", V5))

# Get union of all variants
all_variant_keys <- unique(c(vcf_488B$variant_key, vcf_489$variant_key))
n_variants_union <- length(all_variant_keys)

message(sprintf("[%s] Union of variants: %d unique variants", Sys.time(), n_variants_union))

# Create index mapping from union variants to tissue-specific variants
idx_488B <- match(all_variant_keys, vcf_488B$variant_key)
idx_489 <- match(all_variant_keys, vcf_489$variant_key)

# Get matrices and cell counts
alt_488B_mat <- tissue_data[["488B"]]$alt_mat
ref_488B_mat <- tissue_data[["488B"]]$ref_mat
alt_489_mat <- tissue_data[["489"]]$alt_mat
ref_489_mat <- tissue_data[["489"]]$ref_mat

n_cells_488B <- ncol(alt_488B_mat)
n_cells_489 <- ncol(alt_489_mat)
n_cells_total <- n_cells_488B + n_cells_489

# Build combined matrices using triplet format for efficiency
message(sprintf("[%s] Building combined matrices from triplet format...", Sys.time()))

# Extract non-zero entries from tissue matrices in triplet format
alt_488B_triplet <- as(alt_488B_mat, "TsparseMatrix")
ref_488B_triplet <- as(ref_488B_mat, "TsparseMatrix")
alt_489_triplet <- as(alt_489_mat, "TsparseMatrix")
ref_489_triplet <- as(ref_489_mat, "TsparseMatrix")

# Collect triplet entries for both tissues
# Filter out NAs that occur when variants don't exist in a tissue
alt_488B_idx <- idx_488B[alt_488B_triplet@i + 1]
alt_488B_keep <- !is.na(alt_488B_idx)
alt_i_488B <- alt_488B_idx[alt_488B_keep]
alt_j_488B <- alt_488B_triplet@j[alt_488B_keep] + 1
alt_x_488B <- alt_488B_triplet@x[alt_488B_keep]

alt_489_idx <- idx_489[alt_489_triplet@i + 1]
alt_489_keep <- !is.na(alt_489_idx)
alt_i_489 <- alt_489_idx[alt_489_keep]
alt_j_489 <- alt_489_triplet@j[alt_489_keep] + 1 + n_cells_488B
alt_x_489 <- alt_489_triplet@x[alt_489_keep]

alt_i <- c(alt_i_488B, alt_i_489)
alt_j <- c(alt_j_488B, alt_j_489)
alt_x <- c(alt_x_488B, alt_x_489)

# Same for ref
ref_488B_idx <- idx_488B[ref_488B_triplet@i + 1]
ref_488B_keep <- !is.na(ref_488B_idx)
ref_i_488B <- ref_488B_idx[ref_488B_keep]
ref_j_488B <- ref_488B_triplet@j[ref_488B_keep] + 1
ref_x_488B <- ref_488B_triplet@x[ref_488B_keep]

ref_489_idx <- idx_489[ref_489_triplet@i + 1]
ref_489_keep <- !is.na(ref_489_idx)
ref_i_489 <- ref_489_idx[ref_489_keep]
ref_j_489 <- ref_489_triplet@j[ref_489_keep] + 1 + n_cells_488B
ref_x_489 <- ref_489_triplet@x[ref_489_keep]

ref_i <- c(ref_i_488B, ref_i_489)
ref_j <- c(ref_j_488B, ref_j_489)
ref_x <- c(ref_x_488B, ref_x_489)

# Create combined sparse matrices from triplet format
alt_combined <- sparseMatrix(i = alt_i, j = alt_j, x = alt_x, 
                             dims = c(n_variants_union, n_cells_total), 
                             dimnames = list(NULL, c(colnames(alt_488B_mat), colnames(alt_489_mat))))
ref_combined <- sparseMatrix(i = ref_i, j = ref_j, x = ref_x,
                             dims = c(n_variants_union, n_cells_total),
                             dimnames = list(NULL, c(colnames(ref_488B_mat), colnames(ref_489_mat))))

message(sprintf(
  "[%s] Combined SNP matrices: %d variants x %d cells (488B=%d + 489=%d cells)",
  Sys.time(),
  nrow(alt_combined), ncol(alt_combined),
  n_cells_488B, n_cells_489
))

# Build combined VCF with proper variant order
vcf_cols <- colnames(vcf_488B)[!colnames(vcf_488B) %in% "variant_key"]
vcf_combined <- data.frame(matrix(NA_character_, nrow = n_variants_union, ncol = length(vcf_cols)),
                          stringsAsFactors = FALSE)
colnames(vcf_combined) <- vcf_cols

message(sprintf("[%s] Building combined VCF (%d variants) using vectorized indexing...", Sys.time(), n_variants_union))

# Use vectorized indexing to select rows from appropriate tissue VCF
# Create mask for which tissue each variant comes from
use_488B <- !is.na(idx_488B)
use_489 <- !use_488B

# Get safe indices (handle NAs by using 1 as placeholder - will be overwritten)
safe_idx_488B <- ifelse(is.na(idx_488B), 1, idx_488B)
safe_idx_489 <- ifelse(is.na(idx_489), 1, idx_489)

# Combine rows using vectorized indexing
vcf_combined[use_488B, ] <- vcf_488B[safe_idx_488B[use_488B], vcf_cols]
vcf_combined[use_489, ] <- vcf_489[safe_idx_489[use_489], vcf_cols]

message(sprintf("[%s] VCF combined: %d variants (488B=%d, 489=%d)", 
  Sys.time(), nrow(vcf_combined), sum(use_488B), sum(use_489)))

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: Combine barcodes (union, preserving tissue labels)
# ─────────────────────────────────────────────────────────────────────────────
message(sprintf("[%s] Step 3: Combining barcodes", Sys.time()))

barcodes_488B <- tissue_data[["488B"]]$barcodes
barcodes_489 <- tissue_data[["489"]]$barcodes

# Check for overlap (shouldn't happen, but good to know)
overlap <- intersect(barcodes_488B, barcodes_489)
if (length(overlap) > 0) {
  message(sprintf("[%s] WARNING: Found %d overlapping barcodes", Sys.time(), length(overlap)))
}

barcodes_combined <- c(barcodes_488B, barcodes_489)
message(sprintf(
  "[%s] Combined barcodes: %d total (488B=%d, 489=%d, overlap=%d)",
  Sys.time(),
  length(barcodes_combined),
  length(barcodes_488B),
  length(barcodes_489),
  length(overlap)
))

# Update matrix column names for combined matrix
colnames(alt_combined) <- barcodes_combined
colnames(ref_combined) <- barcodes_combined

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4: Combine raw fragment counts (cbind, matching barcode order)
# ─────────────────────────────────────────────────────────────────────────────
message(sprintf("[%s] Step 4: Combining fragment counts", Sys.time()))

raw_counts_488B <- tissue_data[["488B"]]$raw_counts
raw_counts_489 <- tissue_data[["489"]]$raw_counts

# Ensure both have same bins (should be identical for same dataset)
if (!identical(rownames(raw_counts_488B), rownames(raw_counts_489))) {
  stop("ERROR: Fragment bin names differ between tissues", call. = FALSE)
}

# Combine by columns (cells), preserving bin order
raw_counts_combined <- cbind(raw_counts_488B, raw_counts_489)
colnames(raw_counts_combined) <- barcodes_combined

message(sprintf(
  "[%s] Combined fragments: %d bins x %d cells (488B=%d + 489=%d)",
  Sys.time(),
  nrow(raw_counts_combined), ncol(raw_counts_combined),
  ncol(raw_counts_488B),
  ncol(raw_counts_489)
))

# ─────────────────────────────────────────────────────────────────────────────
# STEP 5: Combine segmentation tables
# ─────────────────────────────────────────────────────────────────────────────
message(sprintf("[%s] Step 5: Combining segmentation tables", Sys.time()))

seg_488B <- tissue_data[["488B"]]$seg_table
seg_489 <- tissue_data[["489"]]$seg_table

# Check if seg tables are identical
if (identical(seg_488B, seg_489)) {
  message(sprintf("[%s] Segmentation tables are identical (using shared table)", Sys.time()))
  seg_combined <- seg_488B
} else {
  # Merge seg tables (should be identical if from same reference)
  seg_combined <- unique(rbind(seg_488B, seg_489))
  message(sprintf("[%s] Merged unique segs: %d rows", Sys.time(), nrow(seg_combined)))
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 6: Write combined outputs
# ─────────────────────────────────────────────────────────────────────────────
message(sprintf("[%s] Step 6: Writing combined outputs", Sys.time()))

# Write matrices
writeMM(alt_combined, file.path(output_dir, "alt_all.mtx"))
writeMM(ref_combined, file.path(output_dir, "ref_all.mtx"))
message(sprintf("[%s] Wrote alt_all.mtx and ref_all.mtx", Sys.time()))

# Write barcodes
writeLines(barcodes_combined, file.path(output_dir, "barcodes.tsv"))
message(sprintf("[%s] Wrote barcodes.tsv (%d barcodes)", Sys.time(), length(barcodes_combined)))

# Write VCF (using the (identical) combined VCF)
write.table(
  vcf_combined,
  file.path(output_dir, "var_all.vcf"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
message(sprintf("[%s] Wrote var_all.vcf (%d variants)", Sys.time(), nrow(vcf_combined)))

# Write fragment counts
write_counts_tsv <- function(counts_matrix, output_file) {
  df <- as.data.frame(as.matrix(counts_matrix))
  df <- cbind(row.names(df), df)
  colnames(df)[1] <- "bin"
  write.table(df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

write_counts_tsv(raw_counts_combined, file.path(output_dir, "chr1000k_fragments.tsv"))
message(sprintf("[%s] Wrote chr1000k_fragments.tsv", Sys.time()))

# Write segmentation table
saveRDS(seg_combined, file.path(output_dir, "seg_table.rds"))
message(sprintf("[%s] Wrote seg_table.rds", Sys.time()))

message(sprintf("[%s] === Alleloscope combined prep COMPLETE ===", Sys.time()))
message(sprintf("[%s] All inputs ready in: %s", Sys.time(), output_dir))

invisible(
  list(
    alt_all = file.path(output_dir, "alt_all.mtx"),
    ref_all = file.path(output_dir, "ref_all.mtx"),
    barcodes = file.path(output_dir, "barcodes.tsv"),
    var_all = file.path(output_dir, "var_all.vcf"),
    raw_counts = file.path(output_dir, "chr1000k_fragments.tsv"),
    seg_table = file.path(output_dir, "seg_table.rds")
  )
)
