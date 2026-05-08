suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

log_msg <- function(tag, msg) {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), tag, msg))
}

project_root <- Sys.getenv("PROJECT_ROOT", "/projectnb/paxlab/presh/projects/spatial_atac")
tissue_dir <- file.path(project_root, "Data", "tissue_barcodes")

deep_dir <- file.path(project_root, "Data", "alleloscope", "deepseq")
low_dir <- file.path(project_root, "Data", "alleloscope", "lowseq")
out_root <- file.path(project_root, "Data", "alleloscope", "lowseq_tissue_from_existing")
deep_snv_root <- file.path(project_root, "Data", "alleloscope", "deepseq_tissue_snvs")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(deep_snv_root, recursive = TRUE, showWarnings = FALSE)

read_barcodes <- function(path) {
  if (!file.exists(path)) return(character(0))
  dt <- fread(path, header = FALSE, data.table = FALSE)
  unique(sub("-1$", "", as.character(dt[[1]])))
}

read_var_table <- function(path) {
  dt <- fread(cmd = sprintf("grep -v '^##' %s", shQuote(path)), sep = "\t", header = TRUE)
  setnames(dt, old = "#CHROM", new = "CHROM", skip_absent = TRUE)
  dt[, key := paste(CHROM, POS, REF, ALT, sep = ":")]
  dt
}

write_var_table <- function(var_dt, path, sample_name) {
  hdr <- c(
    "##fileformat=VCFv4.2",
    paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_name, sep = "\t")
  )
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(hdr, con)

  fixed_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "key")
  sample_cols <- setdiff(names(var_dt), fixed_cols)
  sample_col <- if (length(sample_cols)) sample_cols[1] else NULL

  if (is.null(sample_col)) {
    out <- var_dt[, .(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)]
    out[, SAMPLE := "."]
  } else {
    out <- var_dt[, .(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE = get(sample_col))]
  }

  write.table(out, con, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

log_msg("start", "Loading existing deepseq/lowseq Alleloscope matrices")

deep_alt <- readMM(file.path(deep_dir, "alt_all.mtx"))
deep_ref <- readMM(file.path(deep_dir, "ref_all.mtx"))
deep_bc <- read_barcodes(file.path(deep_dir, "barcodes.tsv"))
deep_var <- read_var_table(file.path(deep_dir, "var_all.vcf"))

low_alt <- readMM(file.path(low_dir, "alt_all.mtx"))
low_ref <- readMM(file.path(low_dir, "ref_all.mtx"))
low_bc <- read_barcodes(file.path(low_dir, "barcodes.tsv"))
low_var <- read_var_table(file.path(low_dir, "var_all.vcf"))

if (ncol(deep_alt) != length(deep_bc) || ncol(low_alt) != length(low_bc)) {
  stop("Barcode length does not match matrix columns in existing Alleloscope inputs", call. = FALSE)
}

if (nrow(deep_alt) != nrow(deep_var) || nrow(low_alt) != nrow(low_var)) {
  stop("var_all.vcf row count does not match matrix rows", call. = FALSE)
}

raw_hdr <- fread(file.path(low_dir, "chr1000k_fragments.tsv"), nrows = 0)
raw_cols <- colnames(raw_hdr)

summary_rows <- list()

for (tissue in c("488B", "489")) {
  deep_keep_file <- file.path(tissue_dir, sprintf("deepseq_%s.no_edge_effect.barcodes.tsv", tissue))
  low_keep_file <- file.path(tissue_dir, sprintf("lowseq_%s.no_edge_effect.barcodes.tsv", tissue))

  deep_keep <- read_barcodes(deep_keep_file)
  low_keep <- read_barcodes(low_keep_file)

  if (!length(deep_keep)) {
    deep_keep <- setdiff(
      read_barcodes(file.path(tissue_dir, sprintf("deepseq_%s.barcodes.tsv", tissue))),
      read_barcodes(file.path(tissue_dir, sprintf("deepseq_%s.edge_effect.barcodes.tsv", tissue)))
    )
  }
  if (!length(low_keep)) {
    low_keep <- setdiff(
      read_barcodes(file.path(tissue_dir, sprintf("lowseq_%s.barcodes.tsv", tissue))),
      read_barcodes(file.path(tissue_dir, sprintf("lowseq_%s.edge_effect.barcodes.tsv", tissue)))
    )
  }

  deep_idx <- which(deep_bc %in% deep_keep)
  low_idx <- which(low_bc %in% low_keep)

  if (!length(deep_idx) || !length(low_idx)) {
    log_msg("warn", sprintf("Skipping tissue %s due to empty deep/low barcode overlap", tissue))
    summary_rows[[length(summary_rows) + 1L]] <- data.table(
      tissue = tissue,
      deep_cells = length(deep_idx),
      low_cells = length(low_idx),
      deep_snv_rows = 0L,
      low_common_snv_rows = 0L,
      output_dir = file.path(out_root, tissue)
    )
    next
  }

  # Avoid allocating a large intermediate sparse matrix for deep_total.
  deep_cov <- Matrix::rowSums(deep_alt[, deep_idx, drop = FALSE]) +
    Matrix::rowSums(deep_ref[, deep_idx, drop = FALSE])
  deep_snv_idx <- which(deep_cov > 0)
  deep_var_t <- deep_var[deep_snv_idx]

  deep_snv_file <- file.path(deep_snv_root, sprintf("deepseq_%s.var_all.vcf", tissue))
  write_var_table(deep_var_t, deep_snv_file, sprintf("deepseq_%s", tissue))

  common_keys <- intersect(deep_var_t$key, low_var$key)
  if (!length(common_keys)) {
    log_msg("warn", sprintf("No shared SNP keys for tissue %s between deepseq and lowseq", tissue))
    summary_rows[[length(summary_rows) + 1L]] <- data.table(
      tissue = tissue,
      deep_cells = length(deep_idx),
      low_cells = length(low_idx),
      deep_snv_rows = nrow(deep_var_t),
      low_common_snv_rows = 0L,
      output_dir = file.path(out_root, tissue)
    )
    next
  }

  low_row_idx <- match(common_keys, low_var$key)
  low_row_idx <- low_row_idx[!is.na(low_row_idx)]
  low_var_t <- low_var[low_row_idx]
  low_alt_t <- low_alt[low_row_idx, low_idx, drop = FALSE]
  low_ref_t <- low_ref[low_row_idx, low_idx, drop = FALSE]
  low_bc_t <- low_bc[low_idx]

  tissue_out <- file.path(out_root, tissue)
  dir.create(tissue_out, recursive = TRUE, showWarnings = FALSE)

  writeMM(low_alt_t, file.path(tissue_out, "alt_all.mtx"))
  writeMM(low_ref_t, file.path(tissue_out, "ref_all.mtx"))
  writeLines(low_bc_t, file.path(tissue_out, "barcodes.tsv"))
  write_var_table(low_var_t, file.path(tissue_out, "var_all.vcf"), sprintf("lowseq_%s", tissue))

  seg_src <- file.path(low_dir, "seg_table.rds")
  if (file.exists(seg_src)) {
    file.copy(seg_src, file.path(tissue_out, "seg_table.rds"), overwrite = TRUE)
  }

  row_col_name <- raw_cols[1]  # "cell" in the lowseq fragments file
  keep_raw_cols <- c(row_col_name, low_bc_t)
  keep_raw_cols <- keep_raw_cols[keep_raw_cols %in% raw_cols]
  raw_sub <- fread(file.path(low_dir, "chr1000k_fragments.tsv"), select = keep_raw_cols)
  fwrite(raw_sub, file.path(tissue_out, "chr1000k_fragments.tsv"), sep = "\t")

  log_msg("done", sprintf(
    "tissue=%s low_cells=%d deep_snv=%d low_common_snv=%d out=%s",
    tissue, length(low_bc_t), nrow(deep_var_t), nrow(low_var_t), tissue_out
  ))

  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    tissue = tissue,
    deep_cells = length(deep_idx),
    low_cells = length(low_idx),
    deep_snv_rows = nrow(deep_var_t),
    low_common_snv_rows = nrow(low_var_t),
    output_dir = tissue_out
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# COMBINED 488B + 489: merge both tissue barcode sets into a single Alleloscope
# input so we can run one joint analysis and identify tissue-specific CNVs.
# Uses the union of lowseq 488B + 489 barcodes; SNPs are kept from the lowseq
# global matrix (no deepseq filtering needed for the combined case).
# ─────────────────────────────────────────────────────────────────────────────
log_msg("combined", "Building combined 488B+489 lowseq tissue inputs")

low_488B_file <- file.path(tissue_dir, "lowseq_488B.no_edge_effect.barcodes.tsv")
low_489_file  <- file.path(tissue_dir, "lowseq_489.no_edge_effect.barcodes.tsv")
low_488B_bc <- read_barcodes(low_488B_file)
low_489_bc  <- read_barcodes(low_489_file)

if (!length(low_488B_bc))
  low_488B_bc <- setdiff(read_barcodes(file.path(tissue_dir, "lowseq_488B.barcodes.tsv")),
                         read_barcodes(file.path(tissue_dir, "lowseq_488B.edge_effect.barcodes.tsv")))
if (!length(low_489_bc))
  low_489_bc  <- setdiff(read_barcodes(file.path(tissue_dir, "lowseq_489.barcodes.tsv")),
                         read_barcodes(file.path(tissue_dir, "lowseq_489.edge_effect.barcodes.tsv")))

comb_bc  <- union(low_488B_bc, low_489_bc)
comb_idx <- which(low_bc %in% comb_bc)

if (length(comb_idx) > 0) {
  comb_alt <- low_alt[, comb_idx, drop = FALSE]
  comb_ref <- low_ref[, comb_idx, drop = FALSE]
  comb_bc_t <- low_bc[comb_idx]

  # Keep SNPs covered by at least one combined-tissue cell
  snp_cov <- Matrix::rowSums(comb_alt) + Matrix::rowSums(comb_ref)
  keep_snp <- which(snp_cov > 0)
  comb_alt  <- comb_alt[keep_snp, , drop = FALSE]
  comb_ref  <- comb_ref[keep_snp, , drop = FALSE]
  comb_var  <- low_var[keep_snp]

  comb_out <- file.path(out_root, "combined_488B_489")
  dir.create(comb_out, recursive = TRUE, showWarnings = FALSE)

  writeMM(comb_alt, file.path(comb_out, "alt_all.mtx"))
  writeMM(comb_ref, file.path(comb_out, "ref_all.mtx"))
  writeLines(comb_bc_t, file.path(comb_out, "barcodes.tsv"))
  write_var_table(comb_var, file.path(comb_out, "var_all.vcf"), "lowseq_combined")

  seg_src <- file.path(low_dir, "seg_table.rds")
  if (file.exists(seg_src))
    file.copy(seg_src, file.path(comb_out, "seg_table.rds"), overwrite = TRUE)

  row_col_name <- raw_cols[1]
  keep_raw_comb <- c(row_col_name, comb_bc_t)
  keep_raw_comb <- keep_raw_comb[keep_raw_comb %in% raw_cols]
  raw_comb <- fread(file.path(low_dir, "chr1000k_fragments.tsv"), select = keep_raw_comb)
  fwrite(raw_comb, file.path(comb_out, "chr1000k_fragments.tsv"), sep = "\t")

  # Write a tissue label file so the analysis script can annotate cells
  tissue_labels <- data.table(
    barcode = comb_bc_t,
    tissue  = ifelse(comb_bc_t %in% low_488B_bc, "488B",
                     ifelse(comb_bc_t %in% low_489_bc, "489", "unknown"))
  )
  fwrite(tissue_labels, file.path(comb_out, "tissue_labels.tsv"), sep = "\t")

  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    tissue = "combined_488B_489",
    deep_cells = NA_integer_,
    low_cells = length(comb_bc_t),
    deep_snv_rows = NA_integer_,
    low_common_snv_rows = nrow(comb_var),
    output_dir = comb_out
  )
  log_msg("combined", sprintf(
    "combined cells=%d (488B=%d, 489=%d) snps=%d out=%s",
    length(comb_bc_t), sum(comb_bc_t %in% low_488B_bc),
    sum(comb_bc_t %in% low_489_bc), nrow(comb_var), comb_out
  ))
} else {
  log_msg("warn", "No barcodes found for combined 488B+489 -- skipping")
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
summary_file <- file.path(out_root, "lowseq_tissue_from_existing_summary.tsv")
fwrite(summary_dt, summary_file, sep = "\t")
log_msg("done", sprintf("Wrote summary: %s", summary_file))
