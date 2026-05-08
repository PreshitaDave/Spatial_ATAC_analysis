#!/usr/bin/env Rscript
# Rebuild chr1000k_fragments.tsv for tissue 489 from lowseq fragments BED
# Writes to Data/alleloscope/lowseq_tissue_from_existing/489/chr1000k_fragments.tsv

suppressPackageStartupMessages({
  library(data.table)
})

project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
# Use the tissue-specific 489 fragments (975 MB); the combined lowseq BED only has 488B reads
bed_file <- file.path(project_root, "Data/variant_calling/lowseq/tissue/lowseq_489.fragments.tsv.gz")
out_file <- file.path(project_root, "Data/alleloscope/lowseq_tissue_from_existing/489/chr1000k_fragments.tsv")
bc_file <- file.path(project_root, "Data/tissue_barcodes/lowseq_489.no_edge_effect.barcodes.tsv")

if (!file.exists(bed_file)) stop("Fragments file not found: ", bed_file)
if (!file.exists(bc_file)) stop("Barcode list not found: ", bc_file)

bc <- fread(bc_file, header = FALSE)[[1]]

# Read BED (cols: chr, start, end, barcode, count?) - detect columns
# Read only barcode column first to decide matching strategy (strip '-N' or keep suffix)
bed_bc_only <- fread(bed_file, header = FALSE, sep = "\t", select = 4, data.table = FALSE)[[1]]
bed_bc_nosuf <- sub("-.*$", "", bed_bc_only)
unique_bed_n <- length(unique(bed_bc_only))
matches_nosuf <- sum(unique(bed_bc_nosuf) %in% bc)
matches_withsuf <- sum(unique(bed_bc_only) %in% paste0(bc, "-1"))
cat(sprintf("[rebuild] bed unique barcodes=%d matches_nosuf=%d matches_withsuf=%d\n",
            unique_bed_n, matches_nosuf, matches_withsuf))

# Decide which barcode form to use for matching
if (matches_nosuf >= matches_withsuf && matches_nosuf > 0) {
  use_suffix_form <- FALSE
  bc_match <- bc
  cat("[rebuild] Using stripped bed barcodes for matching\n")
} else if (matches_withsuf > 0) {
  use_suffix_form <- TRUE
  bc_match <- paste0(bc, "-1")
  cat("[rebuild] Using suffixed barcodes (adding '-1') for matching\n")
} else {
  use_suffix_form <- FALSE
  bc_match <- bc
  cat("[rebuild] No barcode overlap detected between BED and tissue list; proceeding but results may be empty\n")
}

# Now read full BED (cols: chr, start, end, barcode, cnt)
bed_dt <- fread(bed_file, header = FALSE, sep = "\t", select = 1:5)
setnames(bed_dt, c("chr","start","end","barcode","cnt"))
if (!use_suffix_form) {
  # Normalize barcode: strip trailing -NN (eg "-1")
  bed_dt[, barcode := sub("-.*$", "", barcode)]
}
# Filter to tissue barcodes
bed_dt <- bed_dt[barcode %in% bc_match]

# Compute bin: 1Mb bins (1-1e6, 1e6+1-2e6, ...)
bin_size <- 1000000L
bed_dt[, bin_idx := floor(start / bin_size)]
bed_dt[, bin_start := bin_idx * bin_size + 1L]
bed_dt[, bin_end := (bin_idx + 1L) * bin_size]
bed_dt[, bin_label := paste0(chr, "-", bin_start, "-", bin_end)]

# Aggregate counts per bin x barcode
agg <- bed_dt[, .(count = .N), by = .(bin_label, barcode)]

# Build full bin list for autosomes 1:22 using existing global file bins if present
global_bins_file <- file.path(project_root, "Data/alleloscope/lowseq/chr1000k_fragments.tsv")
if (file.exists(global_bins_file)) {
  # read first column (row labels) using fread
  g0 <- fread(global_bins_file, nrows = 0)
  global_bins <- fread(global_bins_file, select = 1)[[1]]
  bins <- global_bins
} else {
  # fallback: construct bins for chr1-22 by scanning agg
  bins <- sort(unique(agg$bin_label))
}

# Ensure bins are ordered as in global_bins if available
bins <- unique(bins)

# Efficient dcast: rows=bins, cols=barcodes
# Restrict to known bins and known barcodes first
agg <- agg[bin_label %in% bins & barcode %in% bc_match]
# If suffix-stripped barcodes, the bc column already stripped above
out_wide <- dcast(agg, bin_label ~ barcode, value.var = "count", fill = 0L, fun.aggregate = sum)

# Add any missing bins as all-zero rows
present_bins <- out_wide$bin_label
missing_bins <- setdiff(bins, present_bins)
if (length(missing_bins) > 0) {
  empty_rows <- data.table(bin_label = missing_bins)
  bc_cols <- setdiff(names(out_wide), "bin_label")
  for (col in bc_cols) empty_rows[, (col) := 0L]
  out_wide <- rbindlist(list(out_wide, empty_rows), fill = TRUE)
}
# Reorder rows to match bins order
out_wide <- out_wide[match(bins, out_wide$bin_label), ]

# If barcodes were matched with '-1' suffix, strip it from column names now
if (use_suffix_form) {
  bc_cols_old <- setdiff(names(out_wide), "bin_label")
  bc_cols_new <- sub("-1$", "", bc_cols_old)
  setnames(out_wide, bc_cols_old, bc_cols_new)
}
# Add any missing barcode columns as zeros
bc_out <- bc  # always use no-suffix form as col names
for (col in bc_out) {
  if (!col %in% names(out_wide)) out_wide[, (col) := 0L]
}
# Reorder cols: bin_label first, then barcodes in tissue order
setcolorder(out_wide, c("bin_label", intersect(bc_out, names(out_wide))))
# Replace any NA with 0
for (j in seq_along(out_wide)) {
  set(out_wide, which(is.na(out_wide[[j]])), j, 0L)
}

fwrite(out_wide, out_file, sep = "\t")
cat(sprintf("Wrote %s  rows=%d  cols(barcodes)=%d\n", out_file, nrow(out_wide), ncol(out_wide)-1))
