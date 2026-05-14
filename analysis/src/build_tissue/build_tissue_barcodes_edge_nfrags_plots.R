suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

log_msg <- function(tag, msg) {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), tag, msg))
}

project_root <- Sys.getenv("PROJECT_ROOT", "/projectnb/paxlab/presh/projects/spatial_atac")
spatial_candidates <- c(
  file.path(project_root, "Data", "tissue_positions_list.csv"),
  file.path(project_root, "Data", "tissue_positions_list.csv.lnk"),
  file.path(project_root, "Data", "01_inputs", "spatial", "tissue_positions_list.csv"),
  file.path(project_root, "Data", "01_inputs", "spatial", "tissue_positions_list.csv.lnk")
)
spatial_file <- spatial_candidates[file.exists(spatial_candidates)][1]
if (is.na(spatial_file) || !nzchar(spatial_file)) {
  spatial_file <- spatial_candidates[1]
}
out_barcode_dir <- Sys.getenv("OUT_BARCODE_DIR", file.path(project_root, "Data", "01_inputs", "barcodes", "tissue_barcodes"))
out_plot_dir <- Sys.getenv("OUT_PLOT_DIR", file.path(project_root, "analysis", "plots", "variant_qc", "edge_effect_nfrags"))

fragment_files <- list(
  deepseq_488B = NULL,
  lowseq_488B = NULL,
  deepseq_489 = NULL,
  lowseq_489 = NULL
)

# Helper: locate fragment file by trying preferred locations and patterns
find_fragment <- function(name_variants) {
  search_dirs <- c(
    file.path(project_root, "Data", "01_inputs", "fragments"),
    file.path(project_root, "Data")
  )
  for (d in search_dirs) {
    for (pat in name_variants) {
      p <- file.path(d, pat)
      # exact match
      if (file.exists(p)) return(p)
      # try glob
      g <- Sys.glob(file.path(d, pat))
      if (length(g) > 0) return(g[1])
    }
  }
  # last resort: try any file in project Data matching the pattern anywhere
  for (pat in name_variants) {
    g <- Sys.glob(file.path(project_root, "**", pat))
    if (length(g) > 0) return(g[1])
  }
  NULL
}

# Preferred name patterns for each object
fragment_files$deepseq_488B <- find_fragment(c(
  "deepseq_488B/deepseq_488B.fragments.sort.filtered.bed.gz",
  "deepseq_488B/deepseq_488B.fragments.sort.filtered.bed",
  "deepseq_488B.fragments.sort.filtered.bed.gz.lnk",
  "deepseq_488B.fragments.sort.filtered.bed.lnk",
  "deepseq_488B.fragments.sort.filtered.bed.gz",
  "deepseq_488B.fragments.sort.filtered.bed",
  "deepseq_488B.fragments.from_bam.*.bed.gz",
  "deepseq_488B/*.fragments*.bed*",
  "deepseq.fragments.sort.filtered.bed.gz",
  "deepseq.fragments.sort.filtered.bed.gz.lnk",
  "deepseq.fragments.sort.filtered.bed.gzip.gz",
  "deepseq.fragments.sort.filtered.bed"
))
fragment_files$lowseq_488B <- find_fragment(c(
  "lowseq_488B/lowseq_488B.fragments.sort.filtered.bed.gz",
  "lowseq_488B/lowseq_488B.fragments.sort.filtered.bed",
  "lowseq_488B.fragments.sort.filtered.bed.gz.lnk",
  "lowseq_488B.fragments.sort.filtered.bed.lnk",
  "lowseq_488B.fragments.sort.filtered.bed.gz",
  "lowseq_488B.fragments.sort.filtered.bed",
  "lowseq_488B.fragments.from_bam.*.bed.gz",
  "lowseq_488B/*.fragments*.bed*",
  "lowseq.fragments.sort.filtered.bed",
  "lowseq.fragments.sort.filtered.bed.lnk",
  "lowseq.fragments.sort.filtered.bed.gz",
  "lowseq.fragments.sort.filtered.bed.gz.lnk",
  "lowseq.fragments.sort.filtered.bed.gzip.gz"
))
fragment_files$deepseq_489 <- find_fragment(c(
  "deepseq_489/deepseq_489.fragments.sort.filtered.bed.gz",
  "deepseq_489/deepseq_489.fragments.sort.filtered.bed",
  "deepseq_489.fragments.sort.filtered.bed.gz.lnk",
  "deepseq_489.fragments.sort.filtered.bed.lnk",
  "deepseq_489.fragments.sort.filtered.bed.gz",
  "deepseq_489.fragments.sort.filtered.bed",
  "deepseq_489.fragments.from_bam.*.bed.gz",
  "deepseq_489/*.fragments*.bed*",
  "deepseq_489.fragments*.bed*",
  "deepseq_489*fragments*"
))
fragment_files$lowseq_489 <- find_fragment(c(
  "lowseq_489/lowseq_489.fragments.sort.filtered.bed.gz",
  "lowseq_489/lowseq_489.fragments.sort.filtered.bed",
  "lowseq_489.fragments.sort.filtered.bed.gz.lnk",
  "lowseq_489.fragments.sort.filtered.bed.lnk",
  "lowseq_489.fragments.sort.filtered.bed.gz",
  "lowseq_489.fragments.sort.filtered.bed",
  "lowseq_489.fragments.from_bam.*.bed.gz",
  "lowseq_489/*.fragments*.bed*",
  "lowseq_489.fragments*.bed*",
  "lowseq_489*fragments*"
))

# Parameters (override via environment variables when needed)
y_tissue_cutoff <- as.numeric(Sys.getenv("Y_TISSUE_CUTOFF", "4000"))
edge_row_fraction <- as.numeric(Sys.getenv("EDGE_ROW_FRACTION", "0.015"))
min_edge_rows <- as.integer(Sys.getenv("MIN_EDGE_ROWS", "2"))
max_edge_rows <- as.integer(Sys.getenv("MAX_EDGE_ROWS", "5"))
upper_quantile <- as.numeric(Sys.getenv("UPPER_QUANTILE", "0.995"))
mad_multiplier <- as.numeric(Sys.getenv("MAD_MULTIPLIER", "8"))
force_recount <- as.integer(Sys.getenv("FORCE_RECOUNT", "0"))
threshold_round_to <- as.numeric(Sys.getenv("THRESHOLD_ROUND_TO", "5000"))
threshold_rule <- Sys.getenv("THRESHOLD_RULE", "quantile")
# Choose axis for edge detection: "row" (top/bottom) or "col" (left/right)
edge_axis <- Sys.getenv("EDGE_AXIS", "row")

log_msg("start", "build_tissue_barcodes_edge_nfrags_plots.R")
log_msg("start", sprintf("PROJECT_ROOT=%s", project_root))
log_msg("start", sprintf("Y_TISSUE_CUTOFF=%s", y_tissue_cutoff))
log_msg("start", sprintf("EDGE_ROW_FRACTION=%s", edge_row_fraction))
log_msg("start", sprintf("MIN_EDGE_ROWS=%s MAX_EDGE_ROWS=%s", min_edge_rows, max_edge_rows))
log_msg("start", sprintf("UPPER_QUANTILE=%s MAD_MULTIPLIER=%s", upper_quantile, mad_multiplier))
log_msg("start", sprintf("FORCE_RECOUNT=%s", force_recount))
log_msg("start", sprintf("THRESHOLD_RULE=%s THRESHOLD_ROUND_TO=%s", threshold_rule, threshold_round_to))

dir.create(out_barcode_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_plot_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(spatial_file)) {
  stop(sprintf("Missing spatial file: %s", spatial_file), call. = FALSE)
}

for (obj in names(fragment_files)) {
  frag <- fragment_files[[obj]]
  if (is.null(frag) || length(frag) == 0L || !file.exists(frag)) {
    stop(sprintf("Missing fragments file for %s: %s", obj, as.character(frag)), call. = FALSE)
  }
}

run_shell <- function(cmd) {
  status <- system2("/bin/bash", c("-o", "pipefail", "-c", cmd))
  if (!identical(status, 0L)) {
    stop(sprintf("Command failed (exit %s): %s", status, cmd), call. = FALSE)
  }
}

spatial_dt <- fread(spatial_file, header = FALSE)
setnames(spatial_dt, c("barcode", "in_tissue", "array_row", "array_col", "x_spatial", "y_spatial"))
spatial_dt <- spatial_dt[in_tissue == 1]
spatial_dt[, barcode := as.character(barcode)]

# User-specified geometry: top tissue is 488B, bottom near origin is 489.
spatial_dt[, tissue := fifelse(y_spatial > y_tissue_cutoff, "488B", "489")]

build_nfrags_from_fragments <- function(object_name, fragments_path) {
  # cache file lives inside a per-object folder under out_barcode_dir
  obj_dir <- file.path(out_barcode_dir, object_name)
  dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)
  cache_file <- file.path(obj_dir, sprintf("%s_nFrags_from_fragments.tsv.gz", object_name))
  if (file.exists(cache_file) && force_recount != 1L) {
    log_msg("resume", sprintf("%s: using cached nFrags table %s", object_name, cache_file))
    dt <- fread(cache_file)
    dt[, barcode := as.character(barcode)]
    dt[, nFrags := as.numeric(nFrags)]
    return(dt)
  }

  log_msg("step", sprintf("%s: counting nFrags per barcode from %s", object_name, fragments_path))
  tmp_file <- tempfile(pattern = sprintf("%s_nfrags_", object_name), fileext = ".tsv")

  # Read fragments directly in R using fread for efficiency
  suppressPackageStartupMessages(library(data.table))
  
  # Use fread with only columns we need (chr, start, end, barcode, count)
  # For gzipped files, fread handles decompression automatically
  dt_frags <- fread(fragments_path, header = FALSE, select = c(4))
  setnames(dt_frags, "barcode")
  dt_frags[, barcode := as.character(barcode)]
  
  # Remove -1 suffix from barcodes
  dt_frags[, barcode := sub("-1$", "", barcode)]
  
  # Count fragments per barcode
  dt_counts <- dt_frags[, .(nFrags = .N), by = barcode]
  
  # Write to temp file
  fwrite(dt_counts, tmp_file, sep = "\t", col.names = FALSE)

  if (!file.exists(tmp_file) || file.info(tmp_file)$size <= 0) {
    stop(sprintf("nFrags counting produced empty output for %s (%s)", object_name, fragments_path), call. = FALSE)
  }

  dt <- dt_counts
  dt <- dt[!is.na(barcode) & nzchar(barcode)]
  dt <- unique(dt, by = "barcode")

  fwrite(dt, cache_file, sep = "\t")
  unlink(tmp_file)
  log_msg("done", sprintf("%s: wrote nFrags cache %s (%d barcodes)", object_name, cache_file, nrow(dt)))

  dt
}

compute_edge_threshold <- function(dt_tissue, axis = "row") {
  if (axis == "col") {
    vals <- dt_tissue$array_col
  } else {
    vals <- dt_tissue$array_row
  }
  rows <- sort(unique(vals))
  if (!length(rows)) {
    return(list(
      edge_n_rows = 0L,
      threshold = NA_real_,
      lower_edge_rows = integer(0),
      upper_edge_rows = integer(0)
    ))
  }

  edge_n_rows <- max(min_edge_rows, round(length(rows) * edge_row_fraction))
  edge_n_rows <- min(max_edge_rows, edge_n_rows)
  edge_n_rows <- max(1L, min(edge_n_rows, length(rows) %/% 2L))

  lower_edge_rows <- rows[seq_len(edge_n_rows)]
  upper_edge_rows <- rows[(length(rows) - edge_n_rows + 1):length(rows)]
  edge_rows <- unique(c(lower_edge_rows, upper_edge_rows))

  # mark edge positions depending on axis
  if (axis == "col") {
    dt_tissue[, is_edge_row := array_col %in% edge_rows]
  } else {
    dt_tissue[, is_edge_row := array_row %in% edge_rows]
  }
  core_vals <- dt_tissue[is_edge_row == FALSE & !is.na(nFrags), nFrags]
  if (length(core_vals) < 50L) {
    core_vals <- dt_tissue[!is.na(nFrags), nFrags]
  }

  q_hi <- as.numeric(quantile(core_vals, probs = upper_quantile, na.rm = TRUE))
  med <- as.numeric(median(core_vals, na.rm = TRUE))
  mad_val <- as.numeric(mad(core_vals, center = med, constant = 1.4826, na.rm = TRUE))
  mad_cut <- med + mad_multiplier * mad_val

  threshold_raw <- switch(
    threshold_rule,
    quantile = q_hi,
    mad = mad_cut,
    max = max(q_hi, mad_cut, na.rm = TRUE),
    min = min(q_hi, mad_cut, na.rm = TRUE),
    q_hi
  )

  if (!is.finite(threshold_raw)) {
    threshold_raw <- as.numeric(quantile(dt_tissue$nFrags, probs = 0.99, na.rm = TRUE))
  }

  threshold <- threshold_raw
  if (is.finite(threshold_round_to) && threshold_round_to > 0) {
    threshold <- ceiling(threshold_raw / threshold_round_to) * threshold_round_to
  }

  list(
    edge_n_rows = edge_n_rows,
    threshold = threshold,
    threshold_raw = threshold_raw,
    q_hi = q_hi,
    mad_cut = mad_cut,
    lower_edge_rows = lower_edge_rows,
    upper_edge_rows = upper_edge_rows
  )
}

plot_before_after <- function(dt_tissue, object_name, dataset, tissue_name, th, out_dir, axis = "row") {
  threshold <- th$threshold
  obj_plot_dir <- file.path(out_dir, object_name)
  dir.create(obj_plot_dir, recursive = TRUE, showWarnings = FALSE)
  before_file <- file.path(obj_plot_dir, sprintf("%s_before_edge_filter.png", object_name))
  after_file <- file.path(obj_plot_dir, sprintf("%s_after_edge_filter.png", object_name))
  hist_file <- file.path(obj_plot_dir, sprintf("%s_nFrags_hist_cutoff.png", object_name))

  # Before plot: show ALL cells, color by nFrags (gray for missing fragments)
  cells_with_frags <- dt_tissue[!is.na(nFrags)]
  cells_without_frags <- dt_tissue[is.na(nFrags)]
  
  p_before <- ggplot(cells_with_frags, aes(x = x_spatial, y = y_spatial, color = nFrags)) +
    geom_point(size = 0.4, alpha = 0.9) +
    # Add cells without fragments in light gray
    {if (nrow(cells_without_frags) > 0) 
      geom_point(data = cells_without_frags, color = "lightgray", size = 0.3, alpha = 0.5)
    } +
    # Emphasize edge-effect cells by plotting them larger and colored by nFrags
    {if (nrow(dt_tissue[is_edge_effect == TRUE & !is.na(nFrags)]) > 0)
      geom_point(data = dt_tissue[is_edge_effect == TRUE & !is.na(nFrags)], aes(x = x_spatial, y = y_spatial, color = nFrags), size = 1.2, stroke = 0)
    } +
    # Overlay hollow black circles to mark edge-effect cells
    geom_point(
      data = dt_tissue[is_edge_effect == TRUE],
      shape = 1,
      color = "black",
      size = 0.9,
      stroke = 0.35
    ) +
    scale_color_viridis_c(option = "C", trans = "log10", na.value = "lightgray") +
    coord_fixed() +
    theme_classic(base_size = 11) +
    labs(
      title = sprintf("%s %s: Before edge filtering (all cells from spatial)", dataset, tissue_name),
      subtitle = sprintf("All %d cells shown; edge %s + nFrags >= %.0f marked", 
                         nrow(dt_tissue), ifelse(axis == "col", "cols", "rows"), threshold),
      x = "x_spatial",
      y = "y_spatial",
      color = "nFrags"
    )

  kept <- dt_tissue[is_edge_effect == FALSE]
  p_after <- ggplot(kept, aes(x = x_spatial, y = y_spatial, color = nFrags)) +
    geom_point(size = 0.4, alpha = 0.9) +
    scale_color_viridis_c(option = "C", trans = "log10") +
    coord_fixed() +
    theme_classic(base_size = 11) +
    labs(
      title = sprintf("%s %s: After edge filtering", dataset, tissue_name),
      subtitle = sprintf("Kept %d / %d cells", nrow(kept), nrow(dt_tissue)),
      x = "x_spatial",
      y = "y_spatial",
      color = "nFrags"
    )

  # Histogram of nFrags with cutoff line
  dt_hist <- dt_tissue[!is.na(nFrags)]
  p_hist <- ggplot(dt_hist, aes(x = nFrags)) +
    geom_histogram(bins = 100, fill = "grey80", color = "grey50") +
    scale_x_continuous(trans = "log10") +
    geom_vline(xintercept = threshold, color = "red", linetype = "dashed", size = 0.7) +
    theme_classic(base_size = 11) +
    labs(title = sprintf("%s %s: nFrags distribution (log10)", dataset, tissue_name), x = "nFrags", y = "count")

  ggsave(before_file, p_before, width = 8.5, height = 6.5, dpi = 220)
  ggsave(after_file, p_after, width = 8.5, height = 6.5, dpi = 220)
  ggsave(hist_file, p_hist, width = 7, height = 5, dpi = 220)
}

process_combo <- function(object_name, dataset, tissue_name, nfrags_dt) {
  log_msg("step", sprintf("Processing %s %s (object=%s)", dataset, tissue_name, object_name))

  # Include ALL cells from spatial file, even those without fragments
  dt <- merge(
    spatial_dt[tissue == tissue_name, .(barcode, tissue, array_row, array_col, x_spatial, y_spatial)],
    nfrags_dt,
    by = "barcode",
    all.x = TRUE
  )

  # Filter to cells with fragment data for threshold computation only
  dt_with_frags <- dt[!is.na(nFrags)]
  
  if (!nrow(dt_with_frags)) {
    log_msg("warn", sprintf("No cells with fragment data for %s %s", dataset, tissue_name))
    return(data.table(
      dataset = dataset,
      tissue = tissue_name,
      total_cells = nrow(dt),
      edge_cells = 0L,
      kept_cells = nrow(dt),
      edge_n_rows = NA_integer_,
      edge_rows_min = NA_integer_,
      edge_rows_max = NA_integer_,
      nFrags_threshold = NA_real_,
      nFrags_threshold_raw = NA_real_,
      q_hi = NA_real_,
      mad_cut = NA_real_,
      all_file = NA_character_,
      edge_file = NA_character_,
      keep_file = NA_character_
    ))
  }

  th <- compute_edge_threshold(copy(dt_with_frags), edge_axis)
  edge_rows <- unique(c(th$lower_edge_rows, th$upper_edge_rows))
  if (edge_axis == "col") {
    dt[, is_edge_row := array_col %in% edge_rows]
  } else {
    dt[, is_edge_row := array_row %in% edge_rows]
  }
  dt[, is_edge_effect := is_edge_row & is.finite(nFrags) & nFrags >= th$threshold]

  all_barcodes <- sort(unique(dt$barcode))
  edge_barcodes <- sort(unique(dt[is_edge_effect == TRUE, barcode]))
  keep_barcodes <- sort(setdiff(all_barcodes, edge_barcodes))

  obj_dir <- file.path(out_barcode_dir, object_name)
  dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)

  all_file <- file.path(obj_dir, sprintf("%s.barcodes.tsv", object_name))
  edge_file <- file.path(obj_dir, sprintf("%s.edge_effect.barcodes.tsv", object_name))
  keep_file <- file.path(obj_dir, sprintf("%s.no_edge_effect.barcodes.tsv", object_name))

  writeLines(all_barcodes, all_file)
  writeLines(edge_barcodes, edge_file)
  writeLines(keep_barcodes, keep_file)

  plot_before_after(dt, object_name, dataset, tissue_name, th, out_plot_dir, edge_axis)

  log_msg("done", sprintf(
    "%s %s: total=%d edge=%d kept=%d edge_rows=%d threshold_raw=%.1f threshold_rounded=%.0f q_hi=%.1f mad_cut=%.1f",
    dataset, tissue_name, length(all_barcodes), length(edge_barcodes), length(keep_barcodes), th$edge_n_rows,
    th$threshold_raw, th$threshold, th$q_hi, th$mad_cut
  ))

  data.table(
    dataset = dataset,
    tissue = tissue_name,
    total_cells = length(all_barcodes),
    edge_cells = length(edge_barcodes),
    kept_cells = length(keep_barcodes),
    edge_n_rows = th$edge_n_rows,
    edge_rows_min = if (length(edge_rows)) min(edge_rows) else NA_integer_,
    edge_rows_max = if (length(edge_rows)) max(edge_rows) else NA_integer_,
    nFrags_threshold = th$threshold,
    nFrags_threshold_raw = th$threshold_raw,
    q_hi = th$q_hi,
    mad_cut = th$mad_cut,
    all_file = all_file,
    edge_file = edge_file,
    keep_file = keep_file
  )
}

# -----------------------------------------------------------------------------
# Build per-object nFrags tables and process each object
# Objects: deepseq_488B, deepseq_489, lowseq_488B, lowseq_489
# Each object gets its own folder under out_barcode_dir
# -----------------------------------------------------------------------------
objects <- list(
  deepseq_488B = list(dataset = "deepseq", tissue = "488B", fragments = fragment_files$deepseq_488B),
  deepseq_489 = list(dataset = "deepseq", tissue = "489", fragments = fragment_files$deepseq_489),
  lowseq_488B = list(dataset = "lowseq", tissue = "488B", fragments = fragment_files$lowseq_488B),
  lowseq_489 = list(dataset = "lowseq", tissue = "489", fragments = fragment_files$lowseq_489)
)

res_list <- list()
for (obj in names(objects)) {
  info <- objects[[obj]]
  dt_nfrags <- build_nfrags_from_fragments(obj, info$fragments)
  res <- process_combo(obj, info$dataset, info$tissue, dt_nfrags)
  res_list[[obj]] <- res
}

s1 <- res_list$deepseq_488B
s2 <- res_list$deepseq_489
s3 <- res_list$lowseq_488B
s4 <- res_list$lowseq_489

summary_dt <- rbindlist(list(s1, s2, s3, s4), fill = TRUE)
summary_file <- file.path(out_barcode_dir, "edge_effect_nfrags_thresholds.tsv")
fwrite(summary_dt, summary_file, sep = "\t")

log_msg("done", sprintf("Wrote summary table: %s", summary_file))
log_msg("done", sprintf("Barcode files in: %s", out_barcode_dir))
log_msg("done", sprintf("Spatial plots in: %s", out_plot_dir))
log_msg("done", "Completed build_tissue_barcodes_edge_nfrags_plots.R")
