suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

log_msg <- function(tag, msg) {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), tag, msg))
}

project_root <- Sys.getenv("PROJECT_ROOT", "/projectnb/paxlab/presh/projects/spatial_atac")
spatial_file <- file.path(project_root, "Data", "tissue_positions_list.csv")
out_barcode_dir <- file.path(project_root, "Data", "tissue_barcodes")
out_plot_dir <- file.path(project_root, "analysis", "plots", "variant_qc", "edge_effect_nfrags")

fragment_files <- list(
  deepseq = file.path(project_root, "Data", "deepseq.fragments.sort.filtered.bed.gz"),
  lowseq = file.path(project_root, "Data", "lowseq.fragments.sort.filtered.bed")
)

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

for (dataset in names(fragment_files)) {
  if (!file.exists(fragment_files[[dataset]])) {
    stop(sprintf("Missing fragments file for %s: %s", dataset, fragment_files[[dataset]]), call. = FALSE)
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

build_nfrags_from_fragments <- function(dataset, fragments_path) {
  cache_file <- file.path(out_barcode_dir, sprintf("%s_nFrags_from_fragments.tsv.gz", dataset))
  if (file.exists(cache_file) && force_recount != 1L) {
    log_msg("resume", sprintf("%s: using cached nFrags table %s", dataset, cache_file))
    dt <- fread(cache_file)
    dt[, barcode := as.character(barcode)]
    dt[, nFrags := as.numeric(nFrags)]
    return(dt)
  }

  log_msg("step", sprintf("%s: counting nFrags per barcode from %s", dataset, fragments_path))
  tmp_file <- tempfile(pattern = sprintf("%s_nfrags_", dataset), fileext = ".tsv")

  if (grepl("\\.gz$", fragments_path)) {
    if (nzchar(Sys.which("gzip"))) {
      reader <- "gzip -dc"
    } else if (nzchar(Sys.which("zcat"))) {
      reader <- "zcat"
    } else {
      stop("Neither gzip nor zcat is available for gzipped fragments input", call. = FALSE)
    }
  } else {
    reader <- "cat"
  }
  cmd <- sprintf(
    "%s %s | awk 'NF>=4{bc=$4; sub(/-1$/, \"\", bc); n[bc]++} END{for(b in n) print b \"\\t\" n[b]}' > %s",
    reader,
    shQuote(fragments_path),
    shQuote(tmp_file)
  )
  run_shell(cmd)

  if (!file.exists(tmp_file) || file.info(tmp_file)$size <= 0) {
    stop(sprintf("nFrags counting produced empty output for %s (%s)", dataset, fragments_path), call. = FALSE)
  }

  dt <- fread(tmp_file, header = FALSE)
  if (ncol(dt) < 2L) {
    stop(sprintf("Malformed nFrags table for %s: expected 2 columns, found %d", dataset, ncol(dt)), call. = FALSE)
  }
  setnames(dt, c("barcode", "nFrags"))
  dt[, barcode := as.character(barcode)]
  dt[, nFrags := as.numeric(nFrags)]
  dt <- dt[!is.na(barcode) & nzchar(barcode)]
  dt <- unique(dt, by = "barcode")

  fwrite(dt, cache_file, sep = "\t")
  unlink(tmp_file)
  log_msg("done", sprintf("%s: wrote nFrags cache %s (%d barcodes)", dataset, cache_file, nrow(dt)))

  dt
}

compute_edge_threshold <- function(dt_tissue) {
  rows <- sort(unique(dt_tissue$array_row))
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

  dt_tissue[, is_edge_row := array_row %in% edge_rows]
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

plot_before_after <- function(dt_tissue, dataset, tissue_name, threshold, out_dir) {
  before_file <- file.path(out_dir, sprintf("%s_%s_before_edge_filter.png", dataset, tissue_name))
  after_file <- file.path(out_dir, sprintf("%s_%s_after_edge_filter.png", dataset, tissue_name))

  p_before <- ggplot(dt_tissue, aes(x = x_spatial, y = y_spatial)) +
    geom_point(aes(color = nFrags), size = 0.4, alpha = 0.9) +
    geom_point(
      data = dt_tissue[is_edge_effect == TRUE],
      shape = 1,
      color = "black",
      size = 0.9,
      stroke = 0.35
    ) +
    scale_color_viridis_c(option = "C", trans = "log10") +
    coord_fixed() +
    theme_classic(base_size = 11) +
    labs(
      title = sprintf("%s %s: Before edge filtering", dataset, tissue_name),
      subtitle = sprintf("Edge-effect: edge rows + nFrags >= %.0f", threshold),
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

  ggsave(before_file, p_before, width = 8.5, height = 6.5, dpi = 220)
  ggsave(after_file, p_after, width = 8.5, height = 6.5, dpi = 220)
}

process_combo <- function(dataset, tissue_name, nfrags_dt) {
  log_msg("step", sprintf("Processing %s %s", dataset, tissue_name))

  dt <- merge(
    spatial_dt[tissue == tissue_name, .(barcode, tissue, array_row, array_col, x_spatial, y_spatial)],
    nfrags_dt,
    by = "barcode",
    all = FALSE
  )

  if (!nrow(dt)) {
    log_msg("warn", sprintf("No overlapping barcodes for %s %s", dataset, tissue_name))
    return(data.table(
      dataset = dataset,
      tissue = tissue_name,
      total_cells = 0L,
      edge_cells = 0L,
      kept_cells = 0L,
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

  th <- compute_edge_threshold(copy(dt))
  edge_rows <- unique(c(th$lower_edge_rows, th$upper_edge_rows))

  dt[, is_edge_row := array_row %in% edge_rows]
  dt[, is_edge_effect := is_edge_row & is.finite(nFrags) & nFrags >= th$threshold]

  all_barcodes <- sort(unique(dt$barcode))
  edge_barcodes <- sort(unique(dt[is_edge_effect == TRUE, barcode]))
  keep_barcodes <- sort(setdiff(all_barcodes, edge_barcodes))

  all_file <- file.path(out_barcode_dir, sprintf("%s_%s.barcodes.tsv", dataset, tissue_name))
  edge_file <- file.path(out_barcode_dir, sprintf("%s_%s.edge_effect.barcodes.tsv", dataset, tissue_name))
  keep_file <- file.path(out_barcode_dir, sprintf("%s_%s.no_edge_effect.barcodes.tsv", dataset, tissue_name))

  writeLines(all_barcodes, all_file)
  writeLines(edge_barcodes, edge_file)
  writeLines(keep_barcodes, keep_file)

  plot_before_after(dt, dataset, tissue_name, th$threshold, out_plot_dir)

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
# Build dataset-level nFrags tables from fragments
# -----------------------------------------------------------------------------
deepseq_nfrags <- build_nfrags_from_fragments("deepseq", fragment_files$deepseq)
lowseq_nfrags <- build_nfrags_from_fragments("lowseq", fragment_files$lowseq)

# -----------------------------------------------------------------------------
# Section 1: deepseq tissue 488B
# -----------------------------------------------------------------------------
s1 <- process_combo("deepseq", "488B", deepseq_nfrags)

# -----------------------------------------------------------------------------
# Section 2: deepseq tissue 489
# -----------------------------------------------------------------------------
s2 <- process_combo("deepseq", "489", deepseq_nfrags)

# -----------------------------------------------------------------------------
# Section 3: lowseq tissue 488B
# -----------------------------------------------------------------------------
s3 <- process_combo("lowseq", "488B", lowseq_nfrags)

# -----------------------------------------------------------------------------
# Section 4: lowseq tissue 489
# -----------------------------------------------------------------------------
s4 <- process_combo("lowseq", "489", lowseq_nfrags)

summary_dt <- rbindlist(list(s1, s2, s3, s4), fill = TRUE)
summary_file <- file.path(out_barcode_dir, "edge_effect_nfrags_thresholds.tsv")
fwrite(summary_dt, summary_file, sep = "\t")

log_msg("done", sprintf("Wrote summary table: %s", summary_file))
log_msg("done", sprintf("Barcode files in: %s", out_barcode_dir))
log_msg("done", sprintf("Spatial plots in: %s", out_plot_dir))
log_msg("done", "Completed build_tissue_barcodes_edge_nfrags_plots.R")
