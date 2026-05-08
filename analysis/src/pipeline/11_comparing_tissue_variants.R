suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
variant_root <- file.path(project_root, "Data/variant_calling")
barcode_root <- file.path(project_root, "Data/alleloscope/barcodes")
out_dir <- Sys.getenv("OUT_DIR", unset = file.path(project_root, "analysis/comparison/tissue_variants"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_info <- function(fmt, ...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%F %T"), sprintf(fmt, ...)))
  flush.console()
}

get_env_int <- function(name, default, min_value = NULL, allow_na = FALSE) {
  raw <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(raw)) {
    value <- default
  } else if (allow_na && toupper(raw) %in% c("NA", "NONE", "NULL")) {
    value <- NA_integer_
  } else {
    value <- suppressWarnings(as.integer(raw))
    if (is.na(value)) {
      stop(sprintf("Environment variable %s must be an integer, got '%s'", name, raw), call. = FALSE)
    }
  }

  if (!is.na(value) && !is.null(min_value) && value < min_value) {
    stop(sprintf("Environment variable %s must be >= %d, got %d", name, min_value, value), call. = FALSE)
  }

  value
}

get_env_csv <- function(name, default) {
  raw <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(raw)) {
    return(default)
  }

  values <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
  values[nzchar(values)]
}

default_workers <- {
  ns <- suppressWarnings(as.integer(Sys.getenv("NSLOTS", unset = "")))
  if (!is.na(ns) && ns > 0L) {
    ns
  } else {
    cores <- parallel::detectCores(logical = FALSE)
    if (is.na(cores) || cores < 1L) 1L else cores
  }
}

datasets <- get_env_csv("DATASETS", c("deepseq", "lowseq"))
chr_start <- get_env_int("CHR_START", 1L, min_value = 1L)
chr_end <- get_env_int("CHR_END", 22L, min_value = 1L)
if (chr_end < chr_start) {
  stop(sprintf("CHR_END (%d) must be >= CHR_START (%d)", chr_end, chr_start), call. = FALSE)
}
chr_numbers <- seq.int(chr_start, chr_end)
progress_every <- get_env_int("PROGRESS_EVERY", 50000L, min_value = 1L)
max_lines_per_chr <- get_env_int("MAX_LINES_PER_CHR", NA_integer_, min_value = 1L, allow_na = TRUE)
n_workers <- get_env_int("N_WORKERS", default_workers, min_value = 1L)

normalize_barcode <- function(x) {
  x <- as.character(x)
  x <- sub("-1$", "", x)
  x
}

make_variant_id <- function(dt) {
  paste(dt$chr, dt$pos, dt$ref, dt$alt, sep = ":")
}

parse_alt_count <- function(token) {
  if (!nzchar(token)) return(0L)
  parts <- strsplit(token, "/", fixed = TRUE)[[1]]
  if (length(parts) < 2) return(0L)
  alt <- suppressWarnings(as.integer(parts[2]))
  if (is.na(alt)) 0L else alt
}

token_has_alt <- function(token) {
  if (!nzchar(token)) return(FALSE)
  slash_pos <- regexpr("/", token, fixed = TRUE)[1]
  if (slash_pos < 1L || slash_pos >= nchar(token)) return(FALSE)
  alt <- suppressWarnings(as.integer(substr(token, slash_pos + 1L, nchar(token))))
  !is.na(alt) && alt > 0L
}

has_alt_in_indices <- function(fields, indices) {
  for (idx in indices) {
    if (idx <= length(fields) && token_has_alt(fields[[idx]])) {
      return(TRUE)
    }
  }
  FALSE
}

read_tissue_barcodes <- function(dataset, tissue_name) {
  f <- file.path(barcode_root, sprintf("%s_%s.barcodes.tsv", dataset, tissue_name))
  if (!file.exists(f)) {
    stop(sprintf("Missing tissue barcode file: %s", f), call. = FALSE)
  }
  bc <- readLines(f, warn = FALSE)
  unique(normalize_barcode(bc[nzchar(bc)]))
}

read_cell_order <- function(dataset, chr = "chr1") {
  f <- file.path(variant_root, dataset, "somatic", sprintf("%s.cell_snv.cellID.filter.csv", chr))
  if (!file.exists(f)) {
    stop(sprintf("Missing cell order file: %s", f), call. = FALSE)
  }
  x <- fread(f)
  normalize_barcode(as.character(x[[ncol(x)]]))
}

collect_tissue_variant_sets <- function(dataset) {
  log_info("Collecting tissue variant sets for %s", dataset)

  bc_488 <- read_tissue_barcodes(dataset, "488B")
  bc_489 <- read_tissue_barcodes(dataset, "489")
  cell_order <- read_cell_order(dataset, "chr1")

  idx_488 <- which(cell_order %in% bc_488)
  idx_489 <- which(cell_order %in% bc_489)

  if (!length(idx_488) || !length(idx_489)) {
    stop(sprintf("No matching barcode indices for %s (488B=%d, 489=%d)", dataset, length(idx_488), length(idx_489)), call. = FALSE)
  }

  process_chromosome <- function(chr_num) {
    chr <- sprintf("chr%d", chr_num)
    snv_file <- file.path(variant_root, dataset, "somatic", sprintf("%s.allSNVs.csv", chr))
    gl_file <- file.path(variant_root, dataset, "somatic", sprintf("%s.gl.filter.hc.cell.mat.gz", chr))

    if (!file.exists(snv_file) || !file.exists(gl_file)) {
      log_info("%s %s skipped (missing file)", dataset, chr)
      return(list(chr = chr, skipped = TRUE, rows = 0L, variants_488 = character(0), variants_489 = character(0)))
    }

    log_info("%s %s started", dataset, chr)

    snv_dt <- fread(snv_file)
    ids <- make_variant_id(snv_dt)
    row_limit <- length(ids)
    if (!is.na(max_lines_per_chr)) {
      row_limit <- min(row_limit, max_lines_per_chr)
    }

    has_alt_488 <- logical(row_limit)
    has_alt_489 <- logical(row_limit)
    n_488 <- 0L
    n_489 <- 0L

    con <- gzfile(gl_file, open = "rt")
    on.exit(close(con), add = TRUE)

    line_i <- 0L
    repeat {
      if (line_i >= row_limit) break
      line <- readLines(con, n = 1L)
      if (!length(line)) break
      line_i <- line_i + 1L
      if (line_i > length(ids)) break

      fields <- strsplit(line, "\t", fixed = TRUE)[[1]]

      if (has_alt_in_indices(fields, idx_488)) {
        has_alt_488[line_i] <- TRUE
        n_488 <- n_488 + 1L
      }

      if (has_alt_in_indices(fields, idx_489)) {
        has_alt_489[line_i] <- TRUE
        n_489 <- n_489 + 1L
      }

      if (line_i %% progress_every == 0L) {
        log_info(
          "%s %s progress: rows=%s/%s, 488B_hits=%s, 489_hits=%s",
          dataset,
          chr,
          format(line_i, big.mark = ","),
          format(row_limit, big.mark = ","),
          format(n_488, big.mark = ","),
          format(n_489, big.mark = ",")
        )
      }
    }

    close(con)
    on.exit(NULL, add = FALSE)

    ids_used <- ids[seq_len(line_i)]
    partial_tag <- if (!is.na(max_lines_per_chr) && row_limit < length(ids)) " [partial]" else ""
    log_info(
      "%s %s finished: rows=%s, 488B_hits=%s, 489_hits=%s%s",
      dataset,
      chr,
      format(line_i, big.mark = ","),
      format(n_488, big.mark = ","),
      format(n_489, big.mark = ","),
      partial_tag
    )

    list(
      chr = chr,
      skipped = FALSE,
      rows = line_i,
      variants_488 = ids_used[has_alt_488[seq_len(line_i)]],
      variants_489 = ids_used[has_alt_489[seq_len(line_i)]]
    )
  }

  run_chromosome_jobs <- function() {
    worker_count <- min(n_workers, length(chr_numbers))
    log_info(
      "%s: scanning %d chromosome(s) with %d worker(s); progress_every=%s%s",
      dataset,
      length(chr_numbers),
      worker_count,
      format(progress_every, big.mark = ","),
      if (!is.na(max_lines_per_chr)) sprintf(", max_lines_per_chr=%s", format(max_lines_per_chr, big.mark = ",")) else ""
    )

    pb <- txtProgressBar(min = 0, max = length(chr_numbers), style = 3)
    on.exit({
      close(pb)
      cat("\n")
    }, add = TRUE)

    if (worker_count <= 1L) {
      results <- vector("list", length(chr_numbers))
      names(results) <- sprintf("chr%d", chr_numbers)
      for (i in seq_along(chr_numbers)) {
        chr_num <- chr_numbers[[i]]
        results[[i]] <- process_chromosome(chr_num)
        setTxtProgressBar(pb, i)
      }
      return(results)
    }

    results <- vector("list", length(chr_numbers))
    names(results) <- sprintf("chr%d", chr_numbers)
    active_jobs <- list()
    active_keys <- list()
    next_idx <- 1L
    completed <- 0L

    launch_job <- function(job_idx) {
      chr_num <- chr_numbers[[job_idx]]
      job <- mcparallel(process_chromosome(chr_num), silent = FALSE)
      key <- as.character(job$pid)
      active_jobs[[key]] <<- job
      active_keys[[key]] <<- sprintf("chr%d", chr_num)
    }

    while (next_idx <= length(chr_numbers) && length(active_jobs) < worker_count) {
      launch_job(next_idx)
      next_idx <- next_idx + 1L
    }

    while (length(active_jobs)) {
      ready <- mccollect(active_jobs, wait = FALSE, timeout = 1)
      if (!length(ready)) {
        next
      }

      for (key in names(ready)) {
        results[[active_keys[[key]]]] <- ready[[key]]
        active_jobs[[key]] <- NULL
        active_keys[[key]] <- NULL
        completed <- completed + 1L
        setTxtProgressBar(pb, completed)

        if (next_idx <= length(chr_numbers)) {
          launch_job(next_idx)
          next_idx <- next_idx + 1L
        }
      }
    }

    results
  }

  chr_results <- run_chromosome_jobs()
  var_488 <- unique(unlist(lapply(chr_results, `[[`, "variants_488"), use.names = FALSE))
  var_489 <- unique(unlist(lapply(chr_results, `[[`, "variants_489"), use.names = FALSE))

  log_info(
    "%s complete: chromosomes=%d, 488B_variants=%s, 489_variants=%s",
    dataset,
    length(chr_results),
    format(length(var_488), big.mark = ","),
    format(length(var_489), big.mark = ",")
  )

  list(
    v488 = var_488,
    v489 = var_489
  )
}

summarize_dataset <- function(dataset, sets) {
  same <- intersect(sets$v488, sets$v489)
  only_488 <- setdiff(sets$v488, sets$v489)
  only_489 <- setdiff(sets$v489, sets$v488)
  any_tissue <- union(sets$v488, sets$v489)

  summary_dt <- data.table(
    dataset = dataset,
    tissue_488B_n = length(sets$v488),
    tissue_489_n = length(sets$v489),
    same_n = length(same),
    only_488B_n = length(only_488),
    only_489_n = length(only_489),
    any_tissue_n = length(any_tissue)
  )

  fwrite(data.table(variant_id = sort(same)), file.path(out_dir, sprintf("%s_same_488B_489.tsv", dataset)), sep = "\t")
  fwrite(data.table(variant_id = sort(only_488)), file.path(out_dir, sprintf("%s_only_488B.tsv", dataset)), sep = "\t")
  fwrite(data.table(variant_id = sort(only_489)), file.path(out_dir, sprintf("%s_only_489.tsv", dataset)), sep = "\t")

  list(summary = summary_dt, same = same, only_488 = only_488, only_489 = only_489, any_tissue = any_tissue)
}

log_info(
  "Starting tissue variant comparison; datasets=%s; chromosomes=%s; out_dir=%s",
  paste(datasets, collapse = ","),
  paste(sprintf("chr%d", chr_numbers), collapse = ","),
  out_dir
)

dataset_sets <- setNames(vector("list", length(datasets)), datasets)
for (dataset in datasets) {
  dataset_sets[[dataset]] <- collect_tissue_variant_sets(dataset)
}

required_datasets <- c("deepseq", "lowseq")
missing_datasets <- setdiff(required_datasets, names(dataset_sets))
if (length(missing_datasets)) {
  stop(sprintf("This script requires results for datasets: %s", paste(required_datasets, collapse = ", ")), call. = FALSE)
}

deep_sets <- dataset_sets[["deepseq"]]
low_sets <- dataset_sets[["lowseq"]]

deep <- summarize_dataset("deepseq", deep_sets)
low <- summarize_dataset("lowseq", low_sets)

overlap_dt <- data.table(
  overlap_type = c(
    "deep_vs_low_any_tissue",
    "deep488_vs_low488",
    "deep489_vs_low489",
    "deep_same_vs_low_same"
  ),
  n_overlap = c(
    length(intersect(deep$any_tissue, low$any_tissue)),
    length(intersect(deep_sets$v488, low_sets$v488)),
    length(intersect(deep_sets$v489, low_sets$v489)),
    length(intersect(deep$same, low$same))
  )
)

summary_all <- rbindlist(list(deep$summary, low$summary))
fwrite(summary_all, file.path(out_dir, "tissue_variant_summary.tsv"), sep = "\t")
fwrite(overlap_dt, file.path(out_dir, "deepseq_lowseq_overlap_counts.tsv"), sep = "\t")

log_info("Wrote outputs to %s", out_dir)
print(summary_all)
print(overlap_dt)
