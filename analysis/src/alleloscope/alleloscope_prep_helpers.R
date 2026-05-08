suppressPackageStartupMessages({
  library(Matrix)
})

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Install it before running this script.", pkg), call. = FALSE)
  }
}

normalize_barcode <- function(barcodes) {
  barcodes <- gsub("^cell-", "", as.character(barcodes))
  sub("-1$", "", barcodes)
}

hg38_chr_sizes <- function(chrom_sizes_file) {
  require_pkg("data.table")

  if (!file.exists(chrom_sizes_file)) {
    stop(sprintf("Chromosome sizes file not found: %s", chrom_sizes_file), call. = FALSE)
  }

  chr_sizes <- data.table::fread(
    chrom_sizes_file,
    sep = "\t",
    header = FALSE,
    data.table = FALSE,
    col.names = c("chr", "size")
  )
  keep_chr <- paste0("chr", 1:22)
  chr_sizes <- chr_sizes[chr_sizes$chr %in% keep_chr, c("chr", "size"), drop = FALSE]
  chr_sizes <- chr_sizes[match(keep_chr, chr_sizes$chr), , drop = FALSE]

  if (nrow(chr_sizes) != length(keep_chr) || anyNA(chr_sizes$chr)) {
    stop("Chromosome sizes file is missing one or more autosomes chr1-chr22.", call. = FALSE)
  }

  chr_sizes$size <- as.integer(chr_sizes$size)
  chr_sizes
}

natural_chr_order <- function(paths) {
  chr_num <- as.integer(sub(".*chr([0-9]+).*", "\\1", basename(paths)))
  paths[order(chr_num)]
}

first_existing <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (!length(existing)) {
    stop("No matching input files were found.", call. = FALSE)
  }
  existing[1]
}

read_barcode_order <- function(somatic_dir) {
  require_pkg("data.table")

  cell_files <- list.files(
    somatic_dir,
    pattern = "\\.cell_snv\\.cellID\\.filter\\.csv$",
    full.names = TRUE
  )
  if (!length(cell_files)) {
    stop("Could not find any *.cell_snv.cellID.filter.csv files.", call. = FALSE)
  }

  # Read and aggregate barcodes from ALL chromosome files
  cell_files <- natural_chr_order(cell_files)
  all_barcodes <- character(0)
  
  for (f in cell_files) {
    cell_dt <- data.table::fread(f, data.table = FALSE, header = TRUE)
    # Extract first non-index column (usually named 'cell'), not last column
    # as last may be numeric ID. Skip first column if it's an index (e.g. row number).
    bc_col <- if (colnames(cell_dt)[1] %in% c("", "V1", "Index")) {
      2L
    } else {
      1L
    }
    barcodes <- normalize_barcode(cell_dt[[bc_col]])
    # Guard against accidental header-like tokens being treated as barcodes.
    barcodes <- barcodes[!is.na(barcodes) & nzchar(barcodes) & !barcodes %in% c("cell", "barcode", "Cell", "Barcode")]
    all_barcodes <- c(all_barcodes, barcodes)
  }
  
  unique(all_barcodes)
}

vcf_header_lines <- function(sample_name) {
  c(
    "##fileformat=VCFv4.2",
    paste(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_name,
      sep = "\t"
    )
  )
}

build_var_table_from_vcf <- function(vcf_file) {
  require_pkg("data.table")

  vcf_dt <- data.table::fread(
    cmd = sprintf("zcat -f %s | grep -v '^##'", shQuote(vcf_file)),
    sep = "\t",
    header = TRUE,
    data.table = FALSE,
    showProgress = FALSE
  )
  needed <- c("#CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "INFO")
  missing_cols <- setdiff(needed, colnames(vcf_dt))
  if (length(missing_cols)) {
    stop(
      sprintf("Missing required columns in %s: %s", vcf_file, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  sample_col <- colnames(vcf_dt)[ncol(vcf_dt)]
  genotype_col <- if (sample_col %in% c("INFO", "FORMAT")) {
    rep(".", nrow(vcf_dt))
  } else {
    sub(":.*$", "", as.character(vcf_dt[[sample_col]]))
  }

  data.frame(
    CHROM = as.character(vcf_dt[["#CHROM"]]),
    POS = as.integer(vcf_dt$POS),
    ID = ".",
    REF = as.character(vcf_dt$REF),
    ALT = as.character(vcf_dt$ALT),
    QUAL = as.character(vcf_dt$QUAL),
    FILTER = as.character(vcf_dt$FILTER),
    INFO = as.character(vcf_dt$INFO),
    FORMAT = "GT",
    SAMPLE = genotype_col,
    stringsAsFactors = FALSE
  )
}

find_chr_bam <- function(bam_dir, sample_name, chr) {
  first_existing(c(
    file.path(bam_dir, sprintf("%s_%s.filter.bam", sample_name, chr)),
    file.path(bam_dir, sprintf("%s.filter.targeted.bam", chr)),
    file.path(bam_dir, sprintf("%s.filter.bam", chr))
  ))
}

build_sparse_pair_from_vartrix <- function(vcf_file, bam_file, barcodes, barcode_file, fasta_file, vartrix_bin,
                                           work_dir, threads = 1L) {
  if (!file.exists(vartrix_bin)) {
    stop(sprintf("VarTrix binary not found: %s", vartrix_bin), call. = FALSE)
  }

  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  alt_matrix_file <- file.path(work_dir, "alt_matrix.mtx")
  ref_matrix_file <- file.path(work_dir, "ref_matrix.mtx")
  variants_file <- file.path(work_dir, "variants.tsv")
  log_file <- file.path(work_dir, "vartrix.log")

  chr_label <- sub(".*(chr[0-9]+).*", "\\1", basename(vcf_file))

  # Resume: if all three output files exist and are non-empty, skip VarTrix entirely.
  outputs_complete <- all(vapply(
    c(alt_matrix_file, ref_matrix_file, variants_file),
    function(f) file.exists(f) && file.size(f) > 0L,
    logical(1)
  ))
  if (outputs_complete) {
    message(sprintf("[%s] VarTrix skip (outputs already exist): %s", Sys.time(), chr_label))
    alt_mat <- Matrix::readMM(alt_matrix_file)
    ref_mat <- Matrix::readMM(ref_matrix_file)

    # Guard resume integrity: zero-row or mismatched cached outputs are stale and
    # must be recomputed to avoid downstream row/object errors.
    if (nrow(alt_mat) == 0L || nrow(ref_mat) == 0L ||
        nrow(alt_mat) != nrow(ref_mat) || ncol(alt_mat) != ncol(ref_mat)) {
      message(sprintf(
        "[%s] [resume-invalid] Cached VarTrix outputs for %s failed integrity checks; recomputing",
        Sys.time(),
        chr_label
      ))
      outputs_complete <- FALSE
    }

    if (outputs_complete && nrow(alt_mat) > 1000L &&
        Matrix::nnzero(alt_mat) == 0L && Matrix::nnzero(ref_mat) == 0L) {
      message(sprintf(
        "[%s] [resume-invalid] Cached VarTrix outputs for %s have zero non-zero entries; recomputing",
        Sys.time(),
        chr_label
      ))
      outputs_complete <- FALSE
    }

    if (outputs_complete) {
      colnames(alt_mat) <- barcodes
      colnames(ref_mat) <- barcodes
      return(list(
        alt = alt_mat,
        ref = ref_mat,
        n_rows = nrow(alt_mat),
        variants_file = variants_file,
        work_dir = work_dir
      ))
    }
  }

  # Clean up any partial outputs before running so VarTrix does not hit
  # "Output path already exists" on restarts.
  for (f in c(alt_matrix_file, ref_matrix_file, variants_file)) {
    if (file.exists(f)) unlink(f)
  }

  t0 <- Sys.time()
  message(sprintf(
    "[%s] VarTrix start: chr=%s threads=%d bam=%s",
    Sys.time(),
    chr_label,
    as.integer(threads),
    basename(bam_file)
  ))

  attempt_threads <- unique(as.integer(c(threads, max(1L, as.integer(threads) %/% 2L))))
  status <- 1L
  attempt_log <- log_file

  for (attempt_i in seq_along(attempt_threads)) {
    th <- attempt_threads[[attempt_i]]
    attempt_log <- if (attempt_i == 1L) {
      log_file
    } else {
      file.path(work_dir, sprintf("vartrix.retry%d.log", attempt_i))
    }

    args <- c(
      "-v", vcf_file,
      "-b", bam_file,
      "-f", fasta_file,
      "-c", barcode_file,
      "-o", alt_matrix_file,
      "--ref-matrix", ref_matrix_file,
      "--out-variants", variants_file,
      "-s", "coverage",
      "--threads", as.character(th),
      "--log-level", "info"
    )

    if (attempt_i > 1L) {
      message(sprintf(
        "[%s] VarTrix retry %d/%d for %s with threads=%d",
        Sys.time(),
        attempt_i,
        length(attempt_threads),
        basename(vcf_file),
        th
      ))
    }

    status <- system2(vartrix_bin, args = args, stdout = attempt_log, stderr = attempt_log)
    if (identical(status, 0L)) {
      break
    }
  }

  if (!identical(status, 0L)) {
    stop(
      sprintf(
        "VarTrix failed for %s after %d attempt(s). Last log: %s",
        basename(vcf_file),
        length(attempt_threads),
        attempt_log
      ),
      call. = FALSE
    )
  }

  alt_mat <- Matrix::readMM(alt_matrix_file)
  ref_mat <- Matrix::readMM(ref_matrix_file)

  if (nrow(alt_mat) > 1000L && Matrix::nnzero(alt_mat) == 0L && Matrix::nnzero(ref_mat) == 0L) {
    stop(
      sprintf(
        "VarTrix produced all-zero matrices for %s (rows=%d, cols=%d). Likely barcode/BAM mismatch or missing CB tags.",
        chr_label,
        nrow(alt_mat),
        ncol(alt_mat)
      ),
      call. = FALSE
    )
  }

  colnames(alt_mat) <- barcodes
  colnames(ref_mat) <- barcodes

  elapsed_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message(sprintf(
    "[%s] VarTrix done: chr=%s rows=%d cols=%d elapsed=%.2f min",
    Sys.time(),
    sub(".*(chr[0-9]+).*", "\\1", basename(vcf_file)),
    nrow(alt_mat),
    ncol(alt_mat),
    elapsed_min
  ))

  list(ref = ref_mat, alt = alt_mat, n_rows = nrow(alt_mat))
}

write_vcf <- function(var_table, out_file, sample_name) {
  con <- file(out_file, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(vcf_header_lines(sample_name), con = con)
  utils::write.table(
    var_table,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

make_bins <- function(size_table, bin_size) {
  bins <- lapply(seq_len(nrow(size_table)), function(idx) {
    chr <- size_table$chr[idx]
    chr_size <- as.integer(size_table$size[idx])
    starts <- seq.int(1L, chr_size, by = bin_size)
    ends <- pmin(starts + bin_size - 1L, chr_size)
    data.frame(
      chr = chr,
      start = starts,
      end = ends,
      row_name = paste(chr, starts, ends, sep = "-"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, bins)
}

build_raw_counts <- function(fragments_file, barcodes, size_table, bin_size = 1000000L, chunk_n = 200000L) {
  require_pkg("data.table")

  message(sprintf(
    "[%s] Raw counts: start reading %s with chunk_n=%d",
    Sys.time(),
    basename(fragments_file),
    as.integer(chunk_n)
  ))

  bins <- make_bins(size_table, bin_size)
  n_bins <- nrow(bins)
  n_cells <- length(barcodes)

  barcodes <- normalize_barcode(barcodes)
  barcode_index <- seq_along(barcodes)
  names(barcode_index) <- barcodes

  chr_bin_count <- table(bins$chr)
  chr_offsets <- cumsum(c(0L, as.integer(chr_bin_count)))[seq_along(chr_bin_count)]
  names(chr_offsets) <- names(chr_bin_count)

  open_fun <- if (grepl("\\.gz$", fragments_file, ignore.case = TRUE)) gzfile else file
  con <- open_fun(fragments_file, open = "rt")
  on.exit(close(con), add = TRUE)

  all_bin_idx <- integer(0L)
  all_cell_idx <- integer(0L)
  all_counts <- integer(0L)
  chunk_i <- 0L
  total_lines <- 0L
  kept_rows <- 0L
  t0 <- Sys.time()

  repeat {
    lines <- readLines(con, n = chunk_n)
    if (!length(lines)) {
      break
    }
    chunk_i <- chunk_i + 1L
    total_lines <- total_lines + length(lines)

    frag_df <- data.table::fread(
      text = paste(lines, collapse = "\n"),
      sep = "\t",
      header = FALSE,
      select = 1:4,
      showProgress = FALSE,
      data.table = FALSE
    )
    colnames(frag_df) <- c("chr", "start", "end", "barcode")

    frag_df$barcode <- normalize_barcode(frag_df$barcode)
    frag_df <- frag_df[frag_df$chr %in% size_table$chr, , drop = FALSE]
    if (!nrow(frag_df)) {
      if (chunk_i %% 10L == 0L) {
        message(sprintf(
          "[%s] Raw counts progress: chunk=%d lines=%d kept=%d (elapsed=%.2f min)",
          Sys.time(),
          chunk_i,
          total_lines,
          kept_rows,
          as.numeric(difftime(Sys.time(), t0, units = "mins"))
        ))
      }
      next
    }

    frag_df$cell_idx <- unname(barcode_index[frag_df$barcode])
    frag_df <- frag_df[!is.na(frag_df$cell_idx), , drop = FALSE]
    if (!nrow(frag_df)) {
      if (chunk_i %% 10L == 0L) {
        message(sprintf(
          "[%s] Raw counts progress: chunk=%d lines=%d kept=%d (elapsed=%.2f min)",
          Sys.time(),
          chunk_i,
          total_lines,
          kept_rows,
          as.numeric(difftime(Sys.time(), t0, units = "mins"))
        ))
      }
      next
    }
    kept_rows <- kept_rows + nrow(frag_df)

    frag_df$midpoint <- as.integer((frag_df$start + frag_df$end) %/% 2L)
    frag_df$chr_bin <- pmin(frag_df$midpoint %/% bin_size + 1L, as.integer(chr_bin_count[frag_df$chr]))
    frag_df$bin_idx <- unname(chr_offsets[frag_df$chr]) + frag_df$chr_bin

    agg_df <- stats::aggregate(
      x = list(N = rep(1L, nrow(frag_df))),
      by = list(bin_idx = frag_df$bin_idx, cell_idx = frag_df$cell_idx),
      FUN = sum
    )
    all_bin_idx  <- c(all_bin_idx,  agg_df$bin_idx)
    all_cell_idx <- c(all_cell_idx, agg_df$cell_idx)
    all_counts   <- c(all_counts,   agg_df$N)

    if (chunk_i %% 10L == 0L) {
      message(sprintf(
        "[%s] Raw counts progress: chunk=%d lines=%d kept=%d nnz_partial=%d (elapsed=%.2f min)",
        Sys.time(),
        chunk_i,
        total_lines,
        kept_rows,
        length(all_counts),
        as.numeric(difftime(Sys.time(), t0, units = "mins"))
      ))
    }
  }

  counts_sparse <- Matrix::sparseMatrix(
    i = all_bin_idx,
    j = all_cell_idx,
    x = all_counts,
    dims = c(n_bins, n_cells),
    dimnames = list(bins$row_name, barcodes)
  )
  message(sprintf(
    "[%s] Raw counts complete: bins=%d cells=%d nnz=%d total_lines=%d kept=%d elapsed=%.2f min",
    Sys.time(),
    n_bins,
    n_cells,
    length(counts_sparse@x),
    total_lines,
    kept_rows,
    as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ))
  counts_sparse
}

write_counts_tsv <- function(mat, file_path, chunk_rows = 200L) {
  n <- nrow(mat)
  first_chunk <- as.matrix(mat[seq_len(min(chunk_rows, n)), , drop = FALSE])
  utils::write.table(
    first_chunk, file = file_path,
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA
  )
  if (n > chunk_rows) {
    con <- file(file_path, open = "at")
    on.exit(close(con), add = TRUE)
    for (start in seq(chunk_rows + 1L, n, by = chunk_rows)) {
      end <- min(start + chunk_rows - 1L, n)
      chunk <- as.matrix(mat[start:end, , drop = FALSE])
      utils::write.table(
        chunk, file = con,
        sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE
      )
    }
  }
  invisible(NULL)
}

build_chromosome_seg_table <- function(size_table, var_table) {
  snp_counts <- table(var_table$CHROM)
  seg_table <- data.frame(
    chr = size_table$chr,
    start = 0L,
    end = as.integer(size_table$size),
    states = 0,
    length = as.integer(size_table$size),
    mean = 0,
    var = 0,
    Var1 = seq_len(nrow(size_table)),
    Freq = as.integer(snp_counts[size_table$chr]),
    stringsAsFactors = FALSE
  )
  seg_table$Freq[is.na(seg_table$Freq)] <- 0L
  seg_table$chrr <- paste0(seg_table$chr, ":", seg_table$start)
  seg_table
}

prepare_alleloscope_inputs <- function(config) {
  require_pkg("data.table")

  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  message(sprintf("[%s] Starting Alleloscope prep for sample '%s'", Sys.time(), config$sample_name))
  message(sprintf("[%s] Output directory: %s", Sys.time(), config$output_dir))

  required_fields <- c("somatic_dir", "vcf_dir", "bam_dir", "fasta_file", "vartrix_bin", "chrom_sizes_file")
  missing_fields <- required_fields[!vapply(required_fields, function(field) !is.null(config[[field]]), logical(1))]
  if (length(missing_fields)) {
    stop(sprintf("Missing required config fields: %s", paste(missing_fields, collapse = ", ")), call. = FALSE)
  }

  size_table <- hg38_chr_sizes(config$chrom_sizes_file)
  barcodes <- read_barcode_order(config$somatic_dir)

  if (!is.null(config$barcode_subset_file)) {
    if (!file.exists(config$barcode_subset_file)) {
      stop(sprintf("barcode_subset_file not found: %s", config$barcode_subset_file), call. = FALSE)
    }
    subset_barcodes <- normalize_barcode(readLines(config$barcode_subset_file, warn = FALSE))
    subset_barcodes <- unique(subset_barcodes[nzchar(subset_barcodes)])
    old_n <- length(barcodes)
    barcodes <- barcodes[barcodes %in% subset_barcodes]
    message(sprintf(
      "[%s] Applied barcode subset: kept %d/%d barcodes from %s",
      Sys.time(), length(barcodes), old_n, config$barcode_subset_file
    ))
    if (!length(barcodes)) {
      stop("No barcodes remaining after applying barcode_subset_file.", call. = FALSE)
    }
  }

  # Barcode translation for datasets where BAM CB tags differ from somatic cell IDs.
  # deepseq BAMs use 8bp CB tags but somatic cell ID files use 16bp barcodes.
  # Setting vartrix_barcode_length=8L in config truncates barcodes to 8bp for VarTrix;
  # column names are translated back to 16bp after VarTrix so all downstream steps
  # receive 16bp barcodes consistently. barcodes.tsv always contains 16bp for analysis.
  vartrix_barcodes <- barcodes
  if (identical(as.integer(config$vartrix_barcode_length), 8L)) {
    vartrix_barcodes <- substr(barcodes, 1L, 8L)
    message(sprintf(
      "[%s] vartrix_barcode_length=8: truncating to 8bp for VarTrix (%d cells)",
      Sys.time(), length(vartrix_barcodes)
    ))
  }

  barcode_file <- file.path(config$output_dir, "barcodes.tsv")
  writeLines(barcodes, con = barcode_file)
  vartrix_barcode_file <- if (identical(vartrix_barcodes, barcodes)) {
    barcode_file
  } else {
    tmp_bc <- file.path(config$output_dir, "barcodes_vartrix8bp.tsv")
    writeLines(vartrix_barcodes, con = tmp_bc)
    tmp_bc
  }
  message(sprintf("[%s] Wrote barcode list (analysis, 16bp): %s", Sys.time(), barcode_file))
  message(sprintf(
    "[%s] VarTrix barcode file: %s (n=%d, unique=%d)",
    Sys.time(),
    vartrix_barcode_file,
    length(vartrix_barcodes),
    length(unique(vartrix_barcodes))
  ))

  vcf_pattern <- if (is.null(config$vcf_pattern)) {
    "^chr([1-9]|1[0-9]|2[0-2])\\.phased\\.vcf\\.gz$"
  } else {
    as.character(config$vcf_pattern)
  }

  vcf_files <- natural_chr_order(list.files(
    config$vcf_dir,
    pattern = vcf_pattern,
    full.names = TRUE
  ))

  if (!length(vcf_files)) {
    stop(
      sprintf("No VCF files found in %s matching pattern: %s", config$vcf_dir, vcf_pattern),
      call. = FALSE
    )
  }

  message(sprintf(
    "[%s] Found %d VCF files in %s matching pattern: %s",
    Sys.time(),
    length(vcf_files),
    config$vcf_dir,
    vcf_pattern
  ))

  vcf_chr <- sub(".*(chr[0-9]+).*", "\\1", basename(vcf_files))
  common_chr <- unique(vcf_chr)
  common_chr <- common_chr[order(as.integer(sub("chr", "", common_chr)))]

  ref_parts <- list()
  alt_parts <- list()
  var_parts <- list()
  total_chr <- length(common_chr)
  prep_t0 <- Sys.time()

  for (chr_idx in seq_along(common_chr)) {
    chr <- common_chr[[chr_idx]]
    message(sprintf("[%s] [%d/%d] Processing %s ...", Sys.time(), chr_idx, total_chr, chr))
    vcf_file <- vcf_files[match(chr, vcf_chr)]
    bam_file <- find_chr_bam(config$bam_dir, config$sample_name, chr)
    vartrix_dir <- file.path(config$output_dir, "vartrix", chr)

    message(sprintf(
      "[%s] [%d/%d] Running VarTrix on %s with BAM %s (threads=%d)",
      Sys.time(),
      chr_idx,
      total_chr,
      chr,
      basename(bam_file),
      if (is.null(config$vartrix_threads)) 1L else as.integer(config$vartrix_threads)
    ))

    chr_t0 <- Sys.time()

    sparse_pair <- build_sparse_pair_from_vartrix(
      vcf_file = vcf_file,
      bam_file = bam_file,
      barcodes = vartrix_barcodes,
      barcode_file = vartrix_barcode_file,
      fasta_file = config$fasta_file,
      vartrix_bin = config$vartrix_bin,
      work_dir = vartrix_dir,
      threads = if (is.null(config$vartrix_threads)) 1L else config$vartrix_threads
    )
    var_part <- build_var_table_from_vcf(vcf_file)

    if (sparse_pair$n_rows != nrow(var_part)) {
      stop(
        sprintf("Row mismatch for %s: matrix rows %s, SNP rows %s.", chr, sparse_pair$n_rows, nrow(var_part)),
        call. = FALSE
      )
    }

    message(sprintf(
      "[%s] [%d/%d] Completed %s: %d variants, elapsed=%.2f min",
      Sys.time(),
      chr_idx,
      total_chr,
      chr,
      nrow(var_part),
      as.numeric(difftime(Sys.time(), chr_t0, units = "mins"))
    ))
    ref_parts[[chr]] <- sparse_pair$ref
    alt_parts[[chr]] <- sparse_pair$alt
    var_parts[[chr]] <- var_part
  }

  ref_all <- if (length(ref_parts) == 1L) ref_parts[[1]] else do.call(rbind, ref_parts)
  alt_all <- if (length(alt_parts) == 1L) alt_parts[[1]] else do.call(rbind, alt_parts)

  # Translate 8bp VarTrix column names back to 16bp barcodes for downstream Alleloscope steps
  if (!identical(vartrix_barcodes, barcodes)) {
    bc_idx <- match(colnames(ref_all), vartrix_barcodes)
    colnames(ref_all) <- barcodes[bc_idx]
    colnames(alt_all) <- barcodes[bc_idx]
    message(sprintf("[%s] Translated VarTrix 8bp column names back to 16bp", Sys.time()))
  }

  var_all <- do.call(rbind, var_parts)

  message(sprintf("[%s] Combined sparse matrices: %d SNP rows x %d cells", Sys.time(), nrow(ref_all), ncol(ref_all)))
  message(sprintf("[%s] Chromosome VarTrix phase done in %.2f min", Sys.time(), as.numeric(difftime(Sys.time(), prep_t0, units = "mins"))))

  message(sprintf("[%s] Building raw counts from fragments ...", Sys.time()))
  raw_counts <- build_raw_counts(
    fragments_file = config$fragments_file,
    barcodes = barcodes,
    size_table = size_table,
    bin_size = config$bin_size
  )
  raw_counts <- raw_counts[, colnames(raw_counts) %in% barcodes, drop = FALSE]
  raw_counts <- raw_counts[, match(barcodes, colnames(raw_counts)), drop = FALSE]

  seg_table <- build_chromosome_seg_table(size_table, var_all)

  message(sprintf("[%s] Writing SNP matrices and VCF ...", Sys.time()))
  Matrix::writeMM(ref_all, file.path(config$output_dir, "ref_all.mtx"))
  Matrix::writeMM(alt_all, file.path(config$output_dir, "alt_all.mtx"))
  write_vcf(var_all, file.path(config$output_dir, "var_all.vcf"), sample_name = config$sample_name)
  message(sprintf("[%s] Writing raw counts matrix ...", Sys.time()))
  write_counts_tsv(
    raw_counts,
    file.path(config$output_dir, sprintf("chr%sk_fragments.tsv", config$bin_size %/% 1000L))
  )
  saveRDS(seg_table, file.path(config$output_dir, "seg_table.rds"))
  message(sprintf("[%s] Alleloscope prep complete for sample '%s'", Sys.time(), config$sample_name))

  invisible(
    list(
      output_dir = config$output_dir,
      barcodes = file.path(config$output_dir, "barcodes.tsv"),
      ref_all = file.path(config$output_dir, "ref_all.mtx"),
      alt_all = file.path(config$output_dir, "alt_all.mtx"),
      var_all = file.path(config$output_dir, "var_all.vcf"),
      raw_counts = file.path(config$output_dir, sprintf("chr%sk_fragments.tsv", config$bin_size %/% 1000L)),
      seg_table = file.path(config$output_dir, "seg_table.rds")
    )
  )
}
