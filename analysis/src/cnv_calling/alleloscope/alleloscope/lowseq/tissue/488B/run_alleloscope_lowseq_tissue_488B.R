#!/usr/bin/env Rscript
# =============================================================================
# Full Alleloscope analysis for lowseq scATAC-seq (no matched DNA)
# Reference: https://github.com/seasoncloud/Alleloscope SU008 scATAC vignette
#
# Steps:
#   1. Load prepared inputs & create Alleloscope object
#   2. Matrix_filter
#   3. Est_regions (parallelized, cont=TRUE) — chromosome level
#   4. Select_normal (pre_sel=TRUE) — estimate normal cells from theta clustering
#   5. Segmentation using normal-cell pseudobulk as reference
#   6. Est_regions (parallelized) — segment level
#   7. Select_normal (pre_sel=FALSE) — finalize normal cells & normal regions
#   8. Genotype_value — compute (rho_hat, theta_hat) per region per cell
#   9. Genotype — scatter plot per region
#  10. plot_scATAC_cnv — Step-6 smoothed coverage CNV heatmap
# =============================================================================

# Alleloscope is installed at a non-default library path
.libPaths(c("/projectnb/paxlab/presh/Rlibs/4.5", .libPaths()))

suppressPackageStartupMessages({
  library(Alleloscope)
  library(Matrix)
  library(parallel)
})

NCORES <- as.integer(Sys.getenv("NSLOTS", "8"))

allelo.path <- "/projectnb/paxlab/presh/software/Alleloscope"
data.path   <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/alleloscope/lowseq_tissue_from_existing/488B"
dir_path    <- file.path(data.path, "output")
dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)

message(sprintf("[%s] === Alleloscope lowseq analysis START (ncores=%d) ===",
                Sys.time(), NCORES))

# ─────────────────────────────────────────────────────────────────────────────
# CHECKPOINT RESUME: load latest available saved state and skip completed steps
# ─────────────────────────────────────────────────────────────────────────────
CKPT_SEG     <- file.path(dir_path, "rds", "Obj_after_seg.rds")
CKPT_EM_SEG  <- file.path(dir_path, "rds", "Obj_after_EM_seg.rds")
CKPT_SEL_NRM <- file.path(dir_path, "rds", "Obj_after_select_normal.rds")
CKPT_GTV     <- file.path(dir_path, "rds", "Obj_after_gtv.rds")

# Always load raw_counts + size (small, needed by downstream steps)
# Determine from which step to resume
if (file.exists(CKPT_GTV)) {
  message(sprintf("[%s] [resume] Loading checkpoint: Obj_after_gtv.rds (resuming from Step 9)", Sys.time()))
  size <- read.table(file.path(allelo.path, "data-raw/sizes.cellranger-GRCh38-1.0.0.txt"), stringsAsFactors = FALSE)[1:22, ]
  raw_counts <- read.table(file.path(data.path, "chr1000k_fragments.tsv"), sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  colnames(raw_counts) <- gsub("[.]", "-", colnames(raw_counts))
  Obj_filtered <- readRDS(CKPT_GTV)
  RESUME_FROM <- 9L
} else if (file.exists(CKPT_SEL_NRM)) {
  message(sprintf("[%s] [resume] Loading checkpoint: Obj_after_select_normal.rds (resuming from Step 8)", Sys.time()))
  size <- read.table(file.path(allelo.path, "data-raw/sizes.cellranger-GRCh38-1.0.0.txt"), stringsAsFactors = FALSE)[1:22, ]
  raw_counts <- read.table(file.path(data.path, "chr1000k_fragments.tsv"), sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  colnames(raw_counts) <- gsub("[.]", "-", colnames(raw_counts))
  Obj_filtered <- readRDS(CKPT_SEL_NRM)
  RESUME_FROM <- 8L
} else if (file.exists(CKPT_EM_SEG)) {
  message(sprintf("[%s] [resume] Loading checkpoint: Obj_after_EM_seg.rds (resuming from Step 7)", Sys.time()))
  size <- read.table(file.path(allelo.path, "data-raw/sizes.cellranger-GRCh38-1.0.0.txt"), stringsAsFactors = FALSE)[1:22, ]
  raw_counts <- read.table(file.path(data.path, "chr1000k_fragments.tsv"), sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  colnames(raw_counts) <- gsub("[.]", "-", colnames(raw_counts))
  Obj_filtered <- readRDS(CKPT_EM_SEG)
  RESUME_FROM <- 7L
} else if (file.exists(CKPT_SEG)) {
  message(sprintf("[%s] [resume] Loading checkpoint: Obj_after_seg.rds (resuming from Step 6)", Sys.time()))
  size <- read.table(file.path(allelo.path, "data-raw/sizes.cellranger-GRCh38-1.0.0.txt"), stringsAsFactors = FALSE)[1:22, ]
  raw_counts <- read.table(file.path(data.path, "chr1000k_fragments.tsv"), sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  colnames(raw_counts) <- gsub("[.]", "-", colnames(raw_counts))
  Obj_filtered <- readRDS(CKPT_SEG)
  RESUME_FROM <- 6L
} else {
  RESUME_FROM <- 1L
}

if (RESUME_FROM > 1L) {
  message(sprintf("[%s] Skipping Steps 1-%d (already completed)", Sys.time(), RESUME_FROM - 1L))
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Load inputs
# ─────────────────────────────────────────────────────────────────────────────
if (RESUME_FROM <= 1L) {
message(sprintf("[%s] Step 1: Loading inputs", Sys.time()))

size <- read.table(
  file.path(allelo.path, "data-raw/sizes.cellranger-GRCh38-1.0.0.txt"),
  stringsAsFactors = FALSE
)
size <- size[1:22, ]

barcodes  <- read.table(file.path(data.path, "barcodes.tsv"),
                        sep = "\t", stringsAsFactors = FALSE, header = FALSE)
alt_all   <- readMM(file.path(data.path, "alt_all.mtx"))
ref_all   <- readMM(file.path(data.path, "ref_all.mtx"))
var_all   <- read.table(file.path(data.path, "var_all.vcf"),
                        header = FALSE, sep = "\t", stringsAsFactors = FALSE)

raw_counts <- read.table(
  file.path(data.path, "chr1000k_fragments.tsv"),
  sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE
)

# fix column name encoding (dots -> dashes sometimes happen with read.table)
colnames(raw_counts) <- gsub("[.]", "-", colnames(raw_counts))

message(sprintf("  barcodes: %d | SNPs: %d | bins: %d x %d",
                nrow(barcodes), nrow(var_all), nrow(raw_counts), ncol(raw_counts)))

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Create object & filter
# ─────────────────────────────────────────────────────────────────────────────
message(sprintf("[%s] Step 2: Createobj + Matrix_filter", Sys.time()))


Obj <- Createobj(
  alt_all          = alt_all,
  ref_all          = ref_all,
  var_all          = var_all,
  samplename       = "lowseq_488B",
  genome_assembly  = "GRCh38",
  dir_path         = dir_path,
  barcodes         = barcodes,
  size             = size,
  assay            = "scATACseq"
)

Obj_filtered <- Matrix_filter(
  Obj         = Obj,
  cell_filter = 5,
  SNP_filter  = 5,
  min_vaf     = 0.1,
  max_vaf     = 0.9
)

message(sprintf("  After filter: %d cells, %d SNPs",
                ncol(Obj_filtered$alt_all), nrow(Obj_filtered$alt_all)))
} # end RESUME_FROM <= 1

# ─────────────────────────────────────────────────────────────────────────────
# Parallel wrapper for Est_regions
#
# Drop-in replacement for Est_regions() that uses mclapply to run chromosome/
# segment EM jobs in parallel. Results are stored in the same format that all
# downstream Alleloscope functions expect: Obj_filtered$rds_list and
# Obj_filtered$seg_table_filtered.
#
# Notes:
#  * Alleloscope stores `size` internally as a NAMED NUMERIC VECTOR
#    (names = chromosome numbers without "chr", e.g. "1","2",...,"22").
#  * RDS files are saved as chr1.rds, chr2.rds, etc. (matching Est_regions).
#  * When seg_table_filtered is NULL the same fallback table that Est_regions
#    creates is used so file names are identical.
# ─────────────────────────────────────────────────────────────────────────────
Est_regions_parallel <- function(Obj_filtered,
                                 max_nSNP  = 30000,
                                 plot_stat = TRUE,
                                 min_cell  = 5,
                                 min_snp   = 0,
                                 cont      = FALSE,
                                 max_iter  = 50,
                                 ncores    = 4) {

  # --- replicate Est_regions internal setup --------------------------------
  assay      <- Obj_filtered$assay
  size       <- Obj_filtered$size          # named numeric vector: "1"->248956422
  samplename <- Obj_filtered$samplename
  plot_path  <- file.path(Obj_filtered$dir_path, "plots")
  rds_path   <- file.path(Obj_filtered$dir_path, "rds")

  if (min_snp == 0) min_snp <- Obj_filtered$SNP_filter

  # Build seg_table if not yet set (mirrors Est_regions fallback exactly)
  filtered_seg_table <- Obj_filtered$seg_table_filtered
  if (is.null(filtered_seg_table)) {
    message("Est_regions_parallel: no seg_table_filtered -- using whole chromosomes")
    filtered_seg_table <- data.frame(
      chr     = names(size),
      start   = rep(0, length(size)),
      end     = as.numeric(size),
      states  = 0,
      length  = as.numeric(size),
      mean    = 0,
      var     = 0,
      Var1    = seq_along(size),
      Freq    = 50000,
      chrr    = names(size),
      stringsAsFactors = FALSE
    )
    Obj_filtered$seg_table_filtered <- filtered_seg_table
  }

  dir.create(rds_path, showWarnings = FALSE)
  dir.create(file.path(rds_path,  "EMresults"), showWarnings = FALSE)
  dir.create(file.path(plot_path, "EMresults"), showWarnings = FALSE)
  em_rds_path  <- file.path(rds_path, "EMresults")
  em_plot_path <- file.path(plot_path, "EMresults")

  # --- pre-compute SNP overlaps once (shared read, each worker only reads) --
  var_list <- Obj_filtered$var_all
  var_chr  <- as.numeric(as.character(var_list[, 1]))
  var_pos  <- as.numeric(as.character(var_list[, 2]))

  query   <- GenomicRanges::GRanges(
    paste0("chr", filtered_seg_table$chr),
    IRanges::IRanges(as.numeric(filtered_seg_table$start) + 1,
                     as.numeric(filtered_seg_table$end))
  )
  subject <- GenomicRanges::GRanges(
    paste0("chr", var_chr),
    IRanges::IRanges(var_pos, var_pos)
  )
  ov <- as.matrix(GenomicRanges::findOverlaps(query, subject))

  selseg <- as.character(filtered_seg_table$chrr)

  # --- worker function for one segment -------------------------------------
  process_one <- function(chrr) {
    out_rds <- file.path(em_rds_path, paste0("chr", chrr, ".rds"))

    if (cont && file.exists(out_rds)) {
      cached <- tryCatch(readRDS(out_rds), error = function(e) NULL)
      if (!is.null(cached) && is.list(cached) && !is.null(cached$theta_hat)) {
        message(sprintf("  [skip-cont] chr%s", chrr))
        return(cached)
      }
      message(sprintf("  [warn] chr%s cached RDS invalid or corrupt, re-running", chrr))
    }

    seg_idx  <- which(filtered_seg_table$chrr == chrr)
    chr_ind  <- ov[ov[, 1] %in% seg_idx, 2, drop = TRUE]

    if (length(chr_ind) == 0) {
      message(sprintf("  [skip-no-SNPs] chr%s", chrr))
      return(NULL)
    }

    alt_sub   <- Obj_filtered$alt_all[chr_ind, , drop = FALSE]
    total_sub <- Obj_filtered$total_all[chr_ind, , drop = FALSE]
    var_sub   <- var_list[chr_ind, ]

    # cell filter
    cc_ind <- which(Matrix::colSums(total_sub) > min_cell)
    if (length(cc_ind) == 0) {
      message(sprintf("  [skip-no-cells] chr%s", chrr))
      return(NULL)
    }
    alt_sub   <- alt_sub[, cc_ind, drop = FALSE]
    total_sub <- total_sub[, cc_ind, drop = FALSE]

    # SNP filter
    rr_keep <- which(Matrix::rowSums(total_sub) >= min_snp)
    if (length(rr_keep) == 0) {
      message(sprintf("  [skip-no-SNPs-filter] chr%s", chrr))
      return(NULL)
    }
    if (length(rr_keep) > max_nSNP) rr_keep <- sort(sample(rr_keep, max_nSNP))

    alt_sub   <- alt_sub[rr_keep, , drop = FALSE]
    total_sub <- total_sub[rr_keep, , drop = FALSE]
    var_sub   <- var_sub[rr_keep, 1:5, drop = FALSE]

    af <- Matrix::rowSums(alt_sub) / Matrix::rowSums(total_sub)
    af[is.na(af)] <- 0

    if (nrow(alt_sub) < 3 || length(unique(af)) < 3) {
      message(sprintf("  [skip-low-diversity] chr%s", chrr))
      return(NULL)
    }

    # optional stats plot
    if (plot_stat) {
      tryCatch({
        pdf(file.path(em_plot_path,
                      sprintf("statistics_%s_chr%s.pdf", assay, chrr)))
        par(mfrow = c(3, 1))
        hist(Matrix::colSums(total_sub),
             main = sprintf("%s %s chr%s (%d cells x %d SNPs)",
                            samplename, assay, chrr,
                            ncol(total_sub), nrow(total_sub)),
             xlab = "per-cell coverage", breaks = 100)
        hist(Matrix::rowSums(total_sub),
             main = sprintf("chr%s SNP coverage", chrr),
             xlab = "per-SNP coverage", xlim = c(0, 100), breaks = 1000)
        hist(af, 100, main = "Histogram of VAF values")
        dev.off()
      }, error = function(e) {
        try(dev.off(), silent = TRUE)
        message(sprintf("  [warn] stats plot failed for chr%s: %s", chrr, e$message))
      })
    }

    # EM
    result <- EM(
      ref_table = as.matrix(total_sub - alt_sub),
      alt_table = as.matrix(alt_sub),
      seed      = 1000,
      max_iter  = max_iter
    )
    result$barcodes <- colnames(total_sub)
    result$SNPs <- paste0("chr", var_sub$V1, ":", var_sub$V2,
                          "_", var_sub$V4, "_", var_sub$V5)

    saveRDS(result, out_rds)

    tryCatch({
      pdf(file.path(em_plot_path, sprintf("EMresult_chr%s.pdf", chrr)))
      par(mfrow = c(2, 1))
      hist(result$I_hat,     100, xlim = c(0, 1), main = paste0("I_hat chr",     chrr))
      hist(result$theta_hat, 100, xlim = c(0, 1), main = paste0("theta_hat chr", chrr))
      dev.off()
    }, error = function(e) {
      try(dev.off(), silent = TRUE)
    })

    message(sprintf("  [done] chr%s (theta range: %.3f-%.3f)",
                    chrr, min(result$theta_hat), max(result$theta_hat)))
    result
  }

  message(sprintf("[%s] Est_regions_parallel: %d regions on %d cores (cont=%s)",
                  Sys.time(), length(selseg), ncores, cont))

  results_list         <- mclapply(selseg, process_one,
                                   mc.cores       = ncores,
                                   mc.preschedule = FALSE)
  names(results_list)  <- paste0("chr", selseg)

  # Filter out NULL, try-error, and any result that is not a valid EM list
  is_valid_em <- function(x) {
    !is.null(x) && !inherits(x, "try-error") && is.list(x) && !is.null(x$theta_hat)
  }
  invalid_names <- names(results_list)[!vapply(results_list, is_valid_em, logical(1))]
  if (length(invalid_names) > 0) {
    message(sprintf(
      "[%s] [warn] %d workers returned errors or invalid results: %s",
      Sys.time(), length(invalid_names), paste(invalid_names, collapse = ", ")
    ))
  }
  results_list <- Filter(is_valid_em, results_list)

  # Retry missing/failed chromosomes serially so downstream steps do not
  # receive a partial rds_list when some parallel workers fail transiently.
  expected_names <- paste0("chr", selseg)
  missing_names <- setdiff(expected_names, names(results_list))
  if (length(missing_names) > 0) {
    message(sprintf(
      "[%s] [warn] Est_regions_parallel missing %d regions after parallel pass: %s",
      Sys.time(), length(missing_names), paste(missing_names, collapse = ", ")
    ))

    retry_res <- lapply(sub("^chr", "", missing_names), function(chrr) {
      tryCatch(process_one(chrr), error = function(e) {
        message(sprintf("  [error] serial retry failed for chr%s: %s", chrr, e$message))
        NULL
      })
    })
    names(retry_res) <- missing_names
    retry_res <- Filter(is_valid_em, retry_res)
    if (length(retry_res) > 0) {
      results_list <- c(results_list, retry_res)
    }
  }

  # Enforce full chromosome-level EM completion before Select_normal.
  final_names <- names(results_list)
  missing_final <- setdiff(expected_names, final_names)
  if (length(missing_final) > 0) {
    stop(sprintf(
      "EM incomplete: %d/%d regions available; missing %s",
      length(final_names),
      length(expected_names),
      paste(missing_final, collapse = ", ")
    ), call. = FALSE)
  }

  Obj_filtered$rds_list <- results_list
  message(sprintf("[%s] Est_regions_parallel: %d/%d regions completed",
                  Sys.time(), length(results_list), length(selseg)))
  Obj_filtered
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: Est_regions -- chromosome level (parallel, resume via cont=TRUE)
# chr1-chr11 already computed; cont=TRUE skips those and runs chr12-chr22.
# ─────────────────────────────────────────────────────────────────────────────
if (RESUME_FROM <= 3L) {
message(sprintf("[%s] Step 3: Est_regions chr-level (parallel, cont=TRUE)", Sys.time()))
# Cap to 1 core for chr-level EM; mcfork OOM persists for large tissue matrices.
CHR_NCORES <- 1L
message(sprintf("  [info] chr-level EM using %d cores (capped from %d to avoid fork OOM)",
                CHR_NCORES, NCORES))

Obj_filtered <- Est_regions_parallel(
  Obj_filtered = Obj_filtered,
  max_nSNP     = 30000,
  plot_stat    = TRUE,
  cont         = TRUE,
  max_iter     = 50,
  ncores       = CHR_NCORES
)

saveRDS(Obj_filtered,
        file.path(dir_path, "rds", "Obj_after_EM_chr.rds"))
message(sprintf("  Checkpoint saved: rds/Obj_after_EM_chr.rds"))
} # end RESUME_FROM <= 3

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4: Select_normal -- pre_sel=TRUE (no matched DNA needed)
# Identifies candidate normal cells from theta_hat clustering alone.
# The cluster whose theta_hat is closest to 0.5 is labelled normal.
# ─────────────────────────────────────────────────────────────────────────────
if (RESUME_FROM <= 4L) {
message(sprintf("[%s] Step 4: Select_normal (pre_sel=TRUE)", Sys.time()))

Obj_filtered <- Select_normal(
  Obj_filtered = Obj_filtered,
  raw_counts   = raw_counts,
  cell_nclust  = 5,
  pre_sel      = TRUE,
  plot_theta   = TRUE
)

barcode_normal_pre <- Obj_filtered$select_normal$barcode_normal
message(sprintf("  Candidate normal cells: %d", length(barcode_normal_pre)))

saveRDS(barcode_normal_pre,
        file.path(dir_path, "rds", "barcode_normal_pre.rds"))
} # end RESUME_FROM <= 4

# ─────────────────────────────────────────────────────────────────────────────
# STEP 5: Segmentation
# Build a pseudobulk reference from the computationally identified normal cells
# and run HMM segmentation on the tumour sample.
# ─────────────────────────────────────────────────────────────────────────────
if (RESUME_FROM <= 5L) {
message(sprintf("[%s] Step 5: Segmentation", Sys.time()))

barcode_normal_pre <- readRDS(file.path(dir_path, "rds", "barcode_normal_pre.rds"))
normal_cols <- intersect(barcode_normal_pre, colnames(raw_counts))
if (length(normal_cols) < 10) {
  warning(sprintf(
    "Only %d normal barcodes found in raw_counts; using all cells as reference",
    length(normal_cols)
  ))
  normal_cols <- colnames(raw_counts)
}
message(sprintf("  Normal cells for pseudobulk ref: %d", length(normal_cols)))

# Pseudobulk reference: same row structure as raw_counts, one column
ref_counts_bulk <- as.data.frame(
  rowSums(raw_counts[, normal_cols, drop = FALSE])
)
rownames(ref_counts_bulk) <- rownames(raw_counts)
colnames(ref_counts_bulk) <- "normal_pseudobulk"

run_segmentation <- function(ref_counts, label) {
  message(sprintf("  Segmentation attempt using %s", label))
  Segmentation(
    Obj_filtered = Obj_filtered,
    raw_counts   = raw_counts,
    ref_counts   = ref_counts,
    plot_seg     = TRUE
  )
}

Obj_filtered <- tryCatch(
  run_segmentation(ref_counts_bulk, "normal-cell pseudobulk"),
  error = function(e) {
    message(sprintf("  [warn] Segmentation failed (%s); retrying with all-cells pseudobulk", e$message))
    ref_all_bulk <- as.data.frame(rowSums(raw_counts))
    rownames(ref_all_bulk) <- rownames(raw_counts)
    colnames(ref_all_bulk) <- "all_cells_pseudobulk"
    tryCatch(
      run_segmentation(ref_all_bulk, "all-cells pseudobulk (fallback)"),
      error = function(e2) {
        message(sprintf("  [warn] Segmentation fallback also failed (%s); proceeding with existing seg_table_filtered", e2$message))
        Obj_filtered
      }
    )
  }
)

n_segs <- if (is.null(Obj_filtered$seg_table_filtered)) 0L else nrow(Obj_filtered$seg_table_filtered)
message(sprintf("  Segments identified: %d", n_segs))
saveRDS(Obj_filtered,
        file.path(dir_path, "rds", "Obj_after_seg.rds"))
} # end RESUME_FROM <= 5

# ─────────────────────────────────────────────────────────────────────────────
# STEP 6: Est_regions -- segment level (parallel, fresh run)
# Use a separate output subdirectory so segment-level EM files do not
# overwrite the chromosome-level EMresults files.
# ─────────────────────────────────────────────────────────────────────────────
if (RESUME_FROM <= 6L) {
message(sprintf("[%s] Step 6: Est_regions seg-level (parallel, cont=TRUE)", Sys.time()))

seg_dir_path <- file.path(dir_path, "seg")
dir.create(seg_dir_path, showWarnings = FALSE)
dir.create(file.path(seg_dir_path, "plots", "EMresults"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(seg_dir_path, "rds",   "EMresults"), showWarnings = FALSE, recursive = TRUE)

Obj_filtered_seg           <- Obj_filtered
Obj_filtered_seg$dir_path  <- seg_dir_path

# Use half the cores to avoid mcfork OOM on large objects
SEG_NCORES <- 1L
message(sprintf("  [info] seg-level EM using %d cores (half of %d to avoid fork OOM)",
                SEG_NCORES, NCORES))

Obj_filtered_seg <- Est_regions_parallel(
  Obj_filtered = Obj_filtered_seg,
  max_nSNP     = 30000,
  plot_stat    = TRUE,
  cont         = TRUE,
  max_iter     = 50,
  ncores       = SEG_NCORES
)

# Merge results back into the main object
Obj_filtered$rds_list           <- Obj_filtered_seg$rds_list
Obj_filtered$seg_table_filtered <- Obj_filtered_seg$seg_table_filtered

saveRDS(Obj_filtered,
        file.path(dir_path, "rds", "Obj_after_EM_seg.rds"))
message(sprintf("  Checkpoint saved: rds/Obj_after_EM_seg.rds"))
} # end RESUME_FROM <= 6

# ─────────────────────────────────────────────────────────────────────────────
# STEP 7: Select_normal -- finalize with coverage + theta (pre_sel=FALSE)
# ─────────────────────────────────────────────────────────────────────────────
if (RESUME_FROM <= 7L) {
message(sprintf("[%s] Step 7: Select_normal (pre_sel=FALSE, post-seg)", Sys.time()))

Obj_filtered <- Select_normal(
  Obj_filtered = Obj_filtered,
  raw_counts   = raw_counts,
  cell_nclust  = 5,
  pre_sel      = FALSE,
  plot_theta   = TRUE
)

barcode_normal_final <- Obj_filtered$select_normal$barcode_normal
region_normal        <- Obj_filtered$select_normal$region_normal
message(sprintf("  Final normal cells: %d", length(barcode_normal_final)))
message(sprintf("  Normal region(s):   %s",
                paste(region_normal[1:min(3, length(region_normal))], collapse = ", ")))

saveRDS(Obj_filtered,
        file.path(dir_path, "rds", "Obj_after_select_normal.rds"))
} # end RESUME_FROM <= 7

# ─────────────────────────────────────────────────────────────────────────────
# STEP 8: Genotype_value -- compute (rho_hat, theta_hat) per region per cell
# Uses the identified normal region for coverage normalisation.
# ─────────────────────────────────────────────────────────────────────────────
if (RESUME_FROM <= 8L) {
message(sprintf("[%s] Step 8: Genotype_value", Sys.time()))

if (RESUME_FROM == 8L) {
  region_normal <- Obj_filtered$select_normal$region_normal
}
Obj_filtered$ref <- region_normal[1]

Obj_filtered <- Genotype_value(
  Obj_filtered = Obj_filtered,
  type         = "tumor",
  raw_counts   = raw_counts,
  cov_adj      = 1,
  qt_filter    = TRUE,
  cell_filter  = TRUE,
  refr         = TRUE
)

saveRDS(Obj_filtered,
        file.path(dir_path, "rds", "Obj_after_gtv.rds"))
message("  Checkpoint saved: rds/Obj_after_gtv.rds")
} # end RESUME_FROM <= 8

# ─────────────────────────────────────────────────────────────────────────────
# STEP 9: Genotype -- per-region (rho_hat, theta_hat) scatter plots
# ─────────────────────────────────────────────────────────────────────────────
if (RESUME_FROM <= 9L) {
message(sprintf("[%s] Step 9: Genotype (scatter plots)", Sys.time()))

gtype_plot <- file.path(dir_path, "plots",
                        paste0("gtype_scatter_ref_", Obj_filtered$ref, ".pdf"))
gtype_rds  <- file.path(dir_path, "rds", "genotypes.rds")

Obj_filtered <- Genotype(
  Obj_filtered = Obj_filtered,
  plot_path    = gtype_plot,
  rds_path     = gtype_rds,
  cell_type    = NULL,
  legend       = TRUE
)

saveRDS(Obj_filtered,
        file.path(dir_path, "rds", "Obj_final.rds"))
message("  Final object saved: rds/Obj_final.rds")
} # end RESUME_FROM <= 9

# ─────────────────────────────────────────────────────────────────────────────
# STEP 10: plot_scATAC_cnv -- Step-6 smoothed coverage CNV heatmap
#
# Cells are labelled normal/tumor from the computationally identified normals.
# The heatmap smooths coverage in 10 Mb windows across all autosomes and
# clusters cells by correlation (ward.D2), annotated by cell type.
# ─────────────────────────────────────────────────────────────────────────────
message(sprintf("[%s] Step 10: plot_scATAC_cnv (Step-6 heatmap)", Sys.time()))

# Ensure barcode_normal_final is set regardless of resume point
if (!exists("barcode_normal_final")) {
  barcode_normal_final <- Obj_filtered$select_normal$barcode_normal
}

# Re-read full size table for plot_scATAC_cnv (needs 2-column data frame)
size_df <- read.table(
  file.path(allelo.path, "data-raw/sizes.cellranger-GRCh38-1.0.0.txt"),
  stringsAsFactors = FALSE
)
size_df <- size_df[1:22, ]

Obj_filtered <- readRDS(file.path(dir_path, "rds", "Obj_final.rds"))
raw_counts <- read.table(file.path(data.path, "chr1000k_fragments.tsv"),
  sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
colnames(raw_counts) <- gsub("[.]", "-", colnames(raw_counts))
barcode_normal_final <- Obj_filtered$select_normal$barcode_normal


# Build cell_type data frame: barcode | label
all_barcodes <- colnames(raw_counts)
cell_type_df <- data.frame(
  barcode   = all_barcodes,
  cell_type = ifelse(all_barcodes %in% barcode_normal_final, "normal", "tumor"),
  stringsAsFactors = FALSE
)

heatmap_path <- file.path(dir_path, "plots", "step6_CNV_coverage_heatmap.png")

# Remove zero-coverage cells before normalization (0/0 → NaN → rowMedians all NA)
keep_cells  <- colSums(raw_counts[, cell_type_df$barcode]) > 0
message(sprintf("  Cells passing coverage filter: %d / %d",
                sum(keep_cells), length(keep_cells)))
cell_type_df_filt <- cell_type_df[keep_cells, , drop = FALSE]


cov_obj <- plot_scATAC_cnv(
  raw_mat        = as.matrix(raw_counts[, cell_type_df_filt$barcode]),
  cell_type      = cell_type_df_filt,
  normal_lab     = "normal",
  size           = size_df,
  window_w       = 10000000,
  window_step    = 2000000,
  plot_path      = heatmap_path,
  nclust         = 3,
  var.filter     = T
)

saveRDS(cov_obj, file.path(dir_path, "rds", "cov_obj.rds"))

# ───────────────────────────────────────────────────────────────────────────
# Filter: Remove windows with ≥50% blacklist overlap
# ───────────────────────────────────────────────────────────────────────────
message("Step: Filter blacklist regions from raw_counts (based on rownames)")

# Load blacklist file
BLACKLIST_FILE <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/cnv/CNV_results/epianeufinder/hg38-blacklist.v2.bed"

if (!file.exists(BLACKLIST_FILE)) {
  stop(paste("Blacklist file not found:", BLACKLIST_FILE))
}

blacklist_df <- data.table::fread(BLACKLIST_FILE, header=FALSE, 
                                  col.names=c("chrom", "start", "end", "name"))
message(paste("Loaded", nrow(blacklist_df), "blacklist regions"))

# Convert blacklist to GRanges
blacklist_gr <- GenomicRanges::GRanges(
  seqnames=blacklist_df$chrom,
  ranges=IRanges::IRanges(start=blacklist_df$start, end=blacklist_df$end),
  name=blacklist_df$name
)

# Parse window coordinates from ROWNAMES (format: chr1-start-end)
message(paste("Parsing", nrow(raw_counts), "windows from raw_counts rownames"))

window_coords <- data.frame(
  window_name = rownames(raw_counts),
  stringsAsFactors = FALSE
)

parts_list <- strsplit(window_coords$window_name, "-")
window_coords$chr <- sapply(parts_list, function(x) x[1])
window_coords$start <- as.numeric(sapply(parts_list, function(x) x[2]))
window_coords$end <- as.numeric(sapply(parts_list, function(x) x[3]))

# Convert windows to GRanges
windows_gr <- GenomicRanges::GRanges(
  seqnames=window_coords$chr,
  ranges=IRanges::IRanges(start=window_coords$start, end=window_coords$end),
  window_name=window_coords$window_name
)

# Calculate overlap percentage
message("Computing overlap percentages with blacklist")

overlap_pct <- numeric(length(windows_gr))
for (i in seq_along(windows_gr)) {
  hits <- GenomicRanges::findOverlaps(windows_gr[i], blacklist_gr, type = "any")
  if (length(hits) == 0) next
  
  subject_idx <- S4Vectors::subjectHits(hits)
  
  # Repeat window to match number of overlapping blacklist regions
  overlapping <- GenomicRanges::pintersect(
    rep(windows_gr[i], length(subject_idx)), 
    blacklist_gr[subject_idx]
  )
  
  overlap_pct[i] <- sum(GenomicRanges::width(overlapping)) / GenomicRanges::width(windows_gr[i]) * 100
}

# Keep rows with <50% overlap
keep_50pct <- overlap_pct < 50

message(paste("Min overlap:", round(min(overlap_pct, na.rm = TRUE), 2), "%"))
message(paste("Max overlap:", round(max(overlap_pct, na.rm = TRUE), 2), "%"))
message(paste("Mean overlap:", round(mean(overlap_pct, na.rm = TRUE), 2), "%"))
message(paste("Windows retained (<50% overlap):", sum(keep_50pct), "/", length(keep_50pct)))
message(paste("Windows removed (>=50% overlap):", sum(!keep_50pct)))

# Apply filter to raw_counts ROWS
raw_counts_filt <- raw_counts[keep_50pct, , drop = FALSE]

message(paste("Filtered matrix:", nrow(raw_counts_filt), "windows x", ncol(raw_counts_filt), "cells"))
message("Blacklist filtering complete")





heatmap_path <- file.path(dir_path, "plots", "step6_CNV_coverage_heatmap_filtered.png")

# Remove zero-coverage cells before normalization (0/0 → NaN → rowMedians all NA)
keep_cells  <- colSums(raw_counts_filt[, cell_type_df$barcode]) > 0
message(sprintf("  Cells passing coverage filter: %d / %d",
                sum(keep_cells), length(keep_cells)))
cell_type_df_filt <- cell_type_df[keep_cells, , drop = FALSE]


cov_obj <- plot_scATAC_cnv(
  raw_mat        = as.matrix(raw_counts_filt[, cell_type_df_filt$barcode]),
  cell_type      = cell_type_df_filt,
  normal_lab     = "normal",
  size           = size_df,
  window_w       = 10000000,
  window_step    = 2000000,
  plot_path      = heatmap_path,
  nclust         = 3,
  var.filter     = T
)


# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
message(sprintf("[%s] === All steps complete ===", Sys.time()))
message(sprintf("Output directory: %s", dir_path))
message("Key outputs:")
message(sprintf("  Step-6 CNV coverage heatmap : plots/step6_CNV_coverage_heatmap.png"))
message(sprintf("  Genotype scatter plots      : plots/gtype_scatter_ref_%s.pdf",
                Obj_filtered$ref))
message(sprintf("  EM results (chr-level)      : rds/EMresults/"))
message(sprintf("  Final Alleloscope object    : rds/Obj_final.rds"))
message(sprintf("  Normal barcodes (final)     : %d cells",
                length(barcode_normal_final)))


