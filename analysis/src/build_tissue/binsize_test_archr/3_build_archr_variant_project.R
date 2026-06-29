#!/usr/bin/env Rscript
# ============================================================================
# build_archr_variant_project_v2.R
#
# Fixed version of build_archr_variant_project.R. The original loaded the
# arrow file as-is (no QC beyond arrow-creation-time minTSS/minFrags), which
# is why deepseq_488B_5000bp_binarizeFALSE ended up with only 7,841 cells
# instead of the ~11k spots the main pipeline (0_create_archr_qc_cluster.R)
# produces for the same tissue.
#
# This version, after loading the (permissively-created, see
# create_arrow_variants_v2.R) arrow file, replicates the main pipeline's
# cell-level QC before running LSI/Clusters/UMAP:
#   1. Filter to the no_edge_effect barcode whitelist
#   2. TSS >= 3 & nFrags >= 1000
#   3. Compute high-nFrags ("doublet", Q3 + 1.5*IQR on nFrags) stats and save
#      them to a CSV for visibility, but do NOT remove any cells based on
#      this - each spot here legitimately covers many cells, so a high
#      fragment count is expected, not evidence of a doublet.
#
# Usage:
#   Rscript build_archr_variant_project_v2.R <tissue> <tilesize> <binarize>
#   Example: Rscript build_archr_variant_project_v2.R deepseq_488B 5000 FALSE
# ============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

log_msg <- function(tag, msg) {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), tag, msg))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  log_msg("error", "Usage: Rscript build_archr_variant_project_v2.R <tissue> <tilesize> <binarize>")
  stop("Missing arguments")
}

tissue <- args[1]
tilesize <- as.integer(args[2])
binarize <- as.logical(args[3])

log_msg("start", sprintf("===== Building ArchR variant project (v2, QC-corrected): %s, %dbp, binarize=%s =====",
                         tissue, tilesize, binarize))

set.seed(42)
addArchRGenome("hg38")
addArchRThreads(threads = 8)

min_tss <- 3
min_frags <- 1000

project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
arrow_dir <- if (binarize) "arrow_binarize_v2" else "arrow_not_binarize_v2"
arrow_file <- file.path(project_root, "Data/01_inputs/arrow", arrow_dir,
                        sprintf("%s_%dbp.arrow", tissue, tilesize))

barcode_file <- file.path(project_root, "Data/01_inputs/barcodes/tissue_barcodes", tissue,
                          sprintf("%s.no_edge_effect.barcodes.tsv", tissue))

output_base_dir <- file.path(project_root, "analysis/binsize_comparison/archr_projects")
output_proj_dir <- file.path(output_base_dir,
                             sprintf("%s_%dbp_binarize%s_v2", tissue, tilesize, binarize))

comparison_dir <- file.path(project_root, "analysis/binsize_comparison")
umap_scores_dir <- file.path(comparison_dir, "umap_scores")
metrics_output_dir <- file.path(comparison_dir, "metrics")

dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(umap_scores_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metrics_output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(arrow_file)) {
  log_msg("error", sprintf("Arrow file not found: %s", arrow_file))
  stop("Arrow file missing - run create_arrow_variants_v2.R first")
}
if (!file.exists(barcode_file)) {
  log_msg("error", sprintf("no_edge_effect barcode file not found: %s", barcode_file))
  stop("Barcode whitelist missing")
}

# ----------------------------------------------------------------------------
# 1. Load ArchR project from the permissively-created arrow file
# ----------------------------------------------------------------------------

log_msg("step", sprintf("Loading arrow file: %s", arrow_file))
proj <- ArchRProject(
  ArrowFiles = arrow_file,
  outputDirectory = output_proj_dir,
  copyArrows = TRUE
)
log_msg("step", sprintf("Loaded ArchR project with %d cells (pre-QC)", ncol(proj)))

# ----------------------------------------------------------------------------
# 2. Filter to no_edge_effect barcode whitelist
# ----------------------------------------------------------------------------

read_barcodes <- function(path) {
  bcs <- readLines(path)
  bcs <- sub("-1$", "", bcs)
  bcs <- bcs[nzchar(bcs)]
  unique(bcs)
}

normalize_barcode <- function(barcode) {
  bc <- as.character(barcode)
  bc <- sub("-1$", "", bc)
  bc <- sub("^.*#", "", bc)
  bc
}

barcodes <- read_barcodes(barcode_file)
log_msg("step", sprintf("Read %d no_edge_effect barcodes", length(barcodes)))

all_cells <- getCellNames(proj)
all_cells_norm <- normalize_barcode(all_cells)
matched_cells <- all_cells[all_cells_norm %in% barcodes]
log_msg("step", sprintf("Matched %d/%d barcodes to cells in arrow", length(matched_cells), length(barcodes)))
if (length(matched_cells) == 0) stop("No matching cells found between arrow and barcode whitelist!")

proj <- proj[matched_cells, ]
log_msg("step", sprintf("After edge-effect filter: %d cells", ncol(proj)))

# ----------------------------------------------------------------------------
# 3. QC: TSS >= 3 & nFrags >= 1000
# ----------------------------------------------------------------------------

before_qc <- ncol(proj)
qc_pass <- proj$TSSEnrichment >= min_tss & proj$nFrags >= min_frags
proj <- proj[qc_pass, ]
log_msg("step", sprintf("After TSS/nFrags QC: %d/%d cells pass (%.1f%% retained)",
                        ncol(proj), before_qc, 100 * ncol(proj) / before_qc))

# ----------------------------------------------------------------------------
# 4. Compute high-nFrags ("doublet") stats ONLY - do NOT remove any cells.
#    Each spot here can legitimately contain many cells, so a high fragment
#    count is expected and is not evidence of a doublet the way it would be
#    in single-cell data. We still report the stats for visibility.
# ----------------------------------------------------------------------------

nfrags <- proj$nFrags  # raw fragment counts
q3 <- quantile(nfrags, 0.75, na.rm = TRUE)
iqr <- IQR(nfrags, na.rm = TRUE)
upper_cutoff <- q3 + 1.5 * iqr
would_flag <- nfrags > upper_cutoff
n_would_flag <- sum(would_flag)

doublet_stats <- data.frame(
  cell = getCellNames(proj),
  nFrags = nfrags,
  high_nFrags_outlier = would_flag
)
doublet_stats_file <- file.path(metrics_output_dir,
                                sprintf("%s_%dbp_binarize%s_v2_high_nfrags_stats.csv", tissue, tilesize, binarize))
write.csv(doublet_stats, doublet_stats_file, row.names = FALSE)

log_msg("step", sprintf(
  "High-nFrags stats (cutoff=%.0f frags): %d/%d cells (%.1f%%) WOULD be flagged - NOT removing them (saved to %s)",
  upper_cutoff, n_would_flag, ncol(proj), 100 * n_would_flag / ncol(proj), doublet_stats_file
))

log_msg("step", sprintf("Final QC'd cell count (no doublet removal): %d", ncol(proj)))

# ----------------------------------------------------------------------------
# 5. Downstream pipeline: LSI -> Clusters -> ImputeWeights -> UMAP (unchanged)
# ----------------------------------------------------------------------------

log_msg("step", "Running addIterativeLSI...")
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  dimsToUse = 1:30,
  varFeatures = 25000,
  clusterParams = list(
    resolution = 0.2,
    sampleCells = 10000,
    n.start = 10
  ),
  force = TRUE
)

log_msg("step", "Running addClusters...")
proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,
  force = TRUE
)

log_msg("step", "Computing imputation weights...")
proj <- addImputeWeights(proj, reducedDims = "IterativeLSI")

log_msg("step", "Running addUMAP...")
proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)

log_msg("step", "Saving ArchR project...")
proj <- saveArchRProject(proj)

log_msg("done", sprintf("Completed successfully (v2): %s %dbp binarize=%s, final n=%d cells",
                        tissue, tilesize, binarize, ncol(proj)))
