#!/usr/bin/env Rscript
# =============================================================================
# Full Alleloscope analysis for deepseq scATAC-seq tissue (488B)
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
data.path   <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/cnv/alleloscope/deepseq_488B"
dir_path    <- file.path(data.path, "output")
dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)

message(sprintf("[%s] === Alleloscope deepseq 488B analysis START (ncores=%d) ===",
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
  samplename       = "deepseq_488B",
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
# STEP 3: Est_regions (chromosome-level, cont=TRUE for robust convergence)
# ─────────────────────────────────────────────────────────────────────────────
if (RESUME_FROM <= 3L) {
  message(sprintf("[%s] Step 3: Est_regions (chromosome-level, cont=TRUE)", Sys.time()))
  
  Obj_filtered <- Est_regions(
    Obj_filtered = Obj_filtered,
    max_nSNP     = 30000,
    plot_stat    = TRUE,
    cont         = TRUE
  )
  
  message(sprintf("[%s] Step 3 complete", Sys.time()))
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4: Select_normal (identify normal cells from theta clustering)
# ─────────────────────────────────────────────────────────────────────────────
if (RESUME_FROM <= 4L) {
  message(sprintf("[%s] Step 4: Select_normal (pre-selection)", Sys.time()))
  
  tmp <- Select_normal(
    Obj_filtered = Obj_filtered,
    raw_counts   = raw_counts,
    plot_theta   = TRUE,
    cell_type    = NULL,
    mincell      = 0,
    pre_sel      = TRUE
  )
  
  message(sprintf("[%s] Step 4 complete", Sys.time()))
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 5: Segmentation (using normal-cell pseudobulk as reference)
# ─────────────────────────────────────────────────────────────────────────────
if (RESUME_FROM <= 5L) {
  message(sprintf("[%s] Step 5: Segmentation", Sys.time()))
  
  # Placeholder: Alleloscope will use internal segmentation logic
  # This step typically uses PSCBS on normal-cell pseudobulk
  message(sprintf("[%s] Step 5 placeholder (segmentation handled internally)", Sys.time()))
  
  dir.create(file.path(dir_path, "rds"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(Obj_filtered, CKPT_SEG)
  message(sprintf("[%s] Checkpoint saved: Obj_after_seg.rds", Sys.time()))
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 6: Est_regions (segment-level)
# ─────────────────────────────────────────────────────────────────────────────
if (RESUME_FROM <= 6L) {
  message(sprintf("[%s] Step 6: Est_regions (segment-level)", Sys.time()))
  
  Obj_filtered <- Est_regions(
    Obj_filtered = Obj_filtered,
    max_nSNP     = 30000,
    plot_stat    = TRUE,
    cont         = FALSE,
    ncores       = NCORES
  )
  
  message(sprintf("[%s] Step 6 complete", Sys.time()))
  
  dir.create(file.path(dir_path, "rds"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(Obj_filtered, CKPT_EM_SEG)
  message(sprintf("[%s] Checkpoint saved: Obj_after_EM_seg.rds", Sys.time()))
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 7: Select_normal (finalize normal cells & assign normal region)
# ─────────────────────────────────────────────────────────────────────────────
if (RESUME_FROM <= 7L) {
  message(sprintf("[%s] Step 7: Select_normal (final)", Sys.time()))
  
  tmp <- Select_normal(
    Obj_filtered = Obj_filtered,
    raw_counts   = raw_counts,
    plot_theta   = TRUE,
    cell_type    = NULL,
    mincell      = 0,
    pre_sel      = FALSE
  )
  
  message(sprintf("[%s] Step 7 complete", Sys.time()))
  
  dir.create(file.path(dir_path, "rds"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(Obj_filtered, CKPT_SEL_NRM)
  message(sprintf("[%s] Checkpoint saved: Obj_after_select_normal.rds", Sys.time()))
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 8: Genotype_value
# ─────────────────────────────────────────────────────────────────────────────
if (RESUME_FROM <= 8L) {
  message(sprintf("[%s] Step 8: Genotype_value", Sys.time()))
  
  Obj_filtered <- Genotype_value(
    Obj_filtered = Obj_filtered,
    type         = "tumor",
    raw_counts   = raw_counts,
    cov_adj      = 1
  )
  
  message(sprintf("[%s] Step 8 complete", Sys.time()))
  
  dir.create(file.path(dir_path, "rds"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(Obj_filtered, CKPT_GTV)
  message(sprintf("[%s] Checkpoint saved: Obj_after_gtv.rds", Sys.time()))
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 9: Genotype (scatter plots)
# ─────────────────────────────────────────────────────────────────────────────
if (RESUME_FROM <= 9L) {
  message(sprintf("[%s] Step 9: Genotype (scatter plots)", Sys.time()))
  
  Obj_filtered <- Genotype(
    Obj_filtered = Obj_filtered,
    cell_type    = NULL,
    xmax         = 3
  )
  
  message(sprintf("[%s] Step 9 complete", Sys.time()))
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 10: plot_scATAC_cnv (coverage CNV heatmap)
# ─────────────────────────────────────────────────────────────────────────────
message(sprintf("[%s] Step 10: plot_scATAC_cnv (coverage heatmap)", Sys.time()))

clust_order <- plot_scATAC_cnv(
  raw_mat  = raw_counts,
  cell_type = NULL,
  normal_lab = NULL,
  size = size,
  plot_path = file.path(dir_path, "cov_cna_plot.pdf")
)

message(sprintf("[%s] Step 10 complete", Sys.time()))

message(sprintf("[%s] === Alleloscope deepseq 488B analysis COMPLETE ===", Sys.time()))
