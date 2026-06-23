#!/usr/bin/env Rscript
# ============================================================================
# build_archr_variant_project.R
#
# Build ArchR project from an existing arrow file, run the full downstream
# pipeline (LSI -> Clusters -> UMAP), and generate comparison artifacts
# (UMAP+gene-score CSV and PDF plots).
#
# Per-variant worker script, called via qsub with args: tissue tilesize binarize
#
# Usage:
#   Rscript build_archr_variant_project.R <tissue> <tilesize> <binarize>
#   Example: Rscript build_archr_variant_project.R lowseq_489 500 FALSE
#
# ============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

# Logging function
log_msg <- function(tag, msg) {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), tag, msg))
}

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  log_msg("error", "Usage: Rscript build_archr_variant_project.R <tissue> <tilesize> <binarize>")
  stop("Missing arguments")
}

tissue <- args[1]
tilesize <- as.integer(args[2])
binarize <- as.logical(args[3])

log_msg("start", sprintf("===== Building ArchR variant project: %s, %dbp, binarize=%s =====",
                         tissue, tilesize, binarize))

# MUST SET THESE BEFORE CREATING ArchRProject
set.seed(42)
addArchRGenome("hg38")
addArchRThreads(threads = 8)

# Setup paths
project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
arrow_dir <- if (binarize) "arrow_binarize" else "arrow_not_binarize"
arrow_file <- file.path(project_root, "Data/01_inputs/arrow", arrow_dir,
                        sprintf("%s_%dbp.arrow", tissue, tilesize))

output_base_dir <- file.path(project_root, "analysis/binsize_comparison/archr_projects")
output_proj_dir <- file.path(output_base_dir,
                             sprintf("%s_%dbp_binarize%s", tissue, tilesize, binarize))

comparison_dir <- file.path(project_root, "analysis/binsize_comparison")

dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)

# Verify arrow file exists
if (!file.exists(arrow_file)) {
  log_msg("error", sprintf("Arrow file not found: %s", arrow_file))
  stop("Arrow file missing")
}

log_msg("step", sprintf("Loading arrow file: %s", arrow_file))

# Load project from arrow file with copyArrows=TRUE to protect the original
tryCatch({
  proj <- ArchRProject(
    ArrowFiles = arrow_file,
    outputDirectory = output_proj_dir,
    copyArrows = TRUE
  )
  log_msg("step", sprintf("Loaded ArchR project from %s", arrow_file))
}, error = function(e) {
  log_msg("error", sprintf("Failed to load arrow: %s", e$message))
  stop(e)
})

# House parameters (consistent with existing pipeline)
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

# Define unified marker genes (union of existing panels)
marker_genes <- c(
  "CD19", "MS4A1", "TERT", "FOXA1", "SOX17", "HOXD9", "KLRC1", "GNLY", "TPSAB1",
  "CD34", "GATA1", "PAX5", "MME", "CD14", "MPO", "CD3D", "CD8A", "ESR1", "ERBB2", "PGR"
)

log_msg("step", sprintf("Extracting UMAP and gene scores (%d marker genes)...", length(marker_genes)))

# Get UMAP coordinates
umap_df <- data.frame(
  cellID = rownames(proj@cellColData),
  Clusters = proj@cellColData$Clusters,
  UMAP_1 = proj@embeddings$UMAP$df[, 1],
  UMAP_2 = proj@embeddings$UMAP$df[, 2],
  stringsAsFactors = FALSE
)

# Get GeneScoreMatrix and impute
tryCatch({
  gexp <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")

  if (!is.null(gexp)) {
    gexp_mat <- assay(gexp)

    # Fix missing rownames: get gene annotations from ArchR genome
    gene_names <- rownames(gexp)
    if (is.null(gene_names) || all(is.na(gene_names))) {
      # Rownames are missing; get them from the genome annotation
      tryCatch({
        ga <- getGeneAnnotation(proj)
        gene_ranges <- granges(ga)
        gene_names_from_ga <- mcols(gene_ranges)$symbol
        if (!is.null(gene_names_from_ga) && length(gene_names_from_ga) == nrow(gexp_mat)) {
          gene_names <- gene_names_from_ga
          log_msg("step", "Restored gene names from ArchR genome annotation")
        }
      }, error = function(e) {
        log_msg("warn", sprintf("Could not restore gene names: %s", e$message))
      })
    }

    # Identify which marker genes are in the matrix
    available_markers <- intersect(marker_genes, gene_names)
    log_msg("step", sprintf("Found %d/%d marker genes in GeneScoreMatrix",
                            length(available_markers), length(marker_genes)))

    # Get imputation weights and impute gene scores
    impute_weights <- getImputeWeights(proj)

    for (gene in available_markers) {
      raw_scores <- gexp_mat[gene, ]
      # Impute across cells using ArchR's imputation weights
      # imputeMatrix signature: imputeMatrix(mat, imputeWeights, ...)
      imputed_scores <- imputeMatrix(
        Matrix::Matrix(raw_scores, nrow = 1),
        impute_weights
      )
      umap_df[[gene]] <- as.numeric(imputed_scores)
    }
  } else {
    log_msg("warn", "GeneScoreMatrix not found in project; skipping gene score extraction")
  }
}, error = function(e) {
  log_msg("warn", sprintf("Failed to extract gene scores: %s", e$message))
})

# Write CSV
csv_file <- file.path(comparison_dir,
                      sprintf("%s_%dbp_binarize%s_umap_genescores.csv",
                              tissue, tilesize, binarize))
write.csv(umap_df, csv_file, row.names = FALSE)
log_msg("step", sprintf("Saved UMAP+gene-scores CSV: %s", csv_file))

# Generate plots
log_msg("step", "Generating PDF plots...")
pdf_file <- file.path(comparison_dir,
                      sprintf("%s_%dbp_binarize%s_plots.pdf",
                              tissue, tilesize, binarize))

pdf(pdf_file, width = 14, height = 12, onefile = TRUE)

# Page 1: UMAP colored by clusters
tryCatch({
  p_clusters <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "cellColData",
    name = "Clusters",
    embedding = "UMAP",
    size = 1.5
  )
  print(p_clusters)
}, error = function(e) {
  log_msg("warn", sprintf("Failed to plot clusters: %s", e$message))
})

# Pages 2+: Gene scores (one gene per plot)
# Try to get available markers again in case the extraction failed
if (!exists("available_markers")) {
  available_markers <- c()
}

if (length(available_markers) > 0) {
  tryCatch({
    impute_weights <- getImputeWeights(proj)

    for (gene in available_markers) {
      tryCatch({
        p <- plotEmbedding(
          ArchRProj = proj,
          colorBy = "GeneScoreMatrix",
          name = gene,
          embedding = "UMAP",
          imputeWeights = impute_weights,
          size = 1.5,
          quantCut = c(0.01, 0.95)
        )
        print(p)
      }, error = function(e) {
        log_msg("warn", sprintf("Failed to plot gene %s: %s", gene, e$message))
      })
    }
  }, error = function(e) {
    log_msg("warn", sprintf("Failed to generate gene score plots: %s", e$message))
  })
}

dev.off()
log_msg("step", sprintf("Saved plots PDF: %s", pdf_file))

log_msg("done", sprintf("Completed successfully: %s %dbp binarize=%s", tissue, tilesize, binarize))
