#!/usr/bin/env Rscript
# 00_export_for_stalign.R
# Export GeneIntegrationMatrix (pseudo-RNA per ATAC cell) + spatial coordinates
# from both tissues for use by Python Tangram and STAlign scripts.
#
# Usage: Rscript 00_export_for_stalign.R [488B|489|both]
# Prerequisites: 01_label_transfer_{488B,489}.R must have been run.

suppressPackageStartupMessages({
  library(ArchR)
  library(Matrix)
  library(dplyr)
})

addArchRThreads(threads = 4)
addArchRGenome("hg38")

args <- commandArgs(trailingOnly = TRUE)
tissues <- if (length(args) > 0 && args[1] %in% c("488B","489")) args[1] else c("488B","489")

CONFIGS <- list(
  "488B" = list(
    archr_dir    = "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration/488B/objects/archr_488B_with_integration",
    scrna_rds    = "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration/488B/objects/seurat_scrna_488B_processed.rds",
    spatial_csv  = "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/cnv/calicoST/deepseq_488B/intermediate/spatial_coords.csv",
    out_dir      = "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration/stalign/488B"
  ),
  "489" = list(
    archr_dir    = "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration/489/objects/archr_489_with_integration",
    scrna_rds    = "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration/489/objects/seurat_scrna_489_processed.rds",
    spatial_csv  = "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/cnv/calicoST/lowseq_489/intermediate/spatial_coords.csv",
    out_dir      = "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration/stalign/489"
  )
)

export_tissue <- function(tissue) {
  cfg <- CONFIGS[[tissue]]
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  cat("\n=== Exporting tissue:", tissue, "===\n")

  if (!dir.exists(cfg$archr_dir)) {
    cat("SKIP:", tissue, "— integrated ArchR project not found at", cfg$archr_dir, "\n")
    return(invisible(NULL))
  }

  proj <- loadArchRProject(cfg$archr_dir, showLogo = FALSE)
  meta <- as.data.frame(getCellColData(proj))

  # ── 1. GeneScoreMatrix (gene accessibility scores, genes × cells) ─────────
  # GeneIntegrationMatrix is only available when ArchR addGeneIntegrationMatrix
  # was run; label transfer was done via Seurat so only GeneScoreMatrix exists.
  cat("Extracting GeneScoreMatrix...\n")
  gs     <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
  gs_mat <- assay(gs)
  if (!is.null(rowData(gs)$name)) {
    rownames(gs_mat) <- make.unique(as.character(rowData(gs)$name))
  }
  cat("GeneScore dims:", nrow(gs_mat), "genes ×", ncol(gs_mat), "cells\n")

  # Write as sparse MTX for Python
  writeMM(gs_mat, file.path(cfg$out_dir, paste0("atac_pseudorna_", tissue, ".mtx")))
  write.csv(data.frame(gene=rownames(gs_mat)),
            file.path(cfg$out_dir, paste0("atac_gene_names_", tissue, ".csv")),
            row.names=FALSE)
  write.csv(data.frame(cell=colnames(gs_mat)),
            file.path(cfg$out_dir, paste0("atac_cell_names_", tissue, ".csv")),
            row.names=FALSE)
  cat("Saved GeneScoreMatrix.\n")

  # ── 2. ATAC spatial coordinates (µm) ──────────────────────────────────────
  # Check for spatial cols in cellColData
  sp_cols <- grep("spatial|x_um|y_um", colnames(meta), ignore.case=TRUE, value=TRUE)
  if (length(sp_cols) >= 2) {
    sp_df <- data.frame(
      cell = rownames(meta),
      x_um = meta[[sp_cols[1]]],
      y_um = meta[[sp_cols[2]]]
    )
  } else {
    cat("Spatial coords not in cellColData, reading from CalicoST spatial_coords.csv...\n")
    coords <- read.csv(cfg$spatial_csv, row.names = 1)  # archr_barcode as rowname
    sp_df <- data.frame(
      cell = colnames(gs_mat),
      x_um = coords[colnames(gs_mat), "x_spatial"],
      y_um = coords[colnames(gs_mat), "y_spatial"]
    )
  }
  write.csv(sp_df, file.path(cfg$out_dir, paste0("atac_spatial_coords_", tissue, ".csv")),
            row.names=FALSE)
  cat("Saved spatial coordinates:", nrow(sp_df), "cells\n")

  # ── 3. ATAC metadata (cell type labels, cluster, prediction score) ─────────
  meta_out <- data.frame(
    cell             = rownames(meta),
    atac_cluster     = meta$Clusters,
    predicted_type   = meta$predictedGroup_co,
    prediction_score = meta$predictedScore_co
  )
  write.csv(meta_out, file.path(cfg$out_dir, paste0("atac_metadata_", tissue, ".csv")),
            row.names=FALSE)

  # ── 4. scRNA-seq normalized matrix for Tangram ────────────────────────────
  cat("Exporting scRNA-seq for Tangram...\n")
  scrna <- readRDS(cfg$scrna_rds)
  scrna_mat <- GetAssayData(scrna, layer = "data")  # log-normalized
  writeMM(scrna_mat, file.path(cfg$out_dir, paste0("scrna_normexpr_", tissue, ".mtx")))
  write.csv(data.frame(gene=rownames(scrna_mat)),
            file.path(cfg$out_dir, paste0("scrna_gene_names_", tissue, ".csv")), row.names=FALSE)
  write.csv(data.frame(cell=colnames(scrna_mat)),
            file.path(cfg$out_dir, paste0("scrna_cell_names_", tissue, ".csv")), row.names=FALSE)
  write.csv(scrna@meta.data[, c("seurat_clusters", "cell_type")],
            file.path(cfg$out_dir, paste0("scrna_metadata_", tissue, ".csv")))
  cat("Saved scRNA data:", nrow(scrna_mat), "genes ×", ncol(scrna_mat), "cells\n")

  cat("Export complete for tissue", tissue, "→", cfg$out_dir, "\n")
}

for (t in tissues) export_tissue(t)
cat("\nAll exports done.\n")
