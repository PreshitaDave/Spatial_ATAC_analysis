#!/usr/bin/env Rscript
# ============================================================================
# 0_create_archr_no_doublet.R
# Create ArchR objects from Arrow files with basic QC (NO DOUBLET REMOVAL),
# cluster, add UMAP, and save comprehensive PDFs
# ============================================================================


suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
  library(gridExtra)
})

project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
threads <- 8
min_tss <- 3
min_frags <- 1000

set.seed(42)
addArchRGenome("hg38")
addArchRThreads(threads = threads)

# Logging function
log_msg <- function(tag, msg) {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), tag, msg))
}

# Setup paths
data_dir <- file.path(project_root, "Data", "01_inputs")
output_root <- file.path(project_root, "Data", "01_outputs")
archr_output <- file.path(output_root, "archR_objects")
arrow_dir <- file.path(data_dir, "arrow")
plot_output <- file.path(project_root, "analysis", "plots", "archr_obj")

dir.create(plot_output, recursive = TRUE, showWarnings = FALSE)

log_msg("start", "===== ArchR Project Creation (NO DOUBLET REMOVAL) =====")

# Define objects and Arrow file paths
objects <- list(
  deepseq_488B = list(
    arrow = file.path(arrow_dir, "Deepseq_488B.arrow"),
    barcodes = file.path(data_dir, "barcodes", "tissue_barcodes", "deepseq_488B", "deepseq_488B.no_edge_effect.barcodes.tsv"),
    sample_name = "Deepseq_488B",
    name = "deepseq_488B"
  ),
  deepseq_489 = list(
    arrow = file.path(arrow_dir, "Deepseq_489.arrow"),
    barcodes = file.path(data_dir, "barcodes", "tissue_barcodes", "deepseq_489", "deepseq_489.no_edge_effect.barcodes.tsv"),
    sample_name = "Deepseq_489",
    name = "deepseq_489"
  ),
  lowseq_488B = list(
    arrow = file.path(arrow_dir, "Lowseq_488B.arrow"),
    barcodes = file.path(data_dir, "barcodes", "tissue_barcodes", "lowseq_488B", "lowseq_488B.no_edge_effect.barcodes.tsv"),
    sample_name = "Lowseq_488B",
    name = "lowseq_488B"
  ),
  lowseq_489 = list(
    arrow = file.path(arrow_dir, "Lowseq_489.arrow"),
    barcodes = file.path(data_dir, "barcodes", "tissue_barcodes", "lowseq_489", "lowseq_489.no_edge_effect.barcodes.tsv"),
    sample_name = "Lowseq_489",
    name = "lowseq_489"
  )
)

# Barcode normalization function
normalize_barcode <- function(x) {
  sub(".*#([A-Za-z0-9]+)-1$", "\\1", x)
}

# Read barcode file
read_barcodes <- function(path) {
  if (!file.exists(path)) return(character(0))
  bcs <- readLines(path)
  bcs <- bcs[nzchar(bcs)]
  unique(bcs)
}

# Filter project to barcodes
filter_to_barcodes <- function(proj, barcode_list) {
  all_cells <- getCellNames(proj)
  all_cells_norm <- normalize_barcode(all_cells)
  matched_idx <- which(all_cells_norm %in% barcode_list)
  matched_cells <- all_cells[matched_idx]
  log_msg("filter", sprintf("Matched %d/%d cells to barcode list", length(matched_cells), length(barcode_list)))
  if (length(matched_cells) == 0) stop("No matching cells found!")
  proj2 = proj[matched_cells, ]
  proj2@cellColData$barcode_mapping = normalize_barcode(getCellNames(proj2))
  return(proj2)
}

# Process each object
for (obj_name in names(objects)) {
  tryCatch({
    obj_info <- objects[[obj_name]]
    log_msg("step", sprintf("--- Processing %s ---", obj_name))
    
    # Verify Arrow file exists
    if (!file.exists(obj_info$arrow)) {
      log_msg("warn", sprintf("Arrow file not found: %s", obj_info$arrow))
      next
    }
    
    # Create ArchR project from Arrow file with output directory set to final location
    log_msg("step", sprintf("Creating ArchR project from Arrow: %s", obj_info$arrow))
    
    # Set output directory for this project
    obj_output_dir <- file.path(archr_output, obj_name)
    dir.create(obj_output_dir, recursive = TRUE, showWarnings = FALSE)
    project_dir <- file.path(obj_output_dir, sprintf("%s_archR_project_final", obj_name))
    
    proj <- ArchRProject(
      ArrowFiles = obj_info$arrow,
      outputDirectory = project_dir,
      copyArrows = FALSE
    )
    proj <- addTileMatrix(proj, force = TRUE, tileSize = 500)
    
    a = getMatrixFromArrow(
      ArrowFile = obj_info$arrow,
      useMatrix = "TileMatrix",
      verbose = TRUE,
      binarize = FALSE,
      logFile = createLogFile("getMatrixFromArrow")
    )
    
    log_msg("step", sprintf("Created project with %d cells at: %s", ncol(proj), project_dir))
    
    # Read and filter to edge-effect-filtered barcodes
    if (file.exists(obj_info$barcodes)) {
      barcodes <- read_barcodes(obj_info$barcodes)
      log_msg("step", sprintf("Filtering to %d edge-effect-filtered barcodes", length(barcodes)))
      proj <- filter_to_barcodes(proj, barcodes)
    }
    
    n_initial <- ncol(proj)
    log_msg("step", sprintf("Cells before QC filtering: %d", n_initial))
    
    # Get metadata
    metadata <- getCellColData(proj)
    tss_vec <- metadata$TSSEnrichment
    nfrags_vec <- metadata$nFrags
    
    # Step 1: Apply basic QC filtering ONLY (TSS >= 3, nFrags >= 1000)
    # NO DOUBLET REMOVAL
    log_msg("step", "Applying basic QC filters (TSS >= 3, nFrags >= 1000)...")
    qc_pass <- (tss_vec >= min_tss) & (nfrags_vec >= min_frags)
    qc_idx <- which(qc_pass)
    n_final <- length(qc_idx)
    
    log_msg("step", sprintf("After basic QC: %d cells (removed %d)", n_final, n_initial - n_final))
    
    # Filter to QC-passing cells
    proj <- proj[qc_idx, ]
    
    # Step 2: Add dimensionality reduction (LSI)
    log_msg("step", "Adding LSI (Latent Semantic Indexing)...")
    proj <- addIterativeLSI(
      ArchRProj = proj,
      useMatrix = "TileMatrix",
      name = "IterativeLSI",
      iterations = 2,
      clusterParams = list(resolution = 0.2, sampleCells = 10000, n.start = 10),
      varFeatures = 25000,
      dimsToUse = 1:30,
      force = TRUE
    )
    log_msg("step", "LSI added")
    
    # Step 3: Add clustering
    log_msg("step", "Adding clustering...")
    proj <- addClusters(
      input = proj,
      reduction = "IterativeLSI",
      method = "Seurat",
      name = "Clusters",
      resolution = 0.8,
      force = TRUE
    )
    n_clusters <- length(unique(proj$Clusters))
    log_msg("step", sprintf("Added clustering: %d clusters", n_clusters))
    
    # Step 4: Add UMAP embedding
    log_msg("step", "Adding UMAP embedding...")
    proj <- addUMAP(
      ArchRProj = proj,
      reducedDims = "IterativeLSI",
      name = "UMAP",
      nNeighbors = 30,
      minDist = 0.5,
      metric = "cosine",
      force = TRUE
    )
    log_msg("step", "UMAP added")
    
    # Step 5: Save ArchR project with all artifacts (LSI, Clusters, UMAP)
    # The project already has outputDirectory set, so just finalize the save
    log_msg("step", "Finalizing ArchR project save with all computed artifacts...")
    saveArchRProject(
      ArchRProj = proj,
      load = FALSE,
      overwrite = TRUE
    )
    log_msg("step", sprintf("ArchR project finalized with: LSI, Clusters (%d unique), UMAP, %d cells", 
                           length(unique(proj$Clusters)), ncol(proj)))
    log_msg("step", sprintf("Project path: %s", proj@projectMetadata$outputDirectory))
    
    # # Also save RDS backup for quick loading
    # rds_backup <- file.path(obj_output_dir, sprintf("%s_archR_final.rds", obj_name))
    # saveRDS(proj, rds_backup)
    # log_msg("step", sprintf("RDS backup saved: %s", rds_backup))
    
    # Step 6: Save final cell barcodes
    barcode_output <- file.path(obj_output_dir, sprintf("%s_final_cell_barcodes.txt", obj_name))
    final_barcodes <- proj$cellNames
    writeLines(final_barcodes, barcode_output)
    log_msg("step", sprintf("Saved barcode list (%d cells): %s", length(final_barcodes), barcode_output))
    
    # Step 7: Generate QC plots (4 panels: TSS, nFrags, scatter, UMAP)
    log_msg("step", "Generating QC plots...")
    
    final_metadata <- getCellColData(proj)
    final_tss_vec <- final_metadata$TSSEnrichment
    final_nfrags_vec <- final_metadata$nFrags
    
    # Panel 1: TSS Distribution
    p1 <- ggplot(data.frame(TSS = final_tss_vec), aes(x = TSS)) +
      geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
      geom_vline(xintercept = min_tss, linetype = "dashed", color = "red", size = 1) +
      labs(title = sprintf("%s: TSS Enrichment (n=%d)", obj_name, ncol(proj)),
           x = "TSS Enrichment Score",
           y = "Number of Cells") +
      theme_minimal() +
      theme(plot.title = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 12))
    
    # Panel 2: Fragment Count Distribution
    p2 <- ggplot(data.frame(nFrags = final_nfrags_vec), aes(x = nFrags)) +
      geom_histogram(bins = 50, fill = "darkgreen", alpha = 0.7) +
      geom_vline(xintercept = min_frags, linetype = "dashed", color = "red", size = 1) +
      labs(title = sprintf("%s: Fragment Count (n=%d)", obj_name, ncol(proj)),
           x = "nFragments",
           y = "Number of Cells") +
      scale_x_log10() +
      theme_minimal() +
      theme(plot.title = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 12))
    
    # Panel 3: TSS vs Fragments scatter
    p3 <- ggplot(data.frame(TSS = final_tss_vec, nFrags = final_nfrags_vec), 
                 aes(x = nFrags, y = TSS)) +
      geom_point(alpha = 0.5, size = 1, color = "steelblue") +
      geom_hline(yintercept = min_tss, linetype = "dashed", color = "red") +
      geom_vline(xintercept = min_frags, linetype = "dashed", color = "red") +
      labs(title = sprintf("%s: TSS vs Fragment Count (n=%d)", obj_name, ncol(proj)),
           x = "nFragments",
           y = "TSS Enrichment") +
      scale_x_log10() +
      theme_minimal() +
      theme(plot.title = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 12))
    
    # Panel 4: UMAP Clustering
    p4 <- tryCatch({
      if (!is.null(proj@embeddings$UMAP)) {
        umap_mat <- proj@embeddings$UMAP$df
        cluster_vec <- proj$Clusters
        
        umap_df <- data.frame(
          UMAP1 = umap_mat[, 1],
          UMAP2 = umap_mat[, 2],
          Cluster = cluster_vec
        )
        
        ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
          geom_point(alpha = 0.6, size = 2) +
          labs(title = sprintf("%s: UMAP Clustering (n=%d, k=%d clusters)", obj_name, ncol(proj), length(unique(proj$Clusters))),
               x = "UMAP1",
               y = "UMAP2") +
          theme_minimal() +
          theme(plot.title = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 12),
                legend.position = "right")
      } else {
        stop("No UMAP embedding")
      }
    }, error = function(e) {
      log_msg("warn", sprintf("UMAP plot failed: %s", e$message))
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                label = "UMAP plot\nnot available", 
                size = 3, ha = "center") +
        theme_void()
    })
    
    # Save PDF with 4-panel layout
    pdf(file.path(plot_output, sprintf("archR_qc_%s.pdf", obj_name)), width = 14, height = 10)
    grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
    dev.off()
    
    log_msg("step", sprintf("Saved plot: archR_qc_%s.pdf", obj_name))
    
  }, error = function(e) {
    log_msg("error", sprintf("Failed to process %s: %s", obj_name, e$message))
  })
}

log_msg("success", "===== ArchR Project Creation Complete (NO DOUBLET REMOVAL) =====")
