#!/usr/bin/env Rscript
# ============================================================================
# 0_create_archr_qc_cluster.R
# Create ArchR objects from filtered fragments, QC filter, remove doublets,
# cluster, add UMAP, and save comprehensive PDFs for all objects
# ============================================================================
# Usage:
#   Rscript 0_create_archr_qc_cluster.R [--project-root /path] [--threads N]
# Environment:
#   Set PROJECT_ROOT, NSLOTS for cluster resource control
# Output:
#   - Data/01_outputs/archR_objects/{object}/{object}_archR_project/
#   - Data/01_inputs/arrow/{object}.arrow (Arrow files)
#   - analysis/plots/archr_obj/archR_qc_{object}.pdf (comprehensive QC plots)
# ============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(data.table)
  library(dplyr)
})

# ============================================================================
# CONFIGURATION
# ============================================================================

project_root <- Sys.getenv("PROJECT_ROOT", "/projectnb/paxlab/presh/projects/spatial_atac")
threads <- as.integer(Sys.getenv("NSLOTS", "4"))
min_tss <- 3
min_frags <- 1000
max_frags <- Inf
doublet_cutoff <- 2

set.seed(42)
addArchRGenome("hg38")
addArchRThreads(threads = threads)

# ============================================================================
# SETUP PATHS
# ============================================================================

log_msg <- function(tag, msg) {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), tag, msg))
}

data_dir <- file.path(project_root, "Data", "01_inputs")
output_root <- file.path(project_root, "Data", "01_outputs")
archr_output <- file.path(output_root, "archR_objects")
arrow_output <- file.path(data_dir, "arrow")
plot_output <- file.path(project_root, "analysis", "plots", "archr_obj")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(archr_output, recursive = TRUE, showWarnings = FALSE)
dir.create(arrow_output, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_output, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# DEFINE OBJECTS & PATHS
# ============================================================================

objects <- list(
  deepseq_488B = list(
    fragments = file.path(data_dir, "fragments", "deepseq_488B", "deepseq_488B.fragments.sort.filtered.bed.gz"),
    barcodes = file.path(data_dir, "barcodes", "tissue_barcodes", "deepseq_488B", "deepseq_488B.no_edge_effect.barcodes.tsv"),
    name = "deepseq_488B",
    sample_name = "Deepseq_488B"
  ),
  deepseq_489 = list(
    fragments = file.path(data_dir, "fragments", "deepseq_489", "deepseq_489.fragments.sort.filtered.bed.gz"),
    barcodes = file.path(data_dir, "barcodes", "tissue_barcodes", "deepseq_489", "deepseq_489.no_edge_effect.barcodes.tsv"),
    name = "deepseq_489",
    sample_name = "Deepseq_489"
  ),
  lowseq_488B = list(
    fragments = file.path(data_dir, "fragments", "lowseq_488B", "lowseq_488B.fragments.sort.filtered.bed.gz"),
    barcodes = file.path(data_dir, "barcodes", "tissue_barcodes", "lowseq_488B", "lowseq_488B.no_edge_effect.barcodes.tsv"),
    name = "lowseq_488B",
    sample_name = "Lowseq_488B"
  ),
  lowseq_489 = list(
    fragments = file.path(data_dir, "fragments", "lowseq_489", "lowseq_489.fragments.sort.filtered.bed.gz"),
    barcodes = file.path(data_dir, "barcodes", "tissue_barcodes", "lowseq_489", "lowseq_489.no_edge_effect.barcodes.tsv"),
    name = "lowseq_489",
    sample_name = "Lowseq_489"
  )
)

log_msg("config", sprintf("Threads: %d | Min TSS: %g | Min Frags: %d", threads, min_tss, min_frags))
log_msg("config", sprintf("Output: %s", archr_output))

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

validate_input_files <- function(obj_info) {
  issues <- c()
  if (!file.exists(obj_info$fragments)) {
    issues <- c(issues, sprintf("Missing fragments: %s", obj_info$fragments))
  }
  if (!file.exists(obj_info$barcodes)) {
    issues <- c(issues, sprintf("Missing barcodes: %s", obj_info$barcodes))
  }
  if (length(issues) > 0) {
    stop(paste("Input validation failed for", obj_info$name, ":\n", paste(issues, collapse = "\n")))
  }
}

read_barcodes <- function(path) {
  if (!file.exists(path)) return(character(0))
  bcs <- readLines(path)
  # Remove -1 suffix if present (ArchR convention)
  bcs <- sub("-1$", "", bcs)
  # Clean empty lines
  bcs <- bcs[nzchar(bcs)]
  unique(bcs)
}

normalize_barcode <- function(barcode) {
  # Remove -1 suffix and any prefixes (e.g., Deepseq#)
  bc <- as.character(barcode)
  bc <- sub("-1$", "", bc)
  bc <- sub("^.*#", "", bc)
  bc
}

create_arrow_from_fragments <- function(fragments_path, sample_name, output_dir, threads = 4) {
  # Create Arrow file from fragment BED file
  log_msg("arrow", sprintf("Creating Arrow from: %s", fragments_path))
  
  arrow_output <- file.path(output_dir, paste0(sample_name, ".arrow"))
  
  if (file.exists(arrow_output)) {
    log_msg("arrow", sprintf("Arrow file exists, skipping creation: %s", arrow_output))
    return(arrow_output)
  }
  
  tryCatch({
    # Change to output directory for arrow creation
    original_wd <- getwd()
    setwd(output_dir)
    
    ArrowFiles <- createArrowFiles(
      inputFiles = fragments_path,
      sampleNames = sample_name,
      outputNames = sample_name,  # Don't add .arrow - ArchR adds it automatically
      minTSS = 2,
      minFrags = 100,
      maxFrags = Inf,
      addTileMat = TRUE,
      addGeneScoreMat = TRUE,
      force = TRUE,
      threads = threads,
      verbose = TRUE
    )
    
    # Restore original working directory
    setwd(original_wd)
    
    log_msg("arrow", sprintf("Created Arrow: %s", arrow_output))
    return(arrow_output)
  }, error = function(e) {
    setwd(original_wd)
    log_msg("error", sprintf("Failed to create Arrow: %s", e$message))
    return(NULL)
  })
}

filter_to_barcodes <- function(proj, barcode_list) {
  # Subset ArchR project to specific barcodes
  all_cells <- getCellNames(proj)
  all_cells_norm <- normalize_barcode(all_cells)
  
  matched_idx <- which(all_cells_norm %in% barcode_list)
  matched_cells <- all_cells[matched_idx]
  
  log_msg("filter", sprintf("Matched %d/%d cells to barcode list", length(matched_cells), length(barcode_list)))
  
  if (length(matched_cells) == 0) {
    stop("No matching cells found!")
  }
  
  proj[matched_cells, ]
}

plot_qc_metrics <- function(proj, sample_name) {
  # Generate QC metric plots
  plots <- list()
  
  # TSS vs nFrags
  df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))
  df <- as.data.frame(df)
  colnames(df) <- c("nFrags", "TSS")
  
  plots$tss_nfrags <- ggplot(df, aes(x = nFrags, y = TSS)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_hline(yintercept = 4, linetype = "dashed", color = "red", alpha = 0.5) +
    geom_vline(xintercept = 3, linetype = "dashed", color = "red", alpha = 0.5) +
    labs(title = sprintf("%s: TSS vs nFrags", sample_name),
         x = "Log10(Fragments)", y = "TSS Enrichment") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
  
  # nFrags histogram
  plots$nfrags_hist <- ggplot(df, aes(x = nFrags)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = 3, linetype = "dashed", color = "red", alpha = 0.5) +
    labs(title = sprintf("%s: Fragment Count Distribution", sample_name),
         x = "Log10(Fragments)", y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
  
  # TSS histogram
  plots$tss_hist <- ggplot(df, aes(x = TSS)) +
    geom_histogram(bins = 50, fill = "steelgreen", alpha = 0.7) +
    geom_vline(xintercept = 4, linetype = "dashed", color = "red", alpha = 0.5) +
    labs(title = sprintf("%s: TSS Enrichment Distribution", sample_name),
         x = "TSS Enrichment", y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
  
  # Cell count by pass/fail
  tss_pass <- sum(df$TSS >= 4)
  frags_pass <- sum(df$nFrags >= 3)
  both_pass <- sum(df$TSS >= 4 & df$nFrags >= 3)
  
  qc_data <- data.frame(
    Metric = c("TSS Pass", "Frags Pass", "Both Pass", "Total"),
    Count = c(tss_pass, frags_pass, both_pass, nrow(df))
  )
  
  plots$qc_summary <- ggplot(qc_data, aes(x = Metric, y = Count, fill = Metric)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    geom_text(aes(label = Count), vjust = -0.5, size = 3) +
    labs(title = sprintf("%s: QC Summary", sample_name),
         y = "Cell Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  plots
}

plot_clustering_umap <- function(proj, sample_name) {
  # Plot clustering and UMAP
  plots <- list()
  
  # UMAP by cluster
  if (!is.null(proj@embeddings$UMAP)) {
    plots$umap_clusters <- plotEmbedding(
      ArchRProj = proj,
      colorBy = "cellColData",
      name = "Clusters",
      embedding = "UMAP",
      size = 1.5,
      labelAsFactors = FALSE
    ) + ggtitle(sprintf("%s: UMAP - Clusters", sample_name)) +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
  }
  
  # UMAP by TSS
  if (!is.null(proj@embeddings$UMAP)) {
    plots$umap_tss <- plotEmbedding(
      ArchRProj = proj,
      colorBy = "cellColData",
      name = "TSSEnrichment",
      embedding = "UMAP",
      size = 1.5
    ) + ggtitle(sprintf("%s: UMAP - TSS Enrichment", sample_name)) +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
  }
  
  # UMAP by nFrags
  if (!is.null(proj@embeddings$UMAP)) {
    plots$umap_frags <- plotEmbedding(
      ArchRProj = proj,
      colorBy = "cellColData",
      name = "log10(nFrags)",
      embedding = "UMAP",
      size = 1.5
    ) + ggtitle(sprintf("%s: UMAP - nFrags", sample_name)) +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
  }
  
  plots
}

# ============================================================================
# MAIN PROCESSING LOOP
# ============================================================================

log_msg("start", "===== ArchR Object Creation & QC Pipeline =====")
log_msg("start", sprintf("Processing %d objects", length(objects)))

all_plots <- list()
summary_data <- list()

for (obj_name in names(objects)) {
  obj_info <- objects[[obj_name]]
  log_msg("process", sprintf("--- Processing %s ---", obj_name))
  
  # Validate inputs
  validate_input_files(obj_info)
  
  # Create object output directory
  obj_output_dir <- file.path(archr_output, obj_name)
  dir.create(obj_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  tryCatch({
    # 1. Load or create Arrow file
    # Arrow files are saved to Data/01_inputs/arrow/
    arrow_file <- create_arrow_from_fragments(
      obj_info$fragments,
      obj_info$sample_name,
      arrow_output,
      threads = threads
    )
    
    if (is.null(arrow_file)) {
      log_msg("error", sprintf("Failed to create Arrow for %s", obj_name))
      next
    }
    
    # 2. Create ArchR project
    log_msg("step", sprintf("Creating ArchR project from Arrow: %s", arrow_file))
    proj <- ArchRProject(
      ArrowFiles = arrow_file,
      outputDirectory = file.path(obj_output_dir, paste0(obj_name, "_archR_project")),
      copyArrows = FALSE
    )
    log_msg("step", sprintf("Created project with %d cells", ncol(proj)))
    
    # 3. Filter to no_edge_effect barcodes
    barcodes <- read_barcodes(obj_info$barcodes)
    log_msg("step", sprintf("Read %d no_edge_effect barcodes", length(barcodes)))
    
    proj_filt <- filter_to_barcodes(proj, barcodes)
    log_msg("step", sprintf("Filtered to %d cells", ncol(proj_filt)))
    
    # 4. Apply QC filters
    log_msg("step", "Applying QC filters")
    before_qc <- ncol(proj_filt)
    
    # Remove cells failing QC thresholds
    # Note: nFrags column stores log10-transformed values
    tss_vec <- proj_filt$TSSEnrichment
    nfrags_vec <- proj_filt$nFrags
    qc_pass <- tss_vec >= min_tss & nfrags_vec >= log10(min_frags)
    proj_qc <- proj_filt[qc_pass, ]
    
    log_msg("step", sprintf("After QC: %d/%d cells pass (%.1f%% retained)",
                            ncol(proj_qc), before_qc, 100 * ncol(proj_qc) / before_qc))
    
    # 5. Add iterative LSI for dimensionality reduction
    log_msg("step", "Computing iterative LSI")
    proj_qc <- addIterativeLSI(
      ArchRProj = proj_qc,
      useMatrix = "TileMatrix",
      name = "IterativeLSI",
      iterations = 2,
      clusterParams = list(
        resolution = c(0.2),
        sampleCells = 10000,
        n.start = 10
      ),
      varFeatures = 25000,
      dimsToUse = 1:30,
      force = TRUE
    )
    log_msg("step", "LSI computation complete")
    
    # 6. Cluster cells
    log_msg("step", "Clustering cells")
    proj_qc <- addClusters(
      input = proj_qc,
      reducedDims = "IterativeLSI",
      method = "Seurat",
      name = "Clusters",
      resolution = 0.8,
      force = TRUE
    )
    log_msg("step", sprintf("Found %d clusters", length(unique(proj_qc$Clusters))))
    
    # 7. Compute UMAP
    log_msg("step", "Computing UMAP")
    proj_qc <- addUMAP(
      ArchRProj = proj_qc,
      reducedDims = "IterativeLSI",
      name = "UMAP",
      force = TRUE
    )
    log_msg("step", "UMAP complete")
    
    # 8. Detect and remove doublets
    log_msg("step", "Detecting doublets using fragment-based approach")
    before_doublet <- ncol(proj_qc)
    
    # Method 1: Check if nDoublets column exists (from Arrow file)
    has_doublets_col <- "nDoublets" %in% colnames(getCellColData(proj_qc))
    
    if (has_doublets_col) {
      log_msg("step", "Using nDoublets column from Arrow file")
      doublet_pass <- proj_qc$nDoublets < doublet_cutoff
      doublets_detected <- sum(!doublet_pass)
    } else {
      log_msg("step", "Computing doublet score based on fragment distribution")
      # Method 2: Flag potential doublets as high-fragment outliers
      # Note: nFrags column stores log10-transformed values
      nfrags <- proj_qc$nFrags
      # Identify high-fragment outliers (potential doublets)
      # Use 1.5*IQR + 3rd quartile as cutoff
      q3 <- quantile(nfrags, 0.75, na.rm = TRUE)
      iqr <- IQR(nfrags, na.rm = TRUE)
      upper_cutoff <- q3 + 1.5 * iqr
      
      doublet_pass <- nfrags <= upper_cutoff
      doublets_detected <- sum(!doublet_pass)
      log_msg("step", sprintf("Fragment cutoff for doublets: %.2f (Q3 + 1.5*IQR)", upper_cutoff))
    }
    
    if (doublets_detected > 0) {
      proj_qc <- proj_qc[doublet_pass, ]
      log_msg("step", sprintf("After doublet removal: %d/%d cells removed, %d cells retained (%.1f%%)",
                              doublets_detected, before_doublet, ncol(proj_qc), 
                              100 * ncol(proj_qc) / before_doublet))
    } else {
      log_msg("step", sprintf("No doublets detected - all %d cells retained", ncol(proj_qc)))
    }
    
    # 9. Generate comprehensive plots
    log_msg("step", "Generating plots")
    obj_plots <- list()
    
    # QC plots
    qc_plots <- plot_qc_metrics(proj_qc, obj_name)
    obj_plots <- c(obj_plots, qc_plots)
    
    # Clustering plots
    cluster_plots <- plot_clustering_umap(proj_qc, obj_name)
    obj_plots <- c(obj_plots, cluster_plots)
    
    all_plots[[obj_name]] <- obj_plots
    
    # 10. Save ArchR project
    log_msg("step", sprintf("Saving ArchR project"))
    saveArchRProject(
      ArchRProj = proj_qc,
      outputDirectory = file.path(obj_output_dir, paste0(obj_name, "_archR_project_final")),
      load = FALSE,
      overwrite = TRUE
    )
    log_msg("step", sprintf("ArchR project saved to: %s", file.path(obj_output_dir, paste0(obj_name, "_archR_project_final"))))
    
    # Save as RDS for quick loading
    rds_path <- file.path(obj_output_dir, paste0(obj_name, "_archR_final.rds"))
    saveRDS(proj_qc, rds_path)
    log_msg("step", sprintf("Saved RDS: %s", rds_path))
    
    # Record summary
    summary_data[[obj_name]] <- data.frame(
      Object = obj_name,
      Initial_Cells = before_qc,
      After_QC = ncol(proj_qc),
      QC_Retention_Pct = round(100 * ncol(proj_qc) / before_qc, 1),
      Clusters = length(unique(proj_qc$Clusters)),
      Status = "SUCCESS"
    )
    
    log_msg("complete", sprintf("Successfully processed %s", obj_name))
    
  }, error = function(e) {
    log_msg("error", sprintf("Failed to process %s: %s", obj_name, e$message))
    summary_data[[obj_name]] <<- data.frame(
      Object = obj_name,
      Status = "FAILED",
      Error = e$message
    )
  })
}

# ============================================================================
# SAVE COMPREHENSIVE PDF PLOTS
# ============================================================================

log_msg("step", "Saving comprehensive PDF reports")

for (obj_name in names(all_plots)) {
  pdf_path <- file.path(plot_output, sprintf("archR_qc_%s.pdf", obj_name))
  
  tryCatch({
    pdf(pdf_path, width = 14, height = 10, onefile = TRUE)
    
    # Title page
    plot.new()
    grid.text(sprintf("ArchR QC Report: %s", obj_name),
              gp = gpar(fontsize = 28, fontface = "bold"),
              y = 0.9)
    grid.text(sprintf("Generated: %s", format(Sys.time(), "%F %T")),
              gp = gpar(fontsize = 12),
              y = 0.8)
    
    # Plot all available plots
    plots <- all_plots[[obj_name]]
    for (i in seq_along(plots)) {
      print(plots[[i]])
    }
    
    dev.off()
    log_msg("step", sprintf("Saved PDF: %s", pdf_path))
  }, error = function(e) {
    log_msg("error", sprintf("Failed to save PDF for %s: %s", obj_name, e$message))
  })
}

# ============================================================================
# SUMMARY REPORT
# ============================================================================

log_msg("summary", "===== Processing Summary =====")

if (length(summary_data) > 0) {
  summary_df <- do.call(rbind, summary_data)
  summary_df <- as.data.frame(summary_df)
  
  cat("\n")
  print(summary_df)
  cat("\n")
  
  # Save summary to file
  summary_path <- file.path(plot_output, "archR_processing_summary.tsv")
  write.table(summary_df, summary_path, sep = "\t", row.names = FALSE, quote = FALSE)
  log_msg("summary", sprintf("Saved summary: %s", summary_path))
}

log_msg("complete", "===== ArchR Pipeline Complete =====")
log_msg("complete", sprintf("ArchR objects: %s", archr_output))
log_msg("complete", sprintf("PDF reports: %s", plot_output))

quit("no", 0)
