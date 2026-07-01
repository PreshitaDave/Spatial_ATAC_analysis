#!/bin/bash
set -eo pipefail

# ============================================================================
# DEBUG: Run ArchR QC pipeline for deepseq_488B AND deepseq_489 with doublet removal
# ============================================================================

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
cd "$PROJECT_ROOT"

# Initialize module system
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

echo "[DEBUG STEP 1] Verifying compute node"
HOSTNAME=$(hostname)
echo "[DEBUG] Hostname: $HOSTNAME"
if [[ "$HOSTNAME" == scc1* ]]; then
  echo "[ERROR] Running on login node! Must use compute node."
  exit 1
fi

echo "[DEBUG STEP 2] Loading R module"
module load R
which Rscript
Rscript --version

echo "[DEBUG STEP 3] Verifying input files"
for TISSUE in deepseq_488B deepseq_489; do
  echo "[DEBUG] Checking $TISSUE..."
  ls -lh Data/01_inputs/fragments/$TISSUE/${TISSUE}.fragments.sort.filtered.bed.gz
  ls -lh Data/01_inputs/barcodes/tissue_barcodes/$TISSUE/${TISSUE}.no_edge_effect.barcodes.tsv
done

echo "[DEBUG STEP 4] Creating output directories"
mkdir -p Data/01_outputs/archR_objects
mkdir -p Data/01_inputs/arrow
mkdir -p analysis/plots/archr_obj

# Create the R script that tests both tissues
cat > /tmp/archr_debug_deepseq.R << 'REOF'
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(data.table)
  library(dplyr)
})

project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
threads <- 4
min_tss <- 4
min_frags <- 1000
doublet_cutoff <- 2

set.seed(42)
addArchRGenome("hg38")
addArchRThreads(threads = threads)

log_msg <- function(tag, msg) {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), tag, msg))
}

data_dir <- file.path(project_root, "Data", "01_inputs")
output_root <- file.path(project_root, "Data", "01_outputs")
archr_output <- file.path(output_root, "archR_objects")
arrow_output <- file.path(data_dir, "arrow")

# Define tissues to process
tissues <- list(
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
  )
)

log_msg("debug", "=== Starting debug run for deepseq samples ===")

# Normalize barcode function
normalize_barcode <- function(bc) {
  bc <- as.character(bc)
  bc <- sub("-1$", "", bc)
  bc <- sub("^.*#", "", bc)
  bc
}

# Process each tissue
for (tissue_name in names(tissues)) {
  log_msg("process", sprintf("--- Processing %s ---", tissue_name))
  
  obj_info <- tissues[[tissue_name]]
  
  # Check files exist
  if (!file.exists(obj_info$fragments)) {
    log_msg("error", sprintf("Missing fragments: %s", obj_info$fragments))
    next
  }
  if (!file.exists(obj_info$barcodes)) {
    log_msg("error", sprintf("Missing barcodes: %s", obj_info$barcodes))
    next
  }
  
  # Load Arrow if exists, else create
  arrow_file <- file.path(arrow_output, paste0(obj_info$sample_name, ".arrow"))
  if (!file.exists(arrow_file)) {
    log_msg("step", "Creating Arrow from fragments...")
    original_wd <- getwd()
    setwd(arrow_output)
    
    ArrowFiles <- createArrowFiles(
      inputFiles = obj_info$fragments,
      sampleNames = obj_info$sample_name,
      outputNames = obj_info$sample_name,
      minTSS = 2,
      minFrags = 100,
      addTileMat = TRUE,
      addGeneScoreMat = TRUE,
      force = TRUE,
      threads = threads,
      verbose = FALSE  # Reduce verbosity
    )
    
    setwd(original_wd)
    log_msg("step", sprintf("Arrow created: %s", arrow_file))
  } else {
    log_msg("step", sprintf("Arrow exists, using existing: %s", arrow_file))
  }
  
  # Create ArchR project
  log_msg("step", "Creating ArchR project...")
  obj_output_dir <- file.path(archr_output, tissue_name)
  dir.create(obj_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  proj <- ArchRProject(
    ArrowFiles = arrow_file,
    outputDirectory = file.path(obj_output_dir, paste0(tissue_name, "_archR_project")),
    copyArrows = FALSE
  )
  log_msg("step", sprintf("Project created with %d cells", ncol(proj)))
  
  # Filter to no_edge_effect barcodes
  log_msg("step", "Reading barcodes...")
  barcodes <- readLines(obj_info$barcodes)
  barcodes <- sub("-1$", "", barcodes)
  barcodes <- barcodes[nzchar(barcodes)]
  barcodes <- unique(barcodes)
  log_msg("step", sprintf("Read %d barcodes", length(barcodes)))
  
  # Filter cells
  all_cells <- getCellNames(proj)
  all_cells_norm <- normalize_barcode(all_cells)
  matched_idx <- which(all_cells_norm %in% barcodes)
  matched_cells <- all_cells[matched_idx]
  
  log_msg("step", sprintf("Matched %d/%d cells to barcode list", length(matched_cells), length(barcodes)))
  proj_filt <- proj[matched_cells, ]
  
  # Apply QC filters - WITH FIX
  log_msg("step", "Applying QC filters...")
  before_qc <- ncol(proj_filt)
  
  tss_vec <- proj_filt$TSSEnrichment
  nfrags_vec <- proj_filt$nFrags
  qc_pass <- tss_vec >= min_tss & nfrags_vec >= log10(min_frags)
  proj_qc <- proj_filt[qc_pass, ]
  
  log_msg("step", sprintf("After QC: %d/%d cells pass (%.1f%% retained)", 
                          ncol(proj_qc), before_qc, 100 * ncol(proj_qc) / before_qc))
  
  # Add iterative LSI
  log_msg("step", "Adding iterative LSI...")
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
  log_msg("step", "LSI complete")
  
  # Clustering
  log_msg("step", "Clustering...")
  proj_qc <- addClusters(
    input = proj_qc,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8,
    force = TRUE
  )
  log_msg("step", sprintf("Found %d clusters", length(unique(proj_qc$Clusters))))
  
  # UMAP
  log_msg("step", "Computing UMAP...")
  proj_qc <- addUMAP(
    ArchRProj = proj_qc,
    reducedDims = "IterativeLSI",
    name = "UMAP",
    force = TRUE
  )
  log_msg("step", "UMAP complete")
  
  # Doublet detection and removal - IMPROVED
  log_msg("step", "Detecting doublets...")
  before_doublet <- ncol(proj_qc)
  
  # Check if nDoublets column exists
  has_doublets_col <- "nDoublets" %in% colnames(getCellColData(proj_qc))
  
  if (has_doublets_col) {
    log_msg("step", "Using nDoublets column from Arrow file")
    doublet_pass <- proj_qc$nDoublets < doublet_cutoff
    doublets_detected <- sum(!doublet_pass)
  } else {
    log_msg("step", "Computing doublet score based on fragment distribution")
    # Flag high-fragment outliers as potential doublets
    # Note: nFrags stores log10-transformed values
    nfrags <- proj_qc$nFrags
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
  
  # Save project
  log_msg("step", "Saving project...")
  saveArchRProject(
    ArchRProj = proj_qc,
    outputDirectory = file.path(obj_output_dir, paste0(tissue_name, "_archR_project_final")),
    load = FALSE,
    overwrite = TRUE
  )
  
  rds_path <- file.path(obj_output_dir, paste0(tissue_name, "_archR_final.rds"))
  saveRDS(proj_qc, rds_path)
  log_msg("step", sprintf("RDS saved: %s", rds_path))
  
  log_msg("success", sprintf("Completed %s successfully", tissue_name))
}

log_msg("debug", "=== All samples processed successfully ===")
REOF

echo "[DEBUG STEP 5] Running ArchR pipeline for both deepseq tissues"
export PROJECT_ROOT
Rscript /tmp/archr_debug_deepseq.R
EXIT_CODE=$?

if [[ $EXIT_CODE -eq 0 ]]; then
  echo "[DEBUG SUCCESS] Both deepseq samples processed successfully"
  echo "[DEBUG] Ready to submit parallel jobs for remaining 2 tissues (lowseq_488B, lowseq_489)"
else
  echo "[DEBUG FAILED] Processing failed with exit code $EXIT_CODE"
  exit $EXIT_CODE
fi
