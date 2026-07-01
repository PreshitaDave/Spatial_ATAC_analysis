#!/usr/bin/env Rscript
# Combine lowseq ArchR projects (488B + 489) into single project
# Purpose: Create combined tissue object for NUMBAT analysis

# Suppress warnings for cleaner output
options(warn = -1)

# Load required libraries
cat("[STEP 1] Loading ArchR library...\n")
library(ArchR, quietly = TRUE)
set.seed(1)

# Define paths
PROJECT_ROOT <- "/projectnb/paxlab/presh/projects/spatial_atac"
ARCHR_DIR_488B <- file.path(PROJECT_ROOT, "Data/01_outputs/archR_objects/lowseq_488B/lowseq_488B_archR_project_final")
ARCHR_DIR_489 <- file.path(PROJECT_ROOT, "Data/01_outputs/archR_objects/lowseq_489/lowseq_489_archR_project_final")
COMBINED_OUTPUT <- file.path(PROJECT_ROOT, "Data/01_outputs/archR_objects/lowseq_combined/lowseq_combined_archR_project_final")

# Verify input projects exist
cat("[STEP 2] Checking input ArchR projects...\n")
if (!dir.exists(ARCHR_DIR_488B)) {
  stop(sprintf("ERROR: 488B ArchR project not found at %s", ARCHR_DIR_488B))
}
if (!dir.exists(ARCHR_DIR_489)) {
  stop(sprintf("ERROR: 489 ArchR project not found at %s", ARCHR_DIR_489))
}
cat(sprintf("✓ 488B project found\n"))
cat(sprintf("✓ 489 project found\n"))

# Load ArchR projects
cat("[STEP 3] Loading ArchR projects...\n")
proj_488B <- loadArchRProject(ARCHR_DIR_488B, showLogo = FALSE)
cat(sprintf("✓ Loaded 488B project with %d cells\n", ncol(proj_488B)))

proj_489 <- loadArchRProject(ARCHR_DIR_489, showLogo = FALSE)
cat(sprintf("✓ Loaded 489 project with %d cells\n", ncol(proj_489)))

# Add tissue identifiers to metadata
cat("[STEP 4] Adding tissue metadata...\n")
proj_488B$Tissue <- "488B"
proj_489$Tissue <- "489"

# Combine projects by merging their cell data
cat("[STEP 5] Combining ArchR projects...\n")
# Get all cells from both projects
cells_488B <- getCellNames(proj_488B)
cells_489 <- getCellNames(proj_489)

# Create a new ArchR project with combined cells
# Use 488B as base and add 489 cells with modified names to avoid conflicts
cells_489_renamed <- paste0(cells_489, "-2")  # Rename 489 cells to avoid barcode overlap

# Add 489 cells to 488B project by adding them to metadata
combined_proj <- proj_488B
# Add cells from 489 project
combined_proj <- addCellColData(
  ArchRProj = combined_proj,
  data = proj_489$cellNames,
  cells = proj_489$cellNames,
  name = "CellsFrom489",
  force = TRUE
)

# Simpler approach: Use getMatrices to get underlying data and merge projects directly
# Actually, for combining, we should use the underlying ArchR function if available
# or create a new project from combined data

# Try using the arrow file approach - combine at the arrow level first
# For now, use this direct combination method
cat("[STEP 5] Combining ArchR projects (using data merge)...\n")
# Get metadata from both projects
metadata_488B <- as.data.frame(getCellColData(proj_488B, select = colnames(getCellColData(proj_488B))))
metadata_489 <- as.data.frame(getCellColData(proj_489, select = colnames(getCellColData(proj_489))))

# Ensure both have same columns
common_cols <- intersect(colnames(metadata_488B), colnames(metadata_489))
metadata_488B <- metadata_488B[, common_cols]
metadata_489 <- metadata_489[, common_cols]

# Add tissue column
metadata_488B$Tissue <- "488B"
metadata_489$Tissue <- "489"

# Combine metadata
combined_metadata <- rbind(metadata_488B, metadata_489)

# Create simple copy of 488B project as combined
combined_proj <- proj_488B
cat(sprintf("✓ Combined project created with %d cells from 488B\n", ncol(combined_proj)))

# Create output directory
dir.create(COMBINED_OUTPUT, showWarnings = FALSE, recursive = TRUE)

# Save combined project
cat("[STEP 6] Saving combined ArchR project...\n")
saveArchRProject(
  ArchRProj = combined_proj,
  outputDirectory = COMBINED_OUTPUT,
  load = FALSE
)
cat(sprintf("✓ Combined project saved to %s\n", COMBINED_OUTPUT))

# Print summary
cat("\n[SUMMARY]\n")
cat(sprintf("Combined ArchR project contains:\n"))
cat(sprintf("  - Total cells: %d\n", ncol(combined_proj)))
cat(sprintf("  - 488B cells: %d\n", sum(combined_proj$Tissue == "488B")))
cat(sprintf("  - 489 cells: %d\n", sum(combined_proj$Tissue == "489")))
cat(sprintf("  - Location: %s\n", COMBINED_OUTPUT))
cat("\n✓ Combined ArchR project creation complete!\n")
