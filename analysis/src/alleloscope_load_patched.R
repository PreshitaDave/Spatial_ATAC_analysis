#' Load patched Alleloscope functions with character coercions
#' 
#' This script sources Alleloscope R functions directly, bypassing the cached
#' installed package. All functions include defensive character coercions before
#' strsplit() calls to prevent "non-character argument" errors.
#'
#' Usage: source("analysis/src/alleloscope_load_patched.R")

allelo_base <- "/projectnb/paxlab/presh/software/Alleloscope"

# Load in dependency order
# Core utilities first
source(file.path(allelo_base, "R/genotype_ref.R"), local=FALSE)
source(file.path(allelo_base, "R/genotype_conf.R"), local=FALSE)

# Matrix/filtering operations
source(file.path(allelo_base, "R/Matrix_filter.R"), local=FALSE)
source(file.path(allelo_base, "R/Segments_filter.R"), local=FALSE)

# Core analysis functions (PATCHED with character coercions)
source(file.path(allelo_base, "R/Genotype_value.R"), local=FALSE)
source(file.path(allelo_base, "R/Genotype.R"), local=FALSE)
source(file.path(allelo_base, "R/Cov_value.R"), local=FALSE)

# Segmentation and EM
source(file.path(allelo_base, "R/Segmentation.R"), local=FALSE)
source(file.path(allelo_base, "R/EM.R"), local=FALSE)

# Plotting and downstream
source(file.path(allelo_base, "R/plot_scATAC_cnv.R"), local=FALSE)
source(file.path(allelo_base, "R/Lineage_plot.R"), local=FALSE)

# Supporting functions
source(file.path(allelo_base, "R/Est_regions.R"), local=FALSE)
source(file.path(allelo_base, "R/Createobj.R"), local=FALSE)
source(file.path(allelo_base, "R/Select_normal.R"), local=FALSE)
source(file.path(allelo_base, "R/AssignClones_ref.R"), local=FALSE)
source(file.path(allelo_base, "R/genotype_neighbor.R"), local=FALSE)

message("[OK] Patched Alleloscope functions loaded with character coercions")
