#!/usr/bin/env Rscript
# Test run_numbat() function directly
args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
atac_file <- args[2]
allele_file <- args[3]
ncores <- as.numeric(args[4])

ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

message(sprintf("[%s] Loading libraries...", ts()))
library(numbat)
library(data.table)

message(sprintf("[%s] Loading data...", ts()))
count_mat <- readRDS(atac_file)
message(sprintf("[%s] ATAC: %d × %d", ts(), nrow(count_mat), ncol(count_mat)))

df_allele <- data.table::fread(allele_file, showProgress = FALSE)
message(sprintf("[%s] Alleles: %d rows", ts(), nrow(df_allele)))

data(ref_hca)
message(sprintf("[%s] Reference loaded", ts()))
message(sprintf("[%s]", ts()))

message(sprintf("[%s] ========================================", ts()))
message(sprintf("[%s] CALLING: run_numbat()", ts()))
message(sprintf("[%s] ========================================", ts()))
message(sprintf("[%s]", ts()))

tryCatch({
  # Try with timeout - wrapped in system call timeout
  out_dir <- args[5]
  
  message(sprintf("[%s] Attempting run_numbat with ncores=%d...", ts(), ncores))
  message(sprintf("[%s] Output dir: %s", ts(), out_dir))
  
  # Direct call to run_numbat
  numbat_obj <- run_numbat(
    count_mat = count_mat,
    lambdas_ref = ref_hca,
    df_allele = df_allele,
    genome = "hg38",
    t = 1e-5,
    ncores = ncores,
    plot = TRUE,
    out_dir = out_dir,
    verbose = TRUE
  )
  
  message(sprintf("[%s] ✓ run_numbat() completed successfully!", ts()))
  
}, error = function(e) {
  message(sprintf("[%s] ✗ ERROR in run_numbat():", ts()))
  message(sprintf("[%s]   Message: %s", ts(), e$message))
  message(sprintf("[%s]   Call: %s", ts(), paste(deparse(sys.call(-1)), collapse=" ")))
  
  # Print full traceback
  message(sprintf("[%s]", ts()))
  message(sprintf("[%s] Full traceback:", ts()))
  traceback()
  
  quit(status = 1)
}, finally = {
  message(sprintf("[%s] run_numbat test finished", ts()))
})
