#!/usr/bin/env Rscript

message("[DIAGNOSTIC TEST] Starting...")
message(sprintf("[%s] R version: %s", Sys.time(), R.version$version.string))

message("[DIAGNOSTIC TEST] Loading numbat library...")
tryCatch({
  library(numbat)
  message("[DIAGNOSTIC TEST] ✓ numbat loaded successfully")
}, error = function(e) {
  message(sprintf("[DIAGNOSTIC TEST] ✗ ERROR loading numbat: %s", e$message))
  quit(status = 1)
})

message("[DIAGNOSTIC TEST] Loading data...")
tryCatch({
  data(ref_hca)
  message("[DIAGNOSTIC TEST] ✓ ref_hca loaded successfully")
}, error = function(e) {
  message(sprintf("[DIAGNOSTIC TEST] ✗ ERROR loading ref_hca: %s", e$message))
  quit(status = 1)
})

message("[DIAGNOSTIC TEST] Loading test ATAC matrix...")
tryCatch({
  atac_file <- "Data/04_analysis/cnv/numbat/inputs/lowseq_489/atac_bin/lowseq_489_atac_bin.rds"
  count_mat <- readRDS(atac_file)
  message(sprintf("[DIAGNOSTIC TEST] ✓ ATAC loaded: %d x %d", nrow(count_mat), ncol(count_mat)))
}, error = function(e) {
  message(sprintf("[DIAGNOSTIC TEST] ✗ ERROR loading ATAC: %s", e$message))
  quit(status = 1)
})

message("[DIAGNOSTIC TEST] ✓ All tests passed!")
