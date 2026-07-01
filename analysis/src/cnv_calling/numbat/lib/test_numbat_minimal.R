#!/usr/bin/env Rscript
# Minimal test - just load and try simplest call
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

message(sprintf("[%s] [MINIMAL TEST] Starting", ts()))

message(sprintf("[%s] Loading numbat...", ts()))
library(numbat)

message(sprintf("[%s] Loading data...", ts()))
count_mat <- readRDS("Data/04_analysis/cnv/numbat/inputs/lowseq_489/atac_bin/lowseq_489_atac_bin.rds")
message(sprintf("[%s] ATAC: %d × %d", ts(), nrow(count_mat), ncol(count_mat)))

df_allele <- data.table::fread("Data/04_analysis/cnv/numbat/inputs/lowseq_489/alleles/lowseq_489_atac_allele_counts.tsv.gz")
message(sprintf("[%s] Alleles: %d rows", ts(), nrow(df_allele)))

data(ref_hca)
message(sprintf("[%s] Reference loaded", ts()))

# Check what parameters run_numbat accepts
message(sprintf("[%s] Checking run_numbat signature...", ts()))
sig <- formals(run_numbat)
message(sprintf("[%s] Parameters: %s", ts(), paste(names(sig), collapse=", ")))
