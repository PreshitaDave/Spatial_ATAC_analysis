#!/usr/bin/env Rscript

# Script: get_binned_atac_fixed.R
# Description: Generate cell x counts matrix per genomic tiles from ATAC-seq fragments
# FIX: Handles barcode format mismatch (fragment barcodes may have -1 or -2 suffix)
# Usage: Rscript get_binned_atac_fixed.R --CB barcodes.txt --frag fragments.tsv.gz --binGR bins.rds --outFile output.tsv

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(rtracklayer)
})

ensure_matrix <- function(x) {
  if (inherits(x, "table")) {
    x <- unclass(as.matrix(x))
  } else if (inherits(x, "array")) {
    x <- as.matrix(x)
  }
  x
}

message('Get cell x counts per supplied GRanges tiles. Saves as tsv file.')
option_list = list(
  make_option('--CB', type = "character", default = NULL,
              help = "Cell barcodes file (txt). First column must be barcodes."),
  make_option('--frag', type = "character", default = NULL,
              help = "Path to fragments.tsv[.gz] file"),
  make_option('--binGR', type = "character", default = NULL,
              help = "Binned genome GRanges file (RDS format)"),
  make_option('--outFile', type = "character", default = NULL,
              help = "Output file path (.tsv)"),
  make_option('--generateAggRef', action = "store_true", default = FALSE,
              help = "Indicator for whether we are generating aggregated reference or just a count matrix")
)

args = parse_args(OptionParser(option_list = option_list))
invisible(list2env(args,environment()))

# Main function to process data and generate bin-by-cell matrix
generate_bin_cell_matrix <- function(barcode_file, fragment_file, bin_file, output_file) {
  # Read cell barcodes
  cat("Reading cell barcodes...\n")
  barcodes <- as.list(data.table::fread(barcode_file)[[1]])
  cat(paste0("Found ", length(barcodes), " cell barcodes.\n"))
  
  # Read fragment file
  cat("Reading fragment file...\n")
  fragments <- import.bed(fragment_file, extraCols = c("type" = "character", "type" = "integer"))
  colnames(mcols(fragments)) <- c("barcode", "dup_counts")
  cat(paste0("Total ", length(fragments), " fragments.\n"))
  
  # Filter fragments by cell barcodes
  cat("Filtering fragments by cell barcodes...\n")
  
  # FIX: Strip -1, -2, etc suffix from fragment barcodes to match cell barcode format
  # Fragment barcodes from 10X often have suffix like "ATCG-1", "ATCG-2", etc
  # Cell barcode file typically has just "ATCG"
  fragment_barcodes_raw <- fragments$barcode
  fragment_barcodes_clean <- gsub("-[0-9]+$", "", fragment_barcodes_raw)
  
  # Create mapping for later use
  fragments$barcode_clean <- fragment_barcodes_clean
  
  # Check how many fragments match after cleaning
  cat(paste0("  Fragment barcodes with suffix (sample): "))
  cat(paste0(head(fragment_barcodes_raw, 3), collapse=", "))
  cat("\n")
  cat(paste0("  Fragment barcodes after cleaning: "))
  cat(paste0(head(fragment_barcodes_clean, 3), collapse=", "))
  cat("\n")
  cat(paste0("  Cell barcodes (sample): "))
  cat(paste0(head(unlist(barcodes), 3), collapse=", "))
  cat("\n")
  
  # Filter using cleaned barcodes
  fragments_in_cell <- fragments[fragments$barcode_clean %in% barcodes]
  cat(paste0("Total ", length(fragments_in_cell), " fragments in cells.\n"))
  
  if (length(fragments_in_cell) == 0) {
    warning("NO FRAGMENTS MATCHED! Check barcode format compatibility")
  }
  
  # Read genomic bins
  cat("Reading genomic bins...\n")
  query <- readRDS(bin_file)
  cat(paste0("Using ", length(query), " genomic bins.\n"))
  
  # Find overlaps between bins and fragments
  cat("Finding overlaps between bins and fragments...\n")
  ov <- findOverlaps(query, fragments_in_cell)
  ov <- as.matrix(ov)
  tmp <- fragments_in_cell$barcode_clean[ov[,2]]  # Use cleaned barcodes
  ov <- cbind(ov, match(tmp, barcodes))
  
  # Generate bin-by-cell matrix
  cat("Generating bin-by-cell matrix...\n")
  mm <- table(ov[,1], ov[,3])
  colnames(mm) <- barcodes[as.numeric(colnames(mm))]
  bins <- paste0(as.character(seqnames(query)), ":", as.character(ranges(query)))
  rownames(mm) <- bins[as.numeric(rownames(mm))]
  mm <- as(ensure_matrix(mm), "dgCMatrix")
  if (args$generateAggRef) {
	  message("Generating aggregated reference...\n")
	  # Read annotations and add group column (all cells as one reference group)
	  annot <- data.table::fread(barcode_file)
	  colnames(annot) <- "barcode"
	  annot$group <- "reference"  # All cells belong to single reference group
	  
	  message("  Annotation data: ", nrow(annot), " cells")
	  message("  Matrix dimensions: ", nrow(mm), " bins x ", ncol(mm), " cells")
	  
	   agg_ref_counts = numbat::aggregate_counts(as.matrix(mm), 
						  annot, 
						  normalized = TRUE, 
						  verbose = TRUE)
	  message("Saving aggregated count matrix.")
	   if(endsWith(output_file, ".rds")) {
	     saveRDS(agg_ref_counts, output_file)
	   }else if( endsWith(output_file, ".tsv")) {
	     write.table(agg_ref_counts, output_file)
	   }
	  
  } else {
	  message("Proceeding with count matrix (not aggregated)\n")
	  # Write output
	  cat("Writing output to file...\n")
	  if(endsWith(output_file, ".rds")) {
	    saveRDS(mm, output_file)
	    cat(paste0("Success! Bin-by-cell matrix saved to: ", output_file, "\n"))
	  }else if( endsWith(output_file, ".tsv")) {
	    write.table(mm, output_file, sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
	    cat(paste0("Success! Bin-by-cell matrix saved to: ", output_file, "\n"))
	  }
	  
  }
}

# Run main function
generate_bin_cell_matrix(CB, frag, binGR, outFile)
