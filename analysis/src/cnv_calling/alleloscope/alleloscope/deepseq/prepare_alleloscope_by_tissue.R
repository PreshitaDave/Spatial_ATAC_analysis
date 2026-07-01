source("/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/cnv_calling/alleloscope/alleloscope/alleloscope_prep_helpers.R")

project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1) args[[1]] else "deepseq"

if (!(dataset %in% c("deepseq", "lowseq"))) {
  stop("dataset must be one of: deepseq, lowseq", call. = FALSE)
}

# For deepseq, use retagged BAM directory with full 16bp CB tags
base_config <- list(
  sample_name = dataset,
  somatic_dir = file.path(project_root, "Data/variant_calling", dataset, "somatic"),
  vcf_dir = file.path(project_root, "Data/variant_calling", dataset, "germline"),
  bam_dir = file.path(project_root, "Data/variant_calling", dataset, "Bam_cb16"),
  fasta_file = file.path(project_root, "Data/hg38_resources/Homo_sapiens_assembly38.fasta"),
  vartrix_bin = "/projectnb/paxlab/presh/software/vartrix-1.1.22/vartrix_linux",
  chrom_sizes_file = file.path(project_root, "Data/hg38_resources/hg38.chrom.sizes"),
  bin_size = 1000000L,
  vartrix_threads = as.integer(Sys.getenv("NSLOTS", "4")),
  # BAMs are retagged to include full 16bp CB before VarTrix
  vartrix_barcode_length = 16L
)

for (tissue_name in c("488B", "489")) {
  message(sprintf("[%s] Starting tissue-specific run: %s %s", Sys.time(), dataset, tissue_name))

  cfg <- base_config
  # Keep sample_name as just the dataset (not tissue-specific) for BAM file lookups
  # BAMs are shared across tissues (deepseq_chr1.filter.bam, etc.)
  # Only the output directory, barcode file, and fragments file are tissue-specific
  cfg$output_dir <- file.path(project_root, "Data/04_analysis/cnv/alleloscope", sprintf("%s_%s", dataset, tissue_name))
  cfg$barcode_subset_file <- file.path(project_root, "Data/04_analysis/cnv/alleloscope/barcodes", sprintf("%s_%s.barcodes.tsv", dataset, tissue_name))
  cfg$fragments_file <- file.path(project_root, "Data/01_inputs/fragments", sprintf("%s_%s", dataset, tissue_name), sprintf("%s_%s.fragments.sort.filtered.bed.gz", dataset, tissue_name))

  prepare_alleloscope_inputs(cfg)

  message(sprintf("[%s] Completed tissue-specific run: %s %s", Sys.time(), dataset, tissue_name))
}
