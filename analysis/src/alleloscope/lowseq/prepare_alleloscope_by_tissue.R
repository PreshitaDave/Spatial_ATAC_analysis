source("/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/alleloscope/alleloscope_prep_helpers.R")

project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1) args[[1]] else "deepseq"

if (!(dataset %in% c("deepseq", "lowseq"))) {
  stop("dataset must be one of: deepseq, lowseq", call. = FALSE)
}

base_config <- list(
  sample_name = dataset,
  somatic_dir = file.path(project_root, "Data/variant_calling", dataset, "somatic"),
  vcf_dir = file.path(project_root, "Data/variant_calling", dataset, "germline"),
  bam_dir = file.path(project_root, "Data/variant_calling", dataset, "Bam"),
  fasta_file = file.path(project_root, "Data/hg38_resources/Homo_sapiens_assembly38.fasta"),
  vartrix_bin = "/projectnb/paxlab/presh/software/vartrix-1.1.22/vartrix_linux",
  chrom_sizes_file = file.path(project_root, "Data/hg38_resources/hg38.chrom.sizes"),
  fragments_file = file.path(project_root, "Data", sprintf("%s.fragments.sort.filtered.bed%s", dataset, if (dataset == "deepseq") ".gz" else "")),
  bin_size = 1000000L,
  vartrix_threads = as.integer(Sys.getenv("NSLOTS", "4"))
)

for (tissue_name in c("488B", "489")) {
  message(sprintf("[%s] Starting tissue-specific run: %s %s", Sys.time(), dataset, tissue_name))

  cfg <- base_config
  cfg$sample_name <- sprintf("%s_%s", dataset, tissue_name)
  cfg$output_dir <- file.path(project_root, "Data/alleloscope", sprintf("%s_%s", dataset, tissue_name))
  cfg$barcode_subset_file <- file.path(project_root, "Data/alleloscope/barcodes", sprintf("%s_%s.barcodes.tsv", dataset, tissue_name))

  prepare_alleloscope_inputs(cfg)

  message(sprintf("[%s] Completed tissue-specific run: %s %s", Sys.time(), dataset, tissue_name))
}
