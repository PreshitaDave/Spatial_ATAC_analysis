args <- commandArgs(trailingOnly = TRUE)
dataset_env <- Sys.getenv("DATASET", unset = "")
dataset <- if (nzchar(dataset_env)) {
  dataset_env
} else if (length(args) >= 1L) {
  args[[1]]
} else {
  "deepseq"
}

source("/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/alleloscope/alleloscope_prep_helpers.R")

somatic_dir <- file.path(
  "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling",
  dataset,
  "somatic"
)
vcf_dir <- file.path(
  "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling",
  dataset,
  "germline"
)
bam_dir <- file.path(
  "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling",
  dataset,
  "Bam"
)

barcodes <- read_barcode_order(somatic_dir)

tmp_dir <- tempfile(pattern = "alleloscope_smoke_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

vcf_in <- file.path(vcf_dir, "chr1.phased.vcf.gz")
vcf_sub <- file.path(tmp_dir, "chr1.first20.phased.vcf")
cmd <- sprintf(
  "zcat -f %s | awk 'BEGIN{n=0} /^##/{print; next} /^#CHROM/{print; next} {if(n<20){print; n++}}' > %s",
  shQuote(vcf_in),
  shQuote(vcf_sub)
)
system(cmd)

barcode_file <- file.path(tmp_dir, "barcodes.tsv")
writeLines(barcodes, barcode_file)

pair <- build_sparse_pair_from_vartrix(
  vcf_file = vcf_sub,
  bam_file = find_chr_bam(bam_dir = bam_dir, sample_name = dataset, chr = "chr1"),
  barcodes = barcodes,
  barcode_file = barcode_file,
  fasta_file = "/projectnb/paxlab/presh/projects/spatial_atac/Data/hg38_resources/Homo_sapiens_assembly38.fasta",
  vartrix_bin = "/projectnb/paxlab/presh/software/vartrix-1.1.22/vartrix_linux",
  work_dir = file.path(tmp_dir, "vartrix"),
  threads = 1L
)
vars <- build_var_table_from_vcf(vcf_sub)

cat("dataset", dataset, "\n")
cat("barcodes", length(barcodes), "\n")
cat("pair_dims", dim(pair$ref)[1], dim(pair$ref)[2], "\n")
cat("var_rows", nrow(vars), "\n")
cat("rows_match", pair$n_rows == nrow(vars), "\n")