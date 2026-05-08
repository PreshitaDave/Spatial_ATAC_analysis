source("/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/alleloscope/alleloscope_prep_helpers.R")

config <- list(
  sample_name = "deepseq",
  somatic_dir = "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq/somatic",
  vcf_dir = "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq/germline",
  bam_dir = Sys.getenv(
    "DEEPSEQ_BAM_DIR",
    unset = "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq/Bam_cb16"
  ),
  fasta_file = "/projectnb/paxlab/presh/projects/spatial_atac/Data/hg38_resources/Homo_sapiens_assembly38.fasta",
  vartrix_bin = "/projectnb/paxlab/presh/software/vartrix-1.1.22/vartrix_linux",
  chrom_sizes_file = "/projectnb/paxlab/presh/projects/spatial_atac/Data/hg38_resources/hg38.chrom.sizes",
  fragments_file = "/projectnb/paxlab/presh/projects/spatial_atac/Data/deepseq.fragments.sort.filtered.bed.gz",
  output_dir = "/projectnb/paxlab/presh/projects/spatial_atac/Data/alleloscope/deepseq",
  bin_size = 1000000L,
  vartrix_threads = as.integer(Sys.getenv("NSLOTS", "8")),
  # BAMs are retagged to include full 16bp CB before VarTrix.
  vartrix_barcode_length = 16L
)

message(sprintf("[%s] Using deepseq BAM directory: %s", Sys.time(), config$bam_dir))
result <- prepare_alleloscope_inputs(config)
print(result)
print(result)