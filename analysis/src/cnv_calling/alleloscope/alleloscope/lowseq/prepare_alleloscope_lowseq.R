source("/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/alleloscope/alleloscope_prep_helpers.R")

vartrix_threads <- as.integer(Sys.getenv("VARTRIX_THREADS", Sys.getenv("NSLOTS", "8")))

message(sprintf(
  "[%s] Launching lowseq Alleloscope prep with VarTrix threads=%d",
  Sys.time(),
  vartrix_threads
))

config <- list(
  sample_name = "lowseq",
  somatic_dir = "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/lowseq/somatic",
  vcf_dir = "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/lowseq/gatk_hc",
  vcf_pattern = "^chr([1-9]|1[0-9]|2[0-2])\\.hc\\.vcf\\.gz$",
  bam_dir = "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/lowseq/Bam",
  fasta_file = "/projectnb/paxlab/presh/projects/spatial_atac/Data/hg38_resources/Homo_sapiens_assembly38.fasta",
  vartrix_bin = "/projectnb/paxlab/presh/software/vartrix-1.1.22/vartrix_linux",
  chrom_sizes_file = "/projectnb/paxlab/presh/projects/spatial_atac/Data/hg38_resources/hg38.chrom.sizes",
  fragments_file = "/projectnb/paxlab/presh/projects/spatial_atac/Data/lowseq.fragments.sort.filtered.bed",
  output_dir = "/projectnb/paxlab/presh/projects/spatial_atac/Data/alleloscope/lowseq",
  bin_size = 1000000L,
  vartrix_threads = vartrix_threads
)

result <- prepare_alleloscope_inputs(config)
print(result)

# Now creating the object 
library(Alleloscope) # load
allelo.path = '/projectnb/paxlab/presh/software/Alleloscope'
data.path = '/projectnb/paxlab/presh/projects/spatial_atac/Data/alleloscope/lowseq'
dir_path <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/alleloscope/lowseq/output/"; dir.create(dir_path) # set up output directory

size=read.table(file.path(allelo.path, "data-raw/sizes.cellranger-GRCh38-1.0.0.txt"), stringsAsFactors = F) # read size file
size=size[1:22,]

### Reading data related files 
# SNP by cell matrices for ref and alt alleles
barcodes=read.table(file.path(data.path, "barcodes.tsv"), sep='\t', stringsAsFactors = F, header=F)
alt_all=readMM(file.path(data.path, "alt_all.mtx"))

ref_all=readMM(file.path(data.path, "ref_all.mtx"))
var_all=read.table(file.path(data.path, "var_all.vcf"), header = F, sep='\t', stringsAsFactors = F)

# bin by cell matrices for tumor and normal for segmentation
raw_counts=read.table(file.path(data.path, 'chr1000k_fragments.tsv'), sep='\t', header=T, row.names = 1,stringsAsFactors = F)
# colnames(raw_counts)=gsub("[.]","-", colnames(raw_counts))


Obj=Createobj(alt_all =alt_all, ref_all = ref_all, var_all = var_all ,samplename='lowseq', genome_assembly="GRCh38", dir_path=dir_path, barcodes=barcodes, size=size, assay='scATACseq')
Obj_filtered=Matrix_filter(Obj=Obj, cell_filter=5, SNP_filter=5, min_vaf = 0.1, max_vaf = 0.9) 

# suggest setting min_vaf=0.1 and max_vaf=0.9 when SNPs are called in the tumor sample for higher confident SNPs
Obj_filtered=Est_regions(Obj_filtered = Obj_filtered, max_nSNP = 30000, plot_stat = T,cont = FALSE)

# Recommend max_nSNP <50000
# Regions without allelic imbalence do not coverge (Reach the max number of iterations.)


Obj_filtered$select_normal$barcode_normal=cell_type[which(cell_type[,2]!='tumor'),1]



