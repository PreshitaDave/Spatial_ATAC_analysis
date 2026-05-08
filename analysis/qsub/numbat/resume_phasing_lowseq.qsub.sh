#!/bin/bash -l
#$ -P paxlab
#$ -N numbat_phase_low
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=8G
#$ -l h_rt=12:00:00
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_phase_low.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/numbat_phase_low.$JOB_ID.err
#$ -j n

set -euo pipefail

module load R
export PATH="/projectnb/paxlab/presh/env/calicost_env/bin:/projectnb/paxlab/presh/software/external/Eagle_v2.4.1:${PATH}"

DATASET=lowseq
ALLELE_OUTDIR=/projectnb/paxlab/presh/projects/spatial_atac/Data/numbat/inputs/alleles
ALLELE_DF="${ALLELE_OUTDIR}/${DATASET}_atac_allele_counts.tsv.gz"
NUMBAT_REPO=/projectnb/paxlab/presh/projects/spatial_atac/Data/numbat/numbat_repo
GMAP=/projectnb/paxlab/presh/software/external/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz
PILEUP_DIR="${ALLELE_OUTDIR}/pileup/${DATASET}_atac"
PHASE_SCRIPT="${ALLELE_OUTDIR}/run_phasing.sh"

echo "[$(date '+%F %T')] [start] resume phasing for dataset=${DATASET}"

if [[ -f "${ALLELE_DF}" ]]; then
  echo "[$(date '+%F %T')] [skip] allele counts already exist: ${ALLELE_DF}"
  exit 0
fi

# -- Step 1: Eagle phasing (uses patched run_phasing.sh) --
echo "[$(date '+%F %T')] [step] Running Eagle phasing via ${PHASE_SCRIPT}"
bash "${PHASE_SCRIPT}" 2>&1 | tee "${ALLELE_OUTDIR}/phasing.log"

# Verify all 22 phased VCFs were created
missing=0
for chr in $(seq 1 22); do
  f="${ALLELE_OUTDIR}/phasing/${DATASET}_chr${chr}.phased.vcf.gz"
  if [[ ! -f "${f}" ]]; then
    echo "[$(date '+%F %T')] [error] missing phased VCF for chr${chr}"
    missing=$((missing + 1))
  fi
done
if [[ "${missing}" -gt 0 ]]; then
  echo "[$(date '+%F %T')] [error] ${missing} phased VCFs missing after Eagle run"
  exit 1
fi
echo "[$(date '+%F %T')] [step] All 22 phased VCFs present"

# -- Step 2: Generate allele count dataframe from pileup + phased VCFs --
echo "[$(date '+%F %T')] [step] Generating allele count dataframe"
Rscript - <<'REOF'
.libPaths(c("/projectnb/paxlab/presh/Rlibs/4.5", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(vcfR)
  library(Matrix)
  library(numbat)
})

label      <- "lowseq"
sample     <- "lowseq_atac"
outdir     <- "/projectnb/paxlab/presh/projects/spatial_atac/Data/numbat/inputs/alleles"
pu_dir     <- file.path(outdir, "pileup", sample)
gmap_file  <- "/projectnb/paxlab/presh/software/external/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"

message(sprintf("[%s] Reading genetic map", Sys.time()))
genetic_map <- fread(gmap_file) %>%
  setNames(c("CHROM", "POS", "rate", "cM")) %>%
  group_by(CHROM) %>%
  mutate(start = POS, end = c(POS[2:length(POS)], POS[length(POS)])) %>%
  ungroup()

message(sprintf("[%s] Reading phased VCFs for 22 chromosomes", Sys.time()))
vcf_phased <- lapply(1:22, function(chr) {
  vcf_file <- file.path(outdir, "phasing", sprintf("%s_chr%d.phased.vcf.gz", label, chr))
  fread(vcf_file, skip = "#CHROM") %>%
    rename(CHROM = `#CHROM`) %>%
    mutate(CHROM = str_remove(CHROM, "chr"))
}) %>%
  Reduce(rbind, .) %>%
  mutate(CHROM = factor(CHROM, unique(CHROM)))

message(sprintf("[%s] Reading pileup VCF and count matrices", Sys.time()))
vcf_pu <- fread(file.path(pu_dir, "cellSNP.base.vcf"), skip = "#CHROM") %>%
  rename(CHROM = `#CHROM`) %>%
  mutate(CHROM = str_remove(CHROM, "chr"))

AD <- readMM(file.path(pu_dir, "cellSNP.tag.AD.mtx"))
DP <- readMM(file.path(pu_dir, "cellSNP.tag.DP.mtx"))
cell_barcodes <- fread(file.path(pu_dir, "cellSNP.samples.tsv"), header = FALSE) %>% pull(V1)

message(sprintf("[%s] Running preprocess_allele", Sys.time()))
df <- numbat:::preprocess_allele(
  sample      = label,
  vcf_pu      = vcf_pu,
  vcf_phased  = vcf_phased,
  AD          = AD,
  DP          = DP,
  barcodes    = cell_barcodes,
  gtf         = gtf_hg38,
  gmap        = genetic_map
) %>%
  filter(GT %in% c("1|0", "0|1"))

out_file <- file.path(outdir, sprintf("%s_allele_counts.tsv.gz", sample))
message(sprintf("[%s] Writing %d allele count rows to %s", Sys.time(), nrow(df), out_file))
fwrite(df, out_file, sep = "\t")
message(sprintf("[%s] Done", Sys.time()))
REOF

echo "[$(date '+%F %T')] [done] Allele counts written to ${ALLELE_DF}"

# -- Step 3: Run NUMBAT ATAC-bin analysis --
echo "[$(date '+%F %T')] [step] Running NUMBAT ATAC-bin analysis"
export DATASET ALLELE_DF
export NCORES="${NSLOTS:-8}"
export NUMBAT_REPO
bash /projectnb/paxlab/presh/projects/spatial_atac/analysis/src/numbat/run_numbat_atac_bin.sh

echo "[$(date '+%F %T')] [done] NUMBAT lowseq phasing+analysis complete"
