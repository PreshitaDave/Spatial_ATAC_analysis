#!/bin/bash -l
#$ -P paxlab
#$ -N pileup_deep_488B
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=8G
#$ -l h_rt=24:00:00
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/pileup_deepseq_488B_$JOB_ID.log
#$ -m bea
#$ -M preshita@bu.edu

set -euo pipefail

set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
    if [[ -f "$profile_file" ]]; then . "$profile_file" 2>/dev/null || true; break; fi
done
set -u

module load R
module load samtools

# Eagle must be on PATH for pileup_and_phase.R
export PATH="/projectnb/paxlab/presh/env/calicost_env/bin:/projectnb/paxlab/presh/software/external/Eagle_v2.4.1:${PATH}"

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
BAM_CB16_DIR="${PROJECT_ROOT}/Data/04_analysis/monopogen_variant_calling/deepseq/Bam_cb16"
ALLELE_OUTDIR="${PROJECT_ROOT}/Data/04_analysis/cnv/numbat/inputs/deepseq_488B/alleles"
MERGED_BAM="${PROJECT_ROOT}/Data/04_analysis/cnv/numbat/inputs/deepseq_488B/bam/deepseq_488B_cb16_merged.bam"
BARCODES="${PROJECT_ROOT}/Data/01_inputs/barcodes/tissue_barcodes/deepseq_488B/deepseq_488B.barcodes.tsv"
ALLELE_DF="${ALLELE_OUTDIR}/deepseq_488B_atac_allele_counts.tsv.gz"
PILEUP_PHASE="/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/cnv/numbat/numbat_repo/inst/bin/pileup_and_phase.R"
GMAP="${PROJECT_ROOT}/Data/02_references/genome/hg38_resources/numbat/genetic_map_hg38_withX.txt.gz"
SNPVCF="${PROJECT_ROOT}/Data/02_references/genome/hg38_resources/numbat/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
PANELDIR="${PROJECT_ROOT}/Data/02_references/genome/hg38_resources/numbat/1000G_hg38"

echo "[$(date '+%F %T')] START pileup+phasing for deepseq_488B"
echo "Job ID: ${JOB_ID}, Host: $(hostname)"

# ============================================================
# Step 1: Merge Bam_cb16 chr BAMs → single BAM with 16bp CB tags
# The original deepseq_488B_merged_for_numbat.bam lacks CB tags entirely;
# Bam_cb16 contains all deepseq tissues but we filter by barcodes at pileup time.
# ============================================================
if [[ ! -f "${MERGED_BAM}" ]]; then
    echo "[$(date '+%F %T')] Merging Bam_cb16 chr1-chr22 BAMs..."
    BAM_LIST=$(ls \
        "${BAM_CB16_DIR}/deepseq_chr1.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr2.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr3.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr4.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr5.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr6.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr7.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr8.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr9.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr10.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr11.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr12.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr13.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr14.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr15.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr16.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr17.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr18.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr19.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr20.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr21.filter.bam" \
        "${BAM_CB16_DIR}/deepseq_chr22.filter.bam")
    samtools merge -f -@ "${NSLOTS:-8}" "${MERGED_BAM}" ${BAM_LIST}
    samtools index -@ "${NSLOTS:-8}" "${MERGED_BAM}"
    echo "[$(date '+%F %T')] Merge complete: ${MERGED_BAM}"
else
    echo "[$(date '+%F %T')] Reusing existing merged BAM: ${MERGED_BAM}"
fi

# Sanity check: CB tags present and 16bp
N_CB=$(samtools view "${MERGED_BAM}" chr1 | head -n 10000 | grep -c 'CB:Z:' || true)
echo "CB tag check (first 10k reads on chr1): ${N_CB} reads with CB"
SAMPLE_CB=$(samtools view "${MERGED_BAM}" chr1 2>/dev/null | head -n 10000 | grep -o 'CB:Z:[ACGT]*' | head -1 | cut -d: -f3 || true)
echo "Sample CB barcode: ${SAMPLE_CB} (length: ${#SAMPLE_CB})"
if [[ "${N_CB}" -eq 0 ]]; then
    echo "ERROR: no CB tags in merged BAM — check source BAMs"
    exit 1
fi
if [[ "${#SAMPLE_CB}" -ne 16 ]]; then
    echo "ERROR: CB barcode length is ${#SAMPLE_CB}, expected 16"
    exit 1
fi

# ============================================================
# Step 2: Pileup + phasing via pileup_and_phase.R
# Filters to deepseq_488B spots via barcodes TSV (16bp, no suffix)
# ============================================================
mkdir -p "${ALLELE_OUTDIR}"

if [[ ! -f "${ALLELE_DF}" ]]; then
    echo "[$(date '+%F %T')] Running pileup_and_phase.R..."
    Rscript "${PILEUP_PHASE}" \
        --label "deepseq_488B" \
        --samples "deepseq_488B_atac" \
        --bams "${MERGED_BAM}" \
        --barcodes "${BARCODES}" \
        --gmap "${GMAP}" \
        --snpvcf "${SNPVCF}" \
        --paneldir "${PANELDIR}" \
        --ncores "${NSLOTS:-8}" \
        --cellTAG CB \
        --UMItag None \
        --outdir "${ALLELE_OUTDIR}"
else
    echo "[$(date '+%F %T')] Allele counts already exist, skipping pileup."
fi

if [[ ! -f "${ALLELE_DF}" ]]; then
    echo "ERROR: allele counts not produced at ${ALLELE_DF}"
    exit 1
fi

echo "[$(date '+%F %T')] DONE pileup+phasing for deepseq_488B"
echo "Allele counts: ${ALLELE_DF}"
