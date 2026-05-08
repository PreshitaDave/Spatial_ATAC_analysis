#!/bin/bash -l

set -euo pipefail

DATASET="${DATASET:?DATASET env var is required}"
CHR="${CHR:?CHR env var is required}"
NSLOTS="${NSLOTS:-4}"

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
REF_FASTA="${PROJECT_ROOT}/Data/hg38_resources/Homo_sapiens_assembly38.fasta"
REF_DICT="${REF_FASTA%.fasta}.dict"
BAM="${PROJECT_ROOT}/Data/variant_calling/${DATASET}/Bam/${DATASET}_${CHR}.filter.bam"
OUT_DIR="${PROJECT_ROOT}/Data/variant_calling/${DATASET}/gatk_hc"
OUT_VCF="${OUT_DIR}/${CHR}.hc.vcf.gz"

mkdir -p "${OUT_DIR}"

echo "[$(date '+%F %T')] START GATK HaplotypeCaller"
echo "[$(date '+%F %T')] DATASET=${DATASET} CHR=${CHR} NSLOTS=${NSLOTS}"
echo "[$(date '+%F %T')] BAM=${BAM}"
echo "[$(date '+%F %T')] OUT_VCF=${OUT_VCF}"

if [[ ! -f "${BAM}" ]]; then
  echo "[$(date '+%F %T')] ERROR: BAM not found: ${BAM}" >&2
  exit 1
fi

module load miniconda/24.5.0
module load java/17.0.8
module load gatk/4.6.2.0
module load samtools

if [[ ! -f "${REF_FASTA}.fai" ]]; then
  echo "[$(date '+%F %T')] FASTA index missing; creating ${REF_FASTA}.fai"
  samtools faidx "${REF_FASTA}"
fi

if [[ ! -f "${REF_DICT}" ]]; then
  echo "[$(date '+%F %T')] FASTA dict missing; creating ${REF_DICT}"
  gatk CreateSequenceDictionary -R "${REF_FASTA}" -O "${REF_DICT}"
fi

if [[ ! -f "${BAM}.bai" ]]; then
  echo "[$(date '+%F %T')] BAM index missing; indexing now"
  samtools index -@ "${NSLOTS}" "${BAM}"
fi

TMP_BASE="${TMPDIR:-/scratch/${JOB_ID:-manual}.${SGE_TASK_ID:-1}.cds}"
mkdir -p "${TMP_BASE}"
BAM_FOR_GATK="${TMP_BASE}/${DATASET}_${CHR}.rg.bam"
echo "[$(date '+%F %T')] Creating temporary RG-fixed BAM for GATK: ${BAM_FOR_GATK}"
samtools addreplacerg \
  -m overwrite_all \
  -r "ID:${DATASET}" \
  -r "SM:${DATASET}" \
  -r "LB:${DATASET}" \
  -r "PL:ILLUMINA" \
  -r "PU:${DATASET}_${CHR}" \
  -o "${BAM_FOR_GATK}" \
  "${BAM}"
samtools index -@ "${NSLOTS}" "${BAM_FOR_GATK}"

echo "[$(date '+%F %T')] Running HaplotypeCaller for ${CHR}"
gatk --java-options "-Xms4g -Xmx24g -XX:ParallelGCThreads=${NSLOTS}" HaplotypeCaller \
  -R "${REF_FASTA}" \
  -I "${BAM_FOR_GATK}" \
  -L "${CHR}" \
  -O "${OUT_VCF}" \
  --native-pair-hmm-threads "${NSLOTS}" \
  --minimum-mapping-quality 20

gatk IndexFeatureFile -I "${OUT_VCF}"

echo "[$(date '+%F %T')] DONE ${DATASET} ${CHR}"
