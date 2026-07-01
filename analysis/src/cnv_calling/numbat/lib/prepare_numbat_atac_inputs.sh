#!/usr/bin/env bash
set -euo pipefail

on_err() {
  echo "[$(date '+%F %T')] ERROR: command failed at line ${1}: ${2}" >&2
}

trap 'on_err ${LINENO} "${BASH_COMMAND}"' ERR



# Ensure critical R packages are available before proceeding
_RPKG_SCRIPT="$(mktemp /tmp/numbat_install_XXXXXX.R)"
cat > "${_RPKG_SCRIPT}" << 'RPKGS'
required_packages <- c('GenomeInfoDb', 'IRanges', 'S4Vectors', 'BiocGenerics', 'numbat')
missing <- required_packages[!sapply(required_packages, function(pkg) requireNamespace(pkg, quietly=TRUE))]
if (length(missing) > 0) {
  cat("Installing missing packages:", paste(missing, collapse=', '), "\n")
  if (!requireNamespace('BiocManager', quietly=TRUE)) {
    install.packages('BiocManager', repos='http://cran.rstudio.com/')
  }
  bioc_pkgs <- intersect(missing, c('GenomeInfoDb', 'IRanges', 'S4Vectors', 'BiocGenerics'))
  if (length(bioc_pkgs) > 0) BiocManager::install(bioc_pkgs, ask=FALSE, update=FALSE)
  if ('numbat' %in% missing) {
    # Install GitHub dependencies then numbat from local source
    for (p in c('hahmmr', 'scistreer')) {
      if (!requireNamespace(p, quietly=TRUE)) devtools::install_github(paste0('kharchenkolab/', p), upgrade='never')
    }
    install.packages(
      file.path(Sys.getenv('NUMBAT_REPO', unset='/projectnb/paxlab/presh/projects/spatial_atac/Data/numbat/numbat_repo')),
      repos=NULL, type='source'
    )
  }
}
cat("All required R packages available\n")
RPKGS

if ! command -v Rscript >/dev/null 2>&1; then
  echo "[$(date '+%F %T')] ERROR: Rscript not found in PATH. Load R module first (e.g., module load R)." >&2
  exit 127
fi

Rscript "${_RPKG_SCRIPT}"
rm -f "${_RPKG_SCRIPT}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${SCRIPT_DIR}/numbat_common.sh"

DATASET="${DATASET:-lowseq}"
TISSUE="${TISSUE:-488B}"
NCORES="${NCORES:-4}"
CHROMS="${CHROMS:-${DEFAULT_CHROMS}}"

# Fragment file path is tissue-specific and organized in Data/01_inputs/fragments/{dataset}_{tissue}/
FRAG_FILE_DEFAULT="${PROJECT_ROOT}/Data/01_inputs/fragments/${DATASET}_${TISSUE}/${DATASET}_${TISSUE}.fragments.sort.filtered.bed.gz"
# Handle uncompressed files if gz doesn't exist
if [[ ! -f "${FRAG_FILE_DEFAULT}" && -f "${FRAG_FILE_DEFAULT%.gz}" ]]; then
  FRAG_FILE_DEFAULT="${FRAG_FILE_DEFAULT%.gz}"
fi

FRAG_FILE="${FRAG_FILE:-${FRAG_FILE_DEFAULT}}"
BARCODE_FILE_ATAC="${NUMBAT_INPUT_DIR}/${DATASET}_${TISSUE}_atac_barcodes.tsv"
BARCODE_FILE_PILEUP="${NUMBAT_INPUT_DIR}/${DATASET}_${TISSUE}_atac_barcodes_for_pileup.tsv"
AGG_ANNOT_FILE="${NUMBAT_INPUT_DIR}/${DATASET}_${TISSUE}_atac_barcodes_with_group.tsv"
BIN_GR="${NUMBAT_REF_DIR}/var220kb.rds"
ATAC_BIN_RDS="${NUMBAT_INPUT_DIR}/${DATASET}_${TISSUE}_atac_bin.rds"
ALLELE_OUTDIR="${NUMBAT_INPUT_DIR}/alleles"
ALLELE_DF="${ALLELE_OUTDIR}/${DATASET}_${TISSUE}_atac_allele_counts.tsv.gz"
PAR_FILE="${NUMBAT_REF_DIR}/par_numbatm.rds"
MERGED_BAM="${NUMBAT_INPUT_DIR}/${DATASET}_${TISSUE}_merged_for_numbat.bam"
MERGED_BAI="${MERGED_BAM}.bai"

PHASE_PANEL="${PHASE_PANEL:-${PROJECT_ROOT}/Data/hg38_resources/numbat/1000G_hg38}"
VCF_GENOME1K="${VCF_GENOME1K:-${PROJECT_ROOT}/Data/hg38_resources/numbat/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf}"
GMAP_GZ="${GMAP_GZ:-${PROJECT_ROOT}/Data/hg38_resources/numbat/genetic_map_hg38_withX.txt.gz}"

if [[ ! -d "${PHASE_PANEL}" && -d "${PROJECT_ROOT}/Data/variant_calling/phased_hg38_1000G_ref" ]]; then
  PHASE_PANEL="${PROJECT_ROOT}/Data/variant_calling/phased_hg38_1000G_ref"
fi

if [[ ! -f "${VCF_GENOME1K}" && -f "${PROJECT_ROOT}/Data/variant_calling/${DATASET}/germline/combined.phased.vcf.gz" ]]; then
  VCF_GENOME1K="${PROJECT_ROOT}/Data/variant_calling/${DATASET}/germline/combined.phased.vcf.gz"
fi

if [[ ! -f "${GMAP_GZ}" && -f "/projectnb/paxlab/presh/software/external/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz" ]]; then
  GMAP_GZ="/projectnb/paxlab/presh/software/external/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"
fi

NUMBAT_BIN=$(resolve_numbat_bin_dir)
GET_BINNED_ATAC="${NUMBAT_BIN}/get_binned_atac.R"
PILEUP_PHASE="${NUMBAT_BIN}/pileup_and_phase.R"

require_file "${GET_BINNED_ATAC}"
require_file "${PILEUP_PHASE}"
require_file "${FRAG_FILE}"

mkdir -p "${ALLELE_OUTDIR}"

log "Preparing NUMBAT ATAC-bin inputs for dataset=${DATASET}"
log "Using numbat scripts from ${NUMBAT_BIN}"
log "Phasing resources: panel=${PHASE_PANEL} snpvcf=${VCF_GENOME1K} gmap=${GMAP_GZ}"
log "Chromosomes: ${CHROMS}"

TARGETED_BAM_PATTERN="${PROJECT_ROOT}/Data/variant_calling/variant_calling/${DATASET}/Bam/%s.filter.targeted.bam"
FULL_BAM_PATTERN="${PROJECT_ROOT}/Data/variant_calling/variant_calling/${DATASET}/Bam/${DATASET}_%s.filter.bam"

# Tissue-specific barcode file from barcode-matching pipeline.
# TISSUE defaults to 488B; set TISSUE=489 to run on the other tissue.
TISSUE="${TISSUE:-488B}"

# Prefer BC16 barcodes (from monopogen/variant_calling) if available
# These already have the -1/-2 suffix matching BAM CB tags
BC16_BARCODE_FILE="${PROJECT_ROOT}/Data/01_inputs/barcodes/tissue_barcodes/${DATASET}_${TISSUE}/${DATASET}_${TISSUE}.bc16.no_edge_effect.barcodes.tsv"

if [[ -f "${BC16_BARCODE_FILE}" ]]; then
  TISSUE_BARCODE_FILE="${BC16_BARCODE_FILE}"
  BARCODE_FORMAT="BC16"
  log "Using BC16 barcodes (with -1/-2 suffix) for ${TISSUE}"
else
  # Fall back to BC8 barcodes (standard ATAC barcodes WITHOUT -1)
  TISSUE_BARCODE_FILE="${PROJECT_ROOT}/Data/01_inputs/barcodes/tissue_barcodes/${DATASET}_${TISSUE}/${DATASET}_${TISSUE}.no_edge_effect.barcodes.tsv"
  
  if [[ ! -f "${TISSUE_BARCODE_FILE}" ]]; then
    TISSUE_BARCODE_FILE="${PROJECT_ROOT}/Data/01_inputs/barcodes/tissue_barcodes/${DATASET}_${TISSUE}/${DATASET}_${TISSUE}.barcodes.tsv"
  fi
  BARCODE_FORMAT="BC8"
  log "Using BC8 barcodes (standard ATAC) for ${TISSUE}"
fi

mapfile -t TARGETED_BAMS < <(collect_existing_files "${TARGETED_BAM_PATTERN}")
mapfile -t FULL_BAMS < <(collect_existing_files "${FULL_BAM_PATTERN}")

if [[ -n "${ATAC_BAM:-}" ]]; then
  BAM_FILE="${ATAC_BAM}"
else
  if [[ "${#TARGETED_BAMS[@]}" -gt 0 ]]; then
    BAM_INPUTS=("${TARGETED_BAMS[@]}")
    log "Using chromosome-targeted BAMs for allele pileup"
  else
    BAM_INPUTS=("${FULL_BAMS[@]}")
    log "Targeted BAMs not found; falling back to full chromosome BAMs"
  fi
  require_nonempty_list "BAM inputs" "${BAM_INPUTS[@]}"
  BAM_FILE="${MERGED_BAM}"
fi

if [[ "${BAM_FILE}" == "${MERGED_BAM}" ]]; then
  if [[ ! -f "${MERGED_BAM}" || ! -f "${MERGED_BAI}" ]]; then
    log "Merging ${#BAM_INPUTS[@]} chromosome BAMs into ${MERGED_BAM}"
    bam_list=$(join_by " " "${BAM_INPUTS[@]}")
    run_cmd "samtools merge -f -@ '${NCORES}' '${MERGED_BAM}' ${bam_list}"
    run_cmd "samtools index -@ '${NCORES}' '${MERGED_BAM}'"
  else
    log "Reusing existing merged BAM ${MERGED_BAM}"
  fi
fi

require_file "${BAM_FILE}"

if [[ ! -f "${BIN_GR}" ]]; then
  run_cmd "curl -L -o '${BIN_GR}' 'https://raw.githubusercontent.com/kharchenkolab/numbat/main/inst/extdata/var220kb.rds'"
fi
require_file "${BIN_GR}"

# Build ATAC barcode file for fragment matrix:
# - For BC16: use as-is (already has -1/-2 suffix)
# - For BC8: append -1 to match fragment barcode format
require_file "${TISSUE_BARCODE_FILE}"
log "Using tissue barcode file for ${TISSUE}: ${TISSUE_BARCODE_FILE} ($(wc -l < "${TISSUE_BARCODE_FILE}") barcodes)"

if [[ "${BARCODE_FORMAT}" == "BC16" ]]; then
  # BC16 already has the -1/-2 suffix
  cat "${TISSUE_BARCODE_FILE}" | sort -u > "${BARCODE_FILE_ATAC}"
  log "BC16 barcodes already have suffix format, using as-is"
else
  # BC8 needs -1 suffix appended
  sed 's/$/-1/' "${TISSUE_BARCODE_FILE}" | sort -u > "${BARCODE_FILE_ATAC}"
  log "BC8 barcodes appended with -1 suffix"
fi

require_file "${BARCODE_FILE_ATAC}"

# Build pileup barcode file to match BAM CB tag format. Some BAMs store CB
# without -1, and passing the wrong format to cellsnp-lite yields 0 variants.
if ! command -v samtools >/dev/null 2>&1; then
  log "ERROR: samtools is required to detect BAM CB tag format"
  exit 1
fi

tmp_bam_cb=$(mktemp /tmp/numbat_cb_bam_XXXXXX.txt)
tmp_raw_cb=$(mktemp /tmp/numbat_cb_raw_XXXXXX.txt)
tmp_plus_cb=$(mktemp /tmp/numbat_cb_plus_XXXXXX.txt)

# Sample first 100k alignments to infer CB format.
# With `set -o pipefail`, samtools may return 141 (SIGPIPE) because `head`
# exits early; treat that as expected for this sampling pipeline.
set +o pipefail
samtools view "${BAM_FILE}" | head -n 100000 | awk 'BEGIN{FS="\t"} {for(i=12;i<=NF;i++) if($i ~ /^CB:Z:/){print substr($i,6); break}}' | sort -u > "${tmp_bam_cb}"
samtools_status=${PIPESTATUS[0]}
set -o pipefail

if [[ "${samtools_status}" -ne 0 && "${samtools_status}" -ne 141 ]]; then
  log "ERROR: samtools view failed while sampling CB tags (exit=${samtools_status})"
  rm -f "${tmp_bam_cb}" "${tmp_raw_cb}" "${tmp_plus_cb}"
  exit 1
fi

if [[ ! -s "${tmp_bam_cb}" ]]; then
  log "ERROR: could not sample any CB tags from BAM ${BAM_FILE}"
  rm -f "${tmp_bam_cb}" "${tmp_raw_cb}" "${tmp_plus_cb}"
  exit 1
fi
sort -u "${TISSUE_BARCODE_FILE}" > "${tmp_raw_cb}"
sort -u "${BARCODE_FILE_ATAC}" > "${tmp_plus_cb}"

raw_overlap=$(comm -12 "${tmp_bam_cb}" "${tmp_raw_cb}" | wc -l)
plus_overlap=$(comm -12 "${tmp_bam_cb}" "${tmp_plus_cb}" | wc -l)
log "BAM barcode overlap (raw vs -1): ${raw_overlap} vs ${plus_overlap}"

if [[ "${raw_overlap}" -eq 0 && "${plus_overlap}" -eq 0 ]]; then
  log "ERROR: no overlap between BAM CB tags and tissue barcode list in either format"
  rm -f "${tmp_bam_cb}" "${tmp_raw_cb}" "${tmp_plus_cb}"
  exit 1
fi

if [[ "${raw_overlap}" -ge "${plus_overlap}" ]]; then
  cp "${tmp_raw_cb}" "${BARCODE_FILE_PILEUP}"
  log "Using no-suffix barcodes for pileup (matches BAM CB tags)"
else
  cp "${tmp_plus_cb}" "${BARCODE_FILE_PILEUP}"
  log "Using -1-suffixed barcodes for pileup (matches BAM CB tags)"
fi

rm -f "${tmp_bam_cb}" "${tmp_raw_cb}" "${tmp_plus_cb}"
require_file "${BARCODE_FILE_PILEUP}"

# aggregate_counts() requires a two-column annotation table with columns
# named exactly 'cell' (barcode) and 'group' (cell-type label).
# See numbat R/utils.R line 72: setNames(.$group, .$cell)
printf 'cell\tgroup\n' > "${AGG_ANNOT_FILE}"
awk '{print $1"\tall_cells"}' "${BARCODE_FILE_ATAC}" >> "${AGG_ANNOT_FILE}"
require_file "${AGG_ANNOT_FILE}"

if [[ ! -f "${ATAC_BIN_RDS}" ]]; then
  log "Generating ATAC bin-by-cell matrix..."
  Rscript "${GET_BINNED_ATAC}" --CB "${BARCODE_FILE_ATAC}" --frag "${FRAG_FILE}" --binGR "${BIN_GR}" --outFile "${ATAC_BIN_RDS}"
else
  log "Reusing existing ATAC bin matrix ${ATAC_BIN_RDS}"
fi
require_file "${ATAC_BIN_RDS}"

if [[ ! -f "${NUMBAT_REF_DIR}/lambdas_ATAC_bincnt.rds" ]]; then
  log "Generating aggregated ATAC reference..."
  Rscript "${GET_BINNED_ATAC}" --CB "${AGG_ANNOT_FILE}" --frag "${FRAG_FILE}" --binGR "${BIN_GR}" --outFile "${NUMBAT_REF_DIR}/lambdas_ATAC_bincnt.rds" --generateAggRef
else
  log "Reusing existing aggregated ATAC reference ${NUMBAT_REF_DIR}/lambdas_ATAC_bincnt.rds"
fi
require_file "${NUMBAT_REF_DIR}/lambdas_ATAC_bincnt.rds"

if [[ ! -d "${PHASE_PANEL}" || ! -f "${VCF_GENOME1K}" || ! -f "${GMAP_GZ}" ]]; then
  log "WARNING: one or more phasing resources are missing; allele dataframe step will fail unless you provide PHASE_PANEL, VCF_GENOME1K, GMAP_GZ."
fi

if [[ ! -f "${ALLELE_DF}" ]]; then
  log "Running pileup and phasing..."
  Rscript "${PILEUP_PHASE}" \
    --label "${DATASET}_${TISSUE}" \
    --samples "${DATASET}_${TISSUE}_atac" \
    --bams "${BAM_FILE}" \
    --barcodes "${BARCODE_FILE_PILEUP}" \
    --gmap "${GMAP_GZ}" \
    --snpvcf "${VCF_GENOME1K}" \
    --paneldir "${PHASE_PANEL}" \
    --ncores "${NCORES}" \
    --cellTAG CB \
    --UMItag None \
    --outdir "${ALLELE_OUTDIR}"

  pileup_vcf="${ALLELE_OUTDIR}/pileup/${DATASET}_${TISSUE}_atac/cellSNP.base.vcf"
  require_file "${pileup_vcf}"
  variant_rows=$(grep -vc '^#' "${pileup_vcf}" || true)
  if [[ "${variant_rows}" -eq 0 ]]; then
    log "ERROR: pileup produced 0 variant rows (${pileup_vcf}); barcode/BAM mismatch likely"
    exit 1
  fi
else
  log "Reusing existing allele counts ${ALLELE_DF}"
fi

log "Prepared files:"
log "  ${BAM_FILE}"
log "  ${BARCODE_FILE_ATAC}"
log "  ${BARCODE_FILE_PILEUP}"
log "  ${ATAC_BIN_RDS}"
log "  ${NUMBAT_REF_DIR}/lambdas_ATAC_bincnt.rds"
log "  ${ALLELE_DF}"
