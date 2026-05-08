#!/bin/bash
set -euo pipefail

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
IN_BAM_DIR="${IN_BAM_DIR:-${PROJECT_ROOT}/Data/variant_calling/deepseq/Bam}"
OUT_BAM_DIR="${OUT_BAM_DIR:-${PROJECT_ROOT}/Data/variant_calling/deepseq/Bam_cb16}"
THREADS="${THREADS:-${NSLOTS:-8}}"
CHROMS="${CHROMS:-$(seq -s ' ' 1 22)}"
FORCE_RETAG="${FORCE_RETAG:-0}"

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

find_chr_bam() {
  local chr="$1"
  local candidates=(
    "${IN_BAM_DIR}/deepseq_chr${chr}.filter.bam"
    "${IN_BAM_DIR}/chr${chr}.filter.targeted.bam"
    "${IN_BAM_DIR}/chr${chr}.filter.bam"
  )
  local f
  for f in "${candidates[@]}"; do
    if [[ -f "${f}" ]]; then
      printf '%s\n' "${f}"
      return 0
    fi
  done
  return 1
}

retag_one_chr() {
  local chr="$1"
  local in_bam out_bam tmp_bam

  in_bam="$(find_chr_bam "${chr}")" || {
    log "[warn] No input BAM found for chr${chr} in ${IN_BAM_DIR}; skipping"
    return 0
  }

  out_bam="${OUT_BAM_DIR}/deepseq_chr${chr}.filter.bam"
  tmp_bam="${out_bam}.tmp"

  if [[ "${FORCE_RETAG}" != "1" && -s "${out_bam}" && -s "${out_bam}.bai" ]]; then
    log "[resume] chr${chr}: reusing ${out_bam}"
    return 0
  fi

  mkdir -p "${OUT_BAM_DIR}"
  rm -f "${tmp_bam}" "${tmp_bam}.bai"

  log "[step] chr${chr}: retagging $(basename "${in_bam}") -> $(basename "${out_bam}")"

  samtools view -h "${in_bam}" \
    | awk 'BEGIN{FS=OFS="\t"}
           /^@/ {print; next}
           {
             ba=""; bb=""; cb_i=0
             for(i=12;i<=NF;i++) {
               if ($i ~ /^ba:Z:/) ba=substr($i,6)
               else if ($i ~ /^bb:Z:/) bb=substr($i,6)
               else if ($i ~ /^CB:Z:/) cb_i=i
             }

             if (ba != "" && bb != "") {
               cb="CB:Z:" ba bb
               if (cb_i > 0) {
                 $cb_i = cb
                 replaced++
               } else {
                 ++NF
                 $NF = cb
                 added++
               }
             } else if (cb_i > 0) {
               existing_only++
             } else {
               missing_all++
             }

             print
           }
           END {
             printf("retag_stats added=%d replaced=%d existing_only=%d missing_all=%d\n",
                    added+0, replaced+0, existing_only+0, missing_all+0) > "/dev/stderr"
           }' \
    | samtools view -b -@ "${THREADS}" -o "${tmp_bam}" -

  samtools index -@ "${THREADS}" "${tmp_bam}"
  mv -f "${tmp_bam}" "${out_bam}"
  mv -f "${tmp_bam}.bai" "${out_bam}.bai"

  local cb_count total_count
  cb_count=$(samtools view "${out_bam}" | awk 'BEGIN{n=0;cb=0} {for(i=12;i<=NF;i++){if($i~/^CB:Z:/){cb++;break}} n++; if(n>=50000) exit} END{print cb+0}' || true)
  total_count=$(samtools view "${out_bam}" | awk 'BEGIN{n=0} {n++; if(n>=50000) exit} END{print n+0}' || true)
  log "[done] chr${chr}: CB tags in sample=${cb_count}/${total_count}"
}

main() {
  if ! command -v samtools >/dev/null 2>&1; then
    log "[error] samtools not found in PATH"
    exit 127
  fi

  log "[start] Retag deepseq BAMs with CB=ba+bb"
  log "[start] IN_BAM_DIR=${IN_BAM_DIR}"
  log "[start] OUT_BAM_DIR=${OUT_BAM_DIR}"
  log "[start] THREADS=${THREADS}"
  log "[start] CHROMS=${CHROMS}"

  local chr
  for chr in ${CHROMS}; do
    retag_one_chr "${chr}"
  done

  log "[done] Retagging finished"
}

main "$@"
