#!/bin/bash
# =============================================================================
# build_tissue_files.sh
#
# For each combination of {deepseq, lowseq} x {488B, 489}, produce:
#   1. barcodes file           (skip if exists)
#   2. edge_effect barcodes    (skip if exists; criterion: outermost 3 rows of tissue)
#   3. merged tissue BAM + index
#   4. tissue fragments file (bgzipped + tabix indexed)
#
# All outputs land in Data/tissue_barcodes/ (barcodes) and the per-dataset
# subdirectory under Data/variant_calling/{dataset}/tissue/ (BAMs & fragments).
#
# Environment / overrides
#   PROJECT_ROOT   path to repo root         (default: auto-detected)
#   DATASETS       space-separated list      (default: "deepseq lowseq")
#   TISSUES        space-separated list      (default: "488B 489")
#   FORCE          set to 1 to re-generate all files even if they exist
#   EDGE_N_ROWS    outer rows to treat as edge (default: 3)
#   CHROMS         chromosomes to merge BAMs from (default: 1..22)
#   THREADS        samtools threads           (default: $NSLOTS or 4)
#
# Dependencies
#   samtools 1.12+, bgzip (tabix), awk, python3 (standard library only)
#
# Checkpointing:  each output is skipped if it already exists and is non-empty
#   unless FORCE=1.
#
# Logging tags: [start] [step] [resume] [skip] [warn] [error] [done]
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../../.." && pwd)}"
DATA_DIR="${PROJECT_ROOT}/Data"
BARCODE_DIR="${DATA_DIR}/tissue_barcodes"
VARIANT_DIR="${DATA_DIR}/variant_calling"
TISSUE_POSITIONS="${DATA_DIR}/tissue_positions_list.csv"

# Source fragments from DriesSpatial D1942 by default (override via env if needed)
D1942_ROOT="${D1942_ROOT:-/projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D1942}"
DEEPSEQ_FRAG_SOURCE="${DEEPSEQ_FRAG_SOURCE:-${D1942_ROOT}/D1942_deepseq/D01942_NG06549_ATv008_1/fragments.sort.bed.gz}"
LOWSEQ_FRAG_SOURCE="${LOWSEQ_FRAG_SOURCE:-${D1942_ROOT}/D1942_lowseq/fragments.tsv.gz}"

DATASETS="${DATASETS:-deepseq lowseq}"
TISSUES="${TISSUES:-488B 489}"
FORCE="${FORCE:-0}"
EDGE_N_ROWS="${EDGE_N_ROWS:-3}"
THREADS="${THREADS:-${NSLOTS:-4}}"
CHROMS="${CHROMS:-$(seq 1 22)}"
BUILD_BAMS="${BUILD_BAMS:-1}"
BUILD_FRAGMENTS="${BUILD_FRAGMENTS:-1}"

_HTSLIB_BIN="/share/pkg.7/samtools/1.12/install/bin"
SAMTOOLS="${SAMTOOLS:-${_HTSLIB_BIN}/samtools}"
BGZIP="${BGZIP:-${_HTSLIB_BIN}/bgzip}"
TABIX="${TABIX:-${_HTSLIB_BIN}/tabix}"

LOG_DIR="${PROJECT_ROOT}/analysis/qsub_logs"
LOG_FILE="${LOG_DIR}/build_tissue_files.$(date +%Y%m%d_%H%M%S).log"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
log() {
    local tag="$1"; shift
    local msg="$*"
    local ts; ts="$(date '+%F %T')"
    echo "[${ts}] [${tag}] ${msg}" | tee -a "${LOG_FILE}"
}

skip_if_exists() {
    # Returns 0 (true) if file exists and non-empty (and FORCE!=1)
    local f="$1"
    if [[ "${FORCE}" != "1" && -s "${f}" ]]; then
        return 0
    fi
    return 1
}

require_file() {
    local f="$1" desc="$2"
    if [[ ! -f "${f}" ]]; then
        log "error" "Required file missing: ${f} (${desc})"
        exit 1
    fi
}

# ---------------------------------------------------------------------------
# Step 0: Validate dependencies
# ---------------------------------------------------------------------------
validate_deps() {
    log "start" "Validating dependencies"
    for cmd in "${SAMTOOLS}" "${BGZIP}" "${TABIX}" awk python3 gzip; do
        if ! command -v "${cmd}" >/dev/null 2>&1 && [[ ! -x "${cmd}" ]]; then
            log "error" "Command not found: ${cmd}"
            exit 1
        fi
    done
    require_file "${TISSUE_POSITIONS}" "tissue_positions_list.csv"
    require_file "${DEEPSEQ_FRAG_SOURCE}" "deepseq source fragments"
    require_file "${LOWSEQ_FRAG_SOURCE}" "lowseq source fragments"
    log "done" "All dependencies present"
}

# ---------------------------------------------------------------------------
# Step 1: Edge-effect barcodes via Python (outermost N rows of tissue)
# ---------------------------------------------------------------------------
build_edge_barcodes_python() {
    local bc_file="$1"        # input: tissue barcode list (16bp, no -1)
    local out_file="$2"       # output: edge barcodes
    local edge_n="${3:-${EDGE_N_ROWS}}"

    python3 - "${bc_file}" "${TISSUE_POSITIONS}" "${out_file}" "${edge_n}" <<'PY'
import sys, collections

bc_file, pos_file, out_file, edge_n = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])

# Load tissue barcodes
barcodes = set(open(bc_file).read().split())

# Load spatial positions: barcode -> (in_tissue, row, col, y, x)
positions = {}
for line in open(pos_file):
    p = line.strip().split(',')
    if len(p) >= 6:
        positions[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]), int(p[5]))

# Find min and max row for this tissue's barcodes
rows = [positions[bc][1] for bc in barcodes if bc in positions]
if not rows:
    print(f"WARNING: no barcodes with spatial positions found in {bc_file}", file=sys.stderr)
    open(out_file, 'w').close()
    sys.exit(0)

min_row, max_row = min(rows), max(rows)
edge_rows = set(range(min_row, min_row + edge_n)) | set(range(max_row - edge_n + 1, max_row + 1))

edge_barcodes = sorted(
    bc for bc in barcodes
    if bc in positions and positions[bc][1] in edge_rows
)

with open(out_file, 'w') as fh:
    fh.write('\n'.join(edge_barcodes))
    if edge_barcodes:
        fh.write('\n')

print(f"Edge barcodes: {len(edge_barcodes)} (rows {sorted(edge_rows)[:3]}...{sorted(edge_rows)[-3:]})", file=sys.stderr)
PY
}

# ---------------------------------------------------------------------------
# Step 2: Subset and merge per-chromosome BAMs into one tissue BAM
# ---------------------------------------------------------------------------
build_tissue_bam() {
    local dataset="$1"   # deepseq or lowseq
    local tissue="$2"    # 488B or 489
    local bc_file="$3"   # barcode file (16bp, no -1)

    local bam_dir out_dir out_bam tmp_bam barcode_arg

    # deepseq original BAMs use ba/bb tags — use Bam_cb16 retagged dir
    if [[ "${dataset}" == "deepseq" ]]; then
        bam_dir="${VARIANT_DIR}/deepseq/Bam_cb16"
        # BAM naming in Bam_cb16 uses deepseq_chr*.filter.bam
        bam_pattern="deepseq_chr{CHR}.filter.bam"
    else
        bam_dir="${VARIANT_DIR}/lowseq/Bam"
        bam_pattern="chr{CHR}.filter.targeted.bam"
    fi

    out_dir="${VARIANT_DIR}/${dataset}/tissue"
    mkdir -p "${out_dir}"
    out_bam="${out_dir}/${dataset}_${tissue}.bam"
    tmp_bam="${out_bam}.tmp_merge.bam"

    if skip_if_exists "${out_bam}" && skip_if_exists "${out_bam}.bai"; then
        log "resume" "${dataset} ${tissue}: BAM + index already exist: ${out_bam}"
        return 0
    fi

    log "step" "${dataset} ${tissue}: Building tissue BAM from per-chr BAMs"
    log "step" "  BAM dir: ${bam_dir}"
    log "step" "  Barcodes: $(wc -l < "${bc_file}") cells"

    # Build whitelist file for samtools view -D
    local whitelist; whitelist="${out_dir}/${dataset}_${tissue}.cb_whitelist.txt"
    cp "${bc_file}" "${whitelist}"

    # Collect per-chr BAMs that exist
    local chr_bams=()
    for chr in ${CHROMS}; do
        local bam_name="${bam_pattern//\{CHR\}/${chr}}"
        local bam_path="${bam_dir}/${bam_name}"
        if [[ -s "${bam_path}" ]]; then
            chr_bams+=("${bam_path}")
        else
            log "warn" "${dataset} ${tissue}: chr${chr} BAM not found at ${bam_path}, skipping"
        fi
    done

    if [[ ${#chr_bams[@]} -eq 0 ]]; then
        log "error" "${dataset} ${tissue}: No input BAMs found in ${bam_dir}"
        return 1
    fi

    log "step" "${dataset} ${tissue}: Merging ${#chr_bams[@]} chr BAMs, filtering by CB whitelist"

    # Subset each chr BAM to tissue barcodes then merge
    # samtools view -D CB:whitelist subsets by CB tag
    rm -f "${tmp_bam}"
    "${SAMTOOLS}" merge -f -@ "${THREADS}" --write-index \
        -o "${tmp_bam}" \
        <( for bam in "${chr_bams[@]}"; do
               "${SAMTOOLS}" view -b -@ 1 \
                   -D CB:"${whitelist}" \
                   "${bam}"
           done | cat )  2>>"${LOG_FILE}" || {
        # fallback: merge first then filter (for samtools versions without -D pipe support)
        log "warn" "${dataset} ${tissue}: pipe merge failed, trying sequential approach"
        rm -f "${tmp_bam}"
        local tmp_dir; tmp_dir="${out_dir}/.tmp_${dataset}_${tissue}"
        mkdir -p "${tmp_dir}"
        local filtered_bams=()
        for bam in "${chr_bams[@]}"; do
            local chr_name; chr_name="$(basename "${bam}" .bam)"
            local filt_bam="${tmp_dir}/${chr_name}.filt.bam"
            "${SAMTOOLS}" view -b -@ "${THREADS}" \
                -D CB:"${whitelist}" \
                -o "${filt_bam}" \
                "${bam}" 2>>"${LOG_FILE}"
            filtered_bams+=("${filt_bam}")
        done
        "${SAMTOOLS}" merge -f -@ "${THREADS}" "${tmp_bam}" "${filtered_bams[@]}" 2>>"${LOG_FILE}"
        rm -rf "${tmp_dir}"
    }

    "${SAMTOOLS}" sort -@ "${THREADS}" -o "${out_bam}" "${tmp_bam}" 2>>"${LOG_FILE}"
    "${SAMTOOLS}" index -@ "${THREADS}" "${out_bam}" 2>>"${LOG_FILE}"
    rm -f "${tmp_bam}"

    local read_count; read_count=$("${SAMTOOLS}" view -c "${out_bam}" 2>>"${LOG_FILE}")
    log "done" "${dataset} ${tissue}: BAM done — ${read_count} reads → ${out_bam}"
}

# ---------------------------------------------------------------------------
# Step 3: Subset fragments file to tissue barcodes
# ---------------------------------------------------------------------------
build_tissue_fragments() {
    local dataset="$1"  # deepseq or lowseq
    local tissue="$2"   # 488B or 489
    local bc_file="$3"  # barcode file (16bp, no -1)
    local label="${4:-}"  # optional suffix e.g. "edge_effect"

    local out_dir; out_dir="${VARIANT_DIR}/${dataset}/tissue"
    mkdir -p "${out_dir}"

    local out_suffix
    if [[ -n "${label}" ]]; then
        out_suffix="${label}.fragments.tsv.gz"
    else
        out_suffix="fragments.tsv.gz"
    fi
    local out_frag="${out_dir}/${dataset}_${tissue}.${out_suffix}"

    if skip_if_exists "${out_frag}" && skip_if_exists "${out_frag}.tbi"; then
        log "resume" "${dataset} ${tissue}: fragments${label:+ (${label})} already exist: ${out_frag}"
        return 0
    fi

    # Find full fragments source
    local src_frag
    if [[ "${dataset}" == "deepseq" ]]; then
        src_frag="${DEEPSEQ_FRAG_SOURCE}"
    else
        src_frag="${LOWSEQ_FRAG_SOURCE}"
    fi

    if [[ ! -s "${src_frag}" ]]; then
        log "error" "${dataset} ${tissue}: Source fragments not found: ${src_frag}"
        return 1
    fi

    log "step" "${dataset} ${tissue}: Subsetting fragments${label:+ (${label})} from ${src_frag}"
    log "step" "  Barcodes: $(wc -l < "${bc_file}") cells (fragments have -1 suffix)"

    # Fragments column 4 = barcode with "-1" suffix; bc_file = 16bp no suffix
    # Build awk hash of barcodes then filter col4 = bc"-1"
    local tmp_frag="${out_frag%.gz}.tmp"

    awk '
        BEGIN { OFS="\t" }
        NR==FNR { bc[$1"-1"] = 1; next }
        $4 in bc { print }
    ' "${bc_file}" \
            <(if [[ "${src_frag}" == *.gz ]]; then gzip -dc "${src_frag}"; else cat "${src_frag}"; fi) \
    > "${tmp_frag}"

    local frag_count; frag_count=$(wc -l < "${tmp_frag}")
    log "step" "${dataset} ${tissue}: ${frag_count} fragment lines matched"

    # bgzip and index
    "${BGZIP}" -@ "${THREADS}" -c "${tmp_frag}" > "${out_frag}"
    "${TABIX}" -p bed "${out_frag}"
    rm -f "${tmp_frag}"

    log "done" "${dataset} ${tissue}: Fragments${label:+ (${label})} done → ${out_frag}"
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
main() {
    mkdir -p "${LOG_DIR}"
    log "start" "build_tissue_files.sh"
    log "start" "PROJECT_ROOT=${PROJECT_ROOT}"
    log "start" "DATASETS=${DATASETS}"
    log "start" "TISSUES=${TISSUES}"
    log "start" "FORCE=${FORCE}"
    log "start" "EDGE_N_ROWS=${EDGE_N_ROWS}"
    log "start" "THREADS=${THREADS}"
    log "start" "BUILD_BAMS=${BUILD_BAMS} BUILD_FRAGMENTS=${BUILD_FRAGMENTS}"
    log "start" "DEEPSEQ_FRAG_SOURCE=${DEEPSEQ_FRAG_SOURCE}"
    log "start" "LOWSEQ_FRAG_SOURCE=${LOWSEQ_FRAG_SOURCE}"
    log "start" "LOG=${LOG_FILE}"

    validate_deps

    for dataset in ${DATASETS}; do
        for tissue in ${TISSUES}; do
            log "step" "===== ${dataset} ${tissue} ====="

            local bc_file="${BARCODE_DIR}/${dataset}_${tissue}.barcodes.tsv"
            local edge_file="${BARCODE_DIR}/${dataset}_${tissue}.edge_effect.barcodes.tsv"

            # --- 1. Barcodes (already built by build_tissue_barcode_lists.R) ---
            if [[ ! -s "${bc_file}" ]]; then
                log "error" "${dataset} ${tissue}: Barcodes file missing: ${bc_file}"
                log "error" "Run build_tissue_barcode_lists.R first to generate tissue barcode lists"
                continue
            fi
            log "resume" "${dataset} ${tissue}: Barcodes exist: $(wc -l < "${bc_file}") cells"

            # --- 2. Edge-effect barcodes ---
            if skip_if_exists "${edge_file}"; then
                log "resume" "${dataset} ${tissue}: Edge barcodes exist: $(wc -l < "${edge_file}") cells"
            else
                log "step" "${dataset} ${tissue}: Building edge_effect barcodes (outer ${EDGE_N_ROWS} rows)"
                build_edge_barcodes_python "${bc_file}" "${edge_file}" "${EDGE_N_ROWS}"
                log "done" "${dataset} ${tissue}: Edge barcodes → ${edge_file} ($(wc -l < "${edge_file}") cells)"
            fi

            # --- 3. Tissue BAM + index ---
            if [[ "${BUILD_BAMS}" == "1" ]]; then
                build_tissue_bam "${dataset}" "${tissue}" "${bc_file}"
            else
                log "skip" "${dataset} ${tissue}: Skipping BAM build (BUILD_BAMS=${BUILD_BAMS})"
            fi

            # --- 4. Tissue fragments (full barcodes) ---
            if [[ "${BUILD_FRAGMENTS}" == "1" ]]; then
                build_tissue_fragments "${dataset}" "${tissue}" "${bc_file}" ""
            else
                log "skip" "${dataset} ${tissue}: Skipping fragments build (BUILD_FRAGMENTS=${BUILD_FRAGMENTS})"
            fi


        done
    done

    log "done" "build_tissue_files.sh complete"
    log "done" "Outputs:"
    log "done" "  Barcodes + edge: ${BARCODE_DIR}/"
    log "done" "  BAMs + fragments: ${VARIANT_DIR}/{dataset}/tissue/"
}

main "$@"
