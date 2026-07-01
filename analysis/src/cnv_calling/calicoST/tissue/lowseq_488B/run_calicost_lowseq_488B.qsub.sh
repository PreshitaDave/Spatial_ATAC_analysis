#!/bin/bash -l
#$ -P paxlab
#$ -N calicost_lowseq_488B
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=8G
#$ -l h_rt=24:00:00
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/calicost_lowseq_488B_$JOB_ID.log
#$ -m bea
#$ -M preshita@bu.edu

set -euo pipefail

TISSUE="lowseq_488B"
PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
CALICOST_SRC_DIR="${PROJECT_ROOT}/analysis/src/cnv_calling/calicoST"
OUTPUT_ROOT="${PROJECT_ROOT}/Data/04_analysis/cnv/calicoST/${TISSUE}"
PYTHON="/projectnb/paxlab/presh/env/calicost_env/bin/python3"
NORMAL_BARCODES="${PROJECT_ROOT}/analysis/binsize_comparison/normal_barcodes/${TISSUE}_normal_barcodes.csv"
ALLELE_FILE="${PROJECT_ROOT}/Data/04_analysis/cnv/numbat/inputs/${TISSUE}/alleles/${TISSUE}_atac_allele_counts.tsv.gz"

echo "=== CalicoST pipeline: ${TISSUE} (mosaicfield normals) ==="
echo "Job ID: ${JOB_ID}"
echo "Start: $(date)"
echo "Host: $(hostname)"

# Verify CalicoST environment
${PYTHON} -c "import sys; sys.path.insert(0,'/projectnb/paxlab/presh/software/CalicoST/src'); import calicost; print('CalicoST loaded from', calicost.__file__)"

# Verify normal barcodes file (generated from mosaicfield per_spot_tumor_purity.csv)
if [ ! -f "${NORMAL_BARCODES}" ]; then
    echo "ERROR: Normal barcodes file not found: ${NORMAL_BARCODES}"
    echo "Regenerate with: python tissue/lowseq_488B/extract_488B_normals.py"
    exit 1
fi
echo "Normal barcodes: $(wc -l < ${NORMAL_BARCODES}) spots"

# Change to script directory so relative config paths work
cd "${CALICOST_SRC_DIR}"

# ============================================================
# Step 1: Export ArchR TileMatrix + spatial coords
# (Requires R environment)
# ============================================================
echo ""
echo "--- Step 1: ArchR export (TileMatrix + spatial coords) ---"
echo "Start: $(date)"

INTER_DIR="${OUTPUT_ROOT}/intermediate"
if [ -f "${INTER_DIR}/archr_tilematrix.mtx" ]; then
    echo "Intermediate files already exist — skipping step 1"
else
    set +u
    for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
        if [[ -f "$profile_file" ]]; then . "$profile_file" 2>/dev/null || true; break; fi
    done
    set -u
    module load R
    Rscript --no-save --no-restore 1_export_archr_data.R "${TISSUE}" 5000
fi

echo "Step 1 complete: $(date)"

# ============================================================
# Step 2: Build CalicoST parsed_inputs
# ============================================================
echo ""
echo "--- Step 2: Building CalicoST parsed_inputs ---"
echo "Start: $(date)"

PARSED_INPUTS="${OUTPUT_ROOT}/parsed_inputs/table_bininfo.csv.gz"
if [ -f "${PARSED_INPUTS}" ]; then
    echo "parsed_inputs already exist — skipping step 2"
else
    ${PYTHON} 2_build_calicost_inputs.py "${TISSUE}" \
        --snps-per-bin 200 \
        --allele-file "${ALLELE_FILE}"
fi

echo "Step 2 complete: $(date)"

# Symlink parsed_inputs into purity/ and cna/ subdirs
mkdir -p "${OUTPUT_ROOT}/purity" "${OUTPUT_ROOT}/cna"

if [ ! -e "${OUTPUT_ROOT}/purity/parsed_inputs" ]; then
    ln -s "${OUTPUT_ROOT}/parsed_inputs" "${OUTPUT_ROOT}/purity/parsed_inputs"
    echo "Symlinked parsed_inputs → purity/"
fi

if [ ! -e "${OUTPUT_ROOT}/cna/parsed_inputs" ]; then
    ln -s "${OUTPUT_ROOT}/parsed_inputs" "${OUTPUT_ROOT}/cna/parsed_inputs"
    echo "Symlinked parsed_inputs → cna/"
fi

# ============================================================
# Step 3: Tumor purity estimation
# ============================================================
echo ""
echo "--- Step 3: Tumor purity estimation (mosaicfield normals) ---"
echo "Start: $(date)"

${PYTHON} 3_run_purity.py tissue/${TISSUE}/config_purity.yaml

echo "Step 3 complete: $(date)"

# ============================================================
# Step 4: CNA + clone calling
# ============================================================
echo ""
echo "--- Step 4: CNA + clone calling ---"
echo "Start: $(date)"

${PYTHON} 4_run_cna.py tissue/${TISSUE}/config_cna.yaml

echo "Step 4 complete: $(date)"

# ============================================================
# Step 5: Postprocessing and visualization
# ============================================================
echo ""
echo "--- Step 5: Postprocessing and visualization ---"
echo "Start: $(date)"

${PYTHON} 5_postprocess_visualize.py "${TISSUE}" --n-clones 3 --n-clones-purity 5

echo "Step 5 complete: $(date)"

echo ""
echo "=== CalicoST pipeline complete: ${TISSUE} ==="
echo "End: $(date)"
echo "Results: ${OUTPUT_ROOT}"
