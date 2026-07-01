#!/bin/bash -l
#$ -P paxlab
#$ -N calicost_lowseq_489_run2
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=16G
#$ -l h_rt=12:00:00
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/calicost_lowseq_489_run2_$JOB_ID.log
#$ -m bea
#$ -M preshita@bu.edu

set -euo pipefail

TISSUE="lowseq_489"
RUN="run2"
PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
CALICOST_SRC_DIR="${PROJECT_ROOT}/analysis/src/cnv_calling/calicoST"
OUTPUT_ROOT="${PROJECT_ROOT}/Data/04_analysis/cnv/calicoST/${TISSUE}/${RUN}"
PARSED_INPUTS="${PROJECT_ROOT}/Data/04_analysis/cnv/calicoST/${TISSUE}/parsed_inputs"
PYTHON="/projectnb/paxlab/presh/env/calicost_env/bin/python3"
NORMAL_BARCODES="${PROJECT_ROOT}/analysis/binsize_comparison/normal_barcodes/${TISSUE}_normal_barcodes.csv"

echo "=== CalicoST pipeline: ${TISSUE} ${RUN} ==="
echo "Job ID: ${JOB_ID}"
echo "Start: $(date)"
echo "Host: $(hostname)"

# Verify environment
${PYTHON} -c "import sys; sys.path.insert(0,'/projectnb/paxlab/presh/software/CalicoST/src'); import calicost; print('CalicoST loaded from', calicost.__file__)"

# Verify normal barcodes file exists (must be generated before submitting this job)
if [ ! -f "${NORMAL_BARCODES}" ]; then
    echo "ERROR: Normal barcodes file not found: ${NORMAL_BARCODES}"
    echo "Run first: Rscript analysis/src/build_tissue/binsize_test_archr/9_extract_normal_barcodes_calicost.R ${TISSUE}"
    exit 1
fi
echo "Normal barcodes: $(wc -l < ${NORMAL_BARCODES}) spots"

# Change to script directory so relative config paths work
cd "${CALICOST_SRC_DIR}"

# ============================================================
# Step 2: parsed_inputs — reuse from run1 via symlink (skip rebuild)
# ============================================================
echo ""
echo "--- Step 2: Symlinking parsed_inputs from run1 ---"
mkdir -p "${OUTPUT_ROOT}/purity" "${OUTPUT_ROOT}/cna"

if [ ! -e "${OUTPUT_ROOT}/purity/parsed_inputs" ]; then
    ln -s "${PARSED_INPUTS}" "${OUTPUT_ROOT}/purity/parsed_inputs"
    echo "Symlinked parsed_inputs → run2/purity/"
fi

if [ ! -e "${OUTPUT_ROOT}/cna/parsed_inputs" ]; then
    ln -s "${PARSED_INPUTS}" "${OUTPUT_ROOT}/cna/parsed_inputs"
    echo "Symlinked parsed_inputs → run2/cna/"
fi

# ============================================================
# Step 3: Tumor purity estimation
# ============================================================
echo ""
echo "--- Step 3: Tumor purity estimation (run2, with normalidx_file) ---"
echo "Start: $(date)"

${PYTHON} 3_run_purity.py tissue/${TISSUE}_run2/config_purity.yaml

echo "Step 3 complete: $(date)"

# ============================================================
# Step 4: CNA + clone calling
# ============================================================
echo ""
echo "--- Step 4: CNA + clone calling ---"
echo "Start: $(date)"

# normalidx_file is already set in config_cna.yaml.
# tumorprop_file is None (ATAC-only: CalicoST uses BAF profiles).
${PYTHON} 4_run_cna.py tissue/${TISSUE}_run2/config_cna.yaml

echo "Step 4 complete: $(date)"

# ============================================================
# Step 5: Postprocess and visualize
# ============================================================
echo ""
echo "--- Step 5: Postprocessing and visualization ---"
echo "Start: $(date)"

${PYTHON} 5_postprocess_visualize.py "${TISSUE}" --n-clones 3 --n-clones-purity 5 \
    --output-base "${OUTPUT_ROOT}"

echo "Step 5 complete: $(date)"

# ============================================================
# Step 6: Enhanced plots (run2 paths)
# ============================================================
echo ""
echo "--- Step 6: Enhanced visualizations ---"
echo "Start: $(date)"

${PYTHON} 6_enhanced_visualize_run2.py "${TISSUE}"

echo "Step 6 complete: $(date)"

echo ""
echo "=== CalicoST pipeline complete: ${TISSUE} ${RUN} ==="
echo "End: $(date)"
echo "Results: ${OUTPUT_ROOT}"
