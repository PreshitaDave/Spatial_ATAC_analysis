#!/bin/bash -l
#$ -P paxlab
#$ -N calicost_lowseq_489
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=16G
#$ -l h_rt=12:00:00
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/calicost_lowseq_489_$JOB_ID.log
#$ -m bea
#$ -M preshita@bu.edu

set -euo pipefail

TISSUE="lowseq_489"
PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
CALICOST_SRC_DIR="${PROJECT_ROOT}/analysis/src/cnv_calling/calicoST"
OUTPUT_ROOT="${PROJECT_ROOT}/Data/04_analysis/cnv/calicoST/${TISSUE}"
PYTHON="/projectnb/paxlab/presh/env/calicost_env/bin/python3"

echo "=== CalicoST pipeline: ${TISSUE} ==="
echo "Job ID: ${JOB_ID}"
echo "Start: $(date)"
echo "Host: $(hostname)"

# Verify environment (conda env; call python binary directly, no activate needed)
${PYTHON} -c "import sys; sys.path.insert(0,'/projectnb/paxlab/presh/software/CalicoST/src'); import calicost; print('CalicoST loaded from', calicost.__file__)"

# Change to script directory so relative config paths work
cd "${CALICOST_SRC_DIR}"

# ============================================================
# Step 2: Build CalicoST parsed_inputs from ATAC + numbat data
# (Script 1 must have been run separately in the R environment)
# ============================================================
echo ""
echo "--- Step 2: Building CalicoST parsed_inputs ---"
echo "Start: $(date)"

${PYTHON} 2_build_calicost_inputs.py "${TISSUE}" --snps-per-bin 200

echo "Step 2 complete: $(date)"

# ============================================================
# Symlink parsed_inputs into purity/ and cna/ output dirs
# (CalicoST expects parsed_inputs/ under its output_dir)
# ============================================================
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
# Step 3: Run tumor purity estimation
# ============================================================
echo ""
echo "--- Step 3: Tumor purity estimation ---"
echo "Start: $(date)"

${PYTHON} 3_run_purity.py tissue/${TISSUE}/config_purity.yaml

echo "Step 3 complete: $(date)"

# ============================================================
# Step 4: Run CNA + clone calling
# ============================================================
echo ""
echo "--- Step 4: CNA + clone calling ---"
echo "Start: $(date)"

# Check that tumorprop file was produced
TUMORPROP="${OUTPUT_ROOT}/purity/clone5_rectangle0_w1.0/tumorprop_spots.tsv"
if [ ! -f "${TUMORPROP}" ]; then
    echo "WARNING: tumorprop_spots.tsv not found at expected path: ${TUMORPROP}"
    echo "Trying to find any tumorprop file..."
    TUMORPROP=$(find "${OUTPUT_ROOT}/purity" -name "tumorprop_spots.tsv" | head -1)
    if [ -z "${TUMORPROP}" ]; then
        echo "ERROR: No tumorprop file found. Running CNA without purity prior."
        # Update config to remove tumorprop_file
        sed -i "s|^tumorprop_file :.*|tumorprop_file : None|" \
            tissue/${TISSUE}/config_cna.yaml
    else
        echo "Found: ${TUMORPROP}"
        # Update config with actual path
        sed -i "s|^tumorprop_file :.*|tumorprop_file : ${TUMORPROP}|" \
            tissue/${TISSUE}/config_cna.yaml
    fi
fi

${PYTHON} 4_run_cna.py tissue/${TISSUE}/config_cna.yaml

echo "Step 4 complete: $(date)"

# ============================================================
# Step 5: Postprocess and visualize
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
