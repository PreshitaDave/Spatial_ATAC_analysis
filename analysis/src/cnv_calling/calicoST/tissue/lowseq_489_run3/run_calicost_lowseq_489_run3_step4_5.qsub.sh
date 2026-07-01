#!/bin/bash -l
#$ -P paxlab
#$ -N calicost_489_r3
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=8G
#$ -l h_rt=24:00:00
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/calicost_lowseq_489_run3_step4_5_$JOB_ID.log
#$ -m bea
#$ -M preshita@bu.edu

set -euo pipefail

TISSUE="lowseq_489"
RUN="run3"
PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
CALICOST_SRC_DIR="${PROJECT_ROOT}/analysis/src/cnv_calling/calicoST"
OUTPUT_ROOT="${PROJECT_ROOT}/Data/04_analysis/cnv/calicoST/${TISSUE}/${RUN}"
PYTHON="/projectnb/paxlab/presh/env/calicost_env/bin/python3"

echo "=== CalicoST ${TISSUE} ${RUN}: steps 4+5 only ==="
echo "Job ID: ${JOB_ID}"
echo "Start: $(date)"
echo "Host: $(hostname)"

${PYTHON} -c "import sys; sys.path.insert(0,'/projectnb/paxlab/presh/software/CalicoST/src'); import calicost; print('CalicoST loaded from', calicost.__file__)"

cd "${CALICOST_SRC_DIR}"

# ============================================================
# Step 4: CNA + clone calling (purity complete; re-run CNA only)
# ============================================================
echo ""
echo "--- Step 4: CNA + clone calling ---"
echo "Start: $(date)"

${PYTHON} 4_run_cna.py tissue/${TISSUE}_run3/config_cna.yaml

echo "Step 4 complete: $(date)"

# ============================================================
# Step 5: Postprocessing and visualization
# ============================================================
echo ""
echo "--- Step 5: Postprocessing and visualization ---"
echo "Start: $(date)"

${PYTHON} 5_postprocess_visualize.py "${TISSUE}" --n-clones 3 --n-clones-purity 5 \
    --output-base "${OUTPUT_ROOT}"

echo "Step 5 complete: $(date)"

echo ""
echo "=== Done: ${TISSUE} ${RUN} ==="
echo "End: $(date)"
echo "Results: ${OUTPUT_ROOT}"
