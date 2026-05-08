#!/bin/bash
# GENERALIZED SUBMIT SCRIPT FOR SOMATIC VARIANT CALLING
# Usage: bash submit_somatic_template.sh <dataset_name> <bam_data_path> <cell_data_csv>
# Example: bash submit_somatic_template.sh deepseq /projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq /projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq_cell_data.csv

# --- Configuration ---
DATASET_NAME="${1:-deepseq}"  # Default: deepseq
BAM_DATA_PATH="${2:-/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq}"
CELL_DATA_CSV="${3:-/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq_cell_data.csv}"

MONOPOGEN_PATH="/projectnb/paxlab/presh/software/Monopogen"
REF_FA="/projectnb/paxlab/presh/projects/spatial_atac/Data/hg38_resources/Homo_sapiens_assembly38.fasta"
REGION_FILE="/projectnb/paxlab/presh/software/Monopogen/resource/GRCh38.region.lst"
MONOPOGEN_ENV="/projectnb/paxlab/presh/env/conda_env/monopgen_env"
SCRIPT_DIR="/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/monopgen/${DATASET_NAME}/somatic_chr"

cd "$SCRIPT_DIR" || { echo "Cannot cd to $SCRIPT_DIR"; exit 1; }

mkdir -p qsub_logs_${DATASET_NAME}

echo "Generating somatic variant calling jobs for dataset: $DATASET_NAME"
echo "BAM path: $BAM_DATA_PATH"
echo "Cell data: $CELL_DATA_CSV"
echo ""

# --- Generate job scripts for each region ---
job_count=0
while IFS=',' read -r chr start end; do
    [[ -z "$chr" ]] && continue
    [[ "$chr" == "chr" ]] && continue  # Skip header if present
    
    region_id="${chr}:${start}-${end}"
    echo "Generating job for region: $region_id"
    
    qsub_script="qsub_${chr}.sh"
    
    cat > "$qsub_script" << EOFQSUB
#!/bin/bash
#\$ -P paxlab
#\$ -l h_rt=48:00:00
#\$ -N mono_${chr}
#\$ -l mem_per_core=8G
#\$ -pe omp 8
#\$ -o ${SCRIPT_DIR}/qsub_logs_${DATASET_NAME}/${region_id}.log
#\$ -e ${SCRIPT_DIR}/qsub_logs_${DATASET_NAME}/${region_id}.error.log

mkdir -p ${SCRIPT_DIR}/qsub_logs_${DATASET_NAME}

echo "========================================" >&2
echo "Job started at \$(date)" >&2
echo "Region: ${region_id}" >&2
echo "Dataset: ${DATASET_NAME}" >&2
echo "Node: \$(hostname)" >&2
echo "========================================" >&2

export PATH="${MONOPOGEN_ENV}/bin:\$PATH"

MONOPOGEN_PATH="${MONOPOGEN_PATH}"
REF_FA="${REF_FA}"
REGION="${region_id}"

# Create region file in proper comma-separated format
REGION_FILE=\$(mktemp)
echo "${chr},${start},${end}" > "\$REGION_FILE"

echo "Starting \$REGION at \$(date)" >&2
echo "Python: \$(which python)" >&2
echo "Python version: \$(python --version)" >&2
echo ""

# Step 1: featureInfo
echo "Running featureInfo..." >&2
python "\${MONOPOGEN_PATH}/src/Monopogen.py" somatic \\
    -a "\${MONOPOGEN_PATH}/apps" \\
    -r "\$REGION_FILE" \\
    -t 8 \\
    -i "${BAM_DATA_PATH}" \\
    -l "${CELL_DATA_CSV}" \\
    -s featureInfo \\
    -g "\${REF_FA}"

if [ \$? -ne 0 ]; then
    echo "featureInfo failed for \$REGION" >&2
    rm -f "\$REGION_FILE"
    exit 1
fi

# Step 2: cellScan
echo "Running cellScan..." >&2
python "\${MONOPOGEN_PATH}/src/Monopogen.py" somatic \\
    -a "\${MONOPOGEN_PATH}/apps" \\
    -r "\$REGION_FILE" \\
    -t 8 \\
    -i "${BAM_DATA_PATH}" \\
    -l "${CELL_DATA_CSV}" \\
    -s cellScan \\
    -g "\${REF_FA}"

if [ \$? -ne 0 ]; then
    echo "cellScan failed for \$REGION" >&2
    rm -f "\$REGION_FILE"
    exit 1
fi

# Step 3: LDrefinement
echo "Running LDrefinement..." >&2
python "\${MONOPOGEN_PATH}/src/Monopogen.py" somatic \\
    -a "\${MONOPOGEN_PATH}/apps" \\
    -r "\$REGION_FILE" \\
    -t 8 \\
    -i "${BAM_DATA_PATH}" \\
    -l "${CELL_DATA_CSV}" \\
    -s LDrefinement \\
    -g "\${REF_FA}"

EXIT_CODE=\$?
echo "" >&2
echo "Python script exit code: \$EXIT_CODE" >&2

rm -f "\$REGION_FILE"

if [ \$EXIT_CODE -eq 0 ]; then
    echo "Finished \$REGION successfully at \$(date)" >&2
else
    echo "ERROR: \$REGION failed with exit code \$EXIT_CODE at \$(date)" >&2
fi

echo "========================================" >&2
exit \$EXIT_CODE
EOFQSUB
    
    chmod +x "$qsub_script"
    ((job_count++))
    
done < "$REGION_FILE"

echo ""
echo "Generated $job_count job scripts in $SCRIPT_DIR"
echo ""
echo "To submit all jobs, run:"
echo "  cd $SCRIPT_DIR"
echo "  for f in qsub_*.sh; do qsub \"\$f\"; done"
