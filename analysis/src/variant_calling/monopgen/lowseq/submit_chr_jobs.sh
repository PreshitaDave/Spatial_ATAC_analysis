#!/bin/bash

MONOPOGEN_PATH="/projectnb/paxlab/presh/software/Monopogen"
REF_FA="/projectnb/paxlab/presh/projects/spatial_atac/Data/hg38_resources/Homo_sapiens_assembly38.fasta"
CHR_LIST="/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/chr_region.lst"
MONOPOGEN_ENV="/projectnb/paxlab/presh/env/conda_env/monopgen_env"
SCRIPT_DIR="/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/variant_calling/monopgen/lowseq/somatic_chr"

cd "$SCRIPT_DIR"

mkdir -p qsub_logs

while read chr; do
    [[ -z "$chr" ]] && continue
    
    echo "Generating job for $chr..."
    
    qsub_script="qsub_${chr}.sh"
    
    cat > "$qsub_script" << EOFQSUB
#!/bin/bash
#\$ -P paxlab
#\$ -l h_rt=48:00:00
#\$ -N mono_${chr}
#\$ -l mem_per_core=8G
#\$ -pe omp 8
#\$ -o ${SCRIPT_DIR}/qsub_logs/${chr}.log
#\$ -e ${SCRIPT_DIR}/qsub_logs/${chr}.error.log

# Create log directory if it doesn't exist
mkdir -p ${SCRIPT_DIR}/qsub_logs

echo "========================================" >&2
echo "Job started at \$(date)" >&2
echo "Chromosome: ${chr}" >&2
echo "Node: \$(hostname)" >&2
echo "========================================" >&2

# Export PATH to use monopgen_env Python
export PATH="${MONOPOGEN_ENV}/bin:\$PATH"

MONOPOGEN_PATH="${MONOPOGEN_PATH}"
REF_FA="${REF_FA}"
CHR="${chr}"

# Create temporary file with chromosome
CHR_FILE=\$(mktemp)
echo "\$CHR" > "\$CHR_FILE"

echo "Starting \$CHR at \$(date)" >&2
echo "Python: \$(which python)" >&2
echo "Python version: \$(python --version)" >&2
echo "Running command..." >&2

python "\${MONOPOGEN_PATH}/src/Monopogen.py" somatic \\
    -a  "\${MONOPOGEN_PATH}/apps" \\
    -r "\$CHR_FILE" \\
    -t 4 \\
    -i  "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/lowseq" \\
    -l  /projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/lowseq_cell_data.csv \\
    -s LDrefinement \\
    -g   "\${REF_FA}"

EXIT_CODE=\$?
echo "" >&2
echo "Python script exit code: \$EXIT_CODE" >&2

rm -f "\$CHR_FILE"

if [ \$EXIT_CODE -eq 0 ]; then
  echo "Finished \$CHR successfully at \$(date)" >&2
else
  echo "ERROR: \$CHR failed with exit code \$EXIT_CODE at \$(date)" >&2
fi

echo "========================================" >&2
exit \$EXIT_CODE
EOFQSUB
    
    chmod +x "$qsub_script"
    
    # COMMENTED OUT: Uncomment the line below to submit jobs to the cluster
    # job_id=$(qsub "$qsub_script")
    # echo "  → Job ID: $job_id"
    
    echo "  ✓ Script generated: $qsub_script (ready to submit)"
    
done < "$CHR_LIST"

echo ""
echo "✓ All jobs submitted!"