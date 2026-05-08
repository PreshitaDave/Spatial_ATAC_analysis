#!/bin/bash

MONOPOGEN_PATH="/projectnb/paxlab/presh/software/Monopogen"
REF_FA="/projectnb/paxlab/presh/projects/spatial_atac/Data/hg38_resources/Homo_sapiens_assembly38.fasta"
MONOPOGEN_ENV="/projectnb/paxlab/presh/env/conda_env/monopgen_env"
SCRIPT_DIR="/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/monopgen/deepseq/somatic_chr"

cd "$SCRIPT_DIR"
mkdir -p qsub_logs_deep

echo "Generating somatic variant calling jobs for full chromosomes"
echo ""

job_count=0
for chr in chr{1..22}; do
    echo "Generating job for chromosome: $chr"
    
    qsub_script="qsub_${chr}.sh"
    
    cat > "$qsub_script" << EOFQSUB
#!/bin/bash
#\$ -P paxlab
#\$ -l h_rt=48:00:00
#\$ -N mono_${chr}
#\$ -l mem_per_core=8G
#\$ -pe omp 8
#\$ -o ${SCRIPT_DIR}/qsub_logs_deep/${chr}.log
#\$ -e ${SCRIPT_DIR}/qsub_logs_deep/${chr}.error.log

mkdir -p ${SCRIPT_DIR}/qsub_logs_deep
echo "Job started at \$(date)" >&2

export PATH="${MONOPOGEN_ENV}/bin:\$PATH"

REGION_FILE=\$(mktemp)
echo "${chr}" > "\$REGION_FILE"

python "${MONOPOGEN_PATH}/src/Monopogen.py" somatic -a "${MONOPOGEN_PATH}/apps" -r "\$REGION_FILE" -t 8 -i "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq" -l "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq_cell_data.csv" -s featureInfo -g "${REF_FA}" && \
python "${MONOPOGEN_PATH}/src/Monopogen.py" somatic -a "${MONOPOGEN_PATH}/apps" -r "\$REGION_FILE" -t 8 -i "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq" -l "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq_cell_data.csv" -s cellScan -g "${REF_FA}" && \
python "${MONOPOGEN_PATH}/src/Monopogen.py" somatic -a "${MONOPOGEN_PATH}/apps" -r "\$REGION_FILE" -t 8 -i "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq" -l "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq_cell_data.csv" -s LDrefinement -g "${REF_FA}"

EXIT_CODE=\$?
rm -f "\$REGION_FILE"
exit \$EXIT_CODE
EOFQSUB
    
    chmod +x "$qsub_script"
    ((job_count++))
done

echo ""
echo "Generated $job_count job scripts (one per chromosome: chr1-chr22)"
echo "To submit all jobs: cd $SCRIPT_DIR && for f in qsub_chr*.sh; do qsub \"\$f\"; done"
