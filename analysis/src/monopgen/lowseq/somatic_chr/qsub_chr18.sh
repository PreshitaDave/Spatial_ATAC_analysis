#!/bin/bash
#$ -P paxlab
#$ -l h_rt=48:00:00
#$ -N mono_chr18
#$ -l mem_per_core=8G
#$ -pe omp 8
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/src/monopgen/lowseq/somatic_chr/qsub_logs/chr18.log
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/src/monopgen/lowseq/somatic_chr/qsub_logs/chr18.error.log

# Create log directory if it doesn't exist
mkdir -p /projectnb/paxlab/presh/projects/spatial_atac/analysis/src/monopgen/lowseq/somatic_chr/qsub_logs

echo "========================================" >&2
echo "Job started at $(date)" >&2
echo "Chromosome: chr18" >&2
echo "Node: $(hostname)" >&2
echo "========================================" >&2

# Export PATH to use monopgen_env Python
export PATH="/projectnb/paxlab/presh/env/conda_env/monopgen_env/bin:$PATH"

MONOPOGEN_PATH="/projectnb/paxlab/presh/software/Monopogen"
REF_FA="/projectnb/paxlab/presh/projects/spatial_atac/Data/hg38_resources/Homo_sapiens_assembly38.fasta"
CHR="chr18"

# Create temporary file with chromosome
CHR_FILE=$(mktemp)
echo "$CHR" > "$CHR_FILE"

echo "Starting $CHR at $(date)" >&2
echo "Python: $(which python)" >&2
echo "Python version: $(python --version)" >&2
echo "Running command..." >&2

python "${MONOPOGEN_PATH}/src/Monopogen.py" somatic \
    -a  "${MONOPOGEN_PATH}/apps" \
    -r "$CHR_FILE" \
    -t 4 \
    -i  "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/lowseq" \
    -l  /projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/lowseq_cell_data.csv \
    -s LDrefinement \
    -g   "${REF_FA}"

EXIT_CODE=$?
echo "" >&2
echo "Python script exit code: $EXIT_CODE" >&2

rm -f "$CHR_FILE"

if [ $EXIT_CODE -eq 0 ]; then
  echo "Finished $CHR successfully at $(date)" >&2
else
  echo "ERROR: $CHR failed with exit code $EXIT_CODE at $(date)" >&2
fi

echo "========================================" >&2
exit $EXIT_CODE
