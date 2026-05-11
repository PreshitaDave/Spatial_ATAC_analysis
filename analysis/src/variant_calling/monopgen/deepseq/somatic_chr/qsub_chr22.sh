#!/bin/bash
#$ -P paxlab
#$ -l h_rt=48:00:00
#$ -N mono_chr22
#$ -l mem_per_core=8G
#$ -pe omp 8
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/src/variant_calling/monopgen/deepseq/somatic_chr/qsub_logs_deep/chr22.log
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/src/variant_calling/monopgen/deepseq/somatic_chr/qsub_logs_deep/chr22.error.log

mkdir -p /projectnb/paxlab/presh/projects/spatial_atac/analysis/src/variant_calling/monopgen/deepseq/somatic_chr/qsub_logs_deep
echo "Job started at $(date)" >&2

export PATH="/projectnb/paxlab/presh/env/conda_env/monopgen_env/bin:$PATH"

REGION_FILE=$(mktemp)
echo "chr22" > "$REGION_FILE"

python "/projectnb/paxlab/presh/software/Monopogen/src/Monopogen.py" somatic -a "/projectnb/paxlab/presh/software/Monopogen/apps" -r "$REGION_FILE" -t 8 -i "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq" -l "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq_cell_data.csv" -s featureInfo -g "/projectnb/paxlab/presh/projects/spatial_atac/Data/hg38_resources/Homo_sapiens_assembly38.fasta" && python "/projectnb/paxlab/presh/software/Monopogen/src/Monopogen.py" somatic -a "/projectnb/paxlab/presh/software/Monopogen/apps" -r "$REGION_FILE" -t 8 -i "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq" -l "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq_cell_data.csv" -s cellScan -g "/projectnb/paxlab/presh/projects/spatial_atac/Data/hg38_resources/Homo_sapiens_assembly38.fasta" && python "/projectnb/paxlab/presh/software/Monopogen/src/Monopogen.py" somatic -a "/projectnb/paxlab/presh/software/Monopogen/apps" -r "$REGION_FILE" -t 8 -i "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq" -l "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq_cell_data.csv" -s LDrefinement -g "/projectnb/paxlab/presh/projects/spatial_atac/Data/hg38_resources/Homo_sapiens_assembly38.fasta"

EXIT_CODE=$?
rm -f "$REGION_FILE"
exit $EXIT_CODE
