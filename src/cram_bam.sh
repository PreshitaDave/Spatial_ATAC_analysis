#!/bin/bash -l
#$ -P  paxlab        # Set SCC project to charge
#$ -wd /projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D1942_deepseq_alignment # Specify the working directory
#$ -pe omp 4       # Request cores
#$ -l h_rt=24:00:00  # Specify hard time limit for the job
#$ -N bam       # Name job
#$ -o /projectnb/paxlab/presh/src/Spatial_ATAC_analysis/qsub_logs/index
#$ -e /projectnb/paxlab/presh/src/Spatial_ATAC_analysis/qsub_logs/index_error
#$ -M preshita@bu.edu
#$ -j y              # Join error and output streams in one file


# Keep track of information related to the current job
echo "# -------------------------------------------------"
echo "Start date: $(date)"
echo "Job name: $JOB_NAME"
echo "Job ID: $JOB_ID  $SGE_TASK_ID"
echo "Running in directory: $PWD"
echo "# -------------------------------------------------"

module load samtools
cd /projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D1942_deepseq_alignment

samtools view -@ 4 -b -T Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -o 419272-D01942_NG06549_AT-Z0005-CATGTATCCTCTGAT.bam 419272-D01942_NG06549_AT-Z0005-CATGTATCCTCTGAT.cram

echo "Job finished: $(date +%F)"

