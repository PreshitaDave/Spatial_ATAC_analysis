#!/bin/bash -l
#$ -P  paxlab        # Set SCC project to charge
#$ -wd /projectnb/paxlab/Projects/CodeSpaghetti/tools/ATX_epigenomics-main/fastq2frags/Refseq/Refdata_scATAC_MAESTRO_GRCh38_1.1.0 # Specify the working directory
#$ -pe omp 16        # Request cores
#$ -l h_rt=24:00:00  # Specify hard time limit for the job
#$ -N index       # Name job
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

module load miniconda
conda activate /projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D01942_NG05827/chromap_env    
chromap -i -r Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -o index

echo "Job finished: $(date +%F)"

# Not running this script since it already has index