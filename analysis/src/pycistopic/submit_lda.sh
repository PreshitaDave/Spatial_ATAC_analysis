#!/bin/bash -l
#$ -P  paxlab        # Set SCC project to charge
#$ -pe omp 16        # Request cores
#$ -l h_rt=24:00:00  # Specify hard time limit for the job
#$ -N pyscistopic       # Name job
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs
#$ -M preshita@bu.edu


# Load environment
# module load GCC/11.2.0
# module load OpenMPI/4.1.1
# module load R/4.1.2

# Activate conda environment
export PATH=/projectnb/paxlab/presh/env/scenicplus/bin:$PATH

# Set working directory
cd /projectnb/paxlab/presh/projects/spatial_atac/Data/pycistopic

mkdir -p outs/mallet_tmp outs/lda_models

echo "Starting LDA model training at $(date)"
python /projectnb/paxlab/presh/projects/spatial_atac/analysis/src/pycistopic/run_lda_models.py 2>&1
echo "Exit code: $?"
echo "Finished at $(date)"

