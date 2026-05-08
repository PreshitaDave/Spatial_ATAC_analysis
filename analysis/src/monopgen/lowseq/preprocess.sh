#!/bin/bash -l
#$ -P paxlab            # Specify the SCC project name you want to use
#$ -l h_rt=48:00:00     # Specify the hard time limit for the job
#$ -N preprocess        # Give job a name
#$ -l mem_per_core=8G  # Request memory per core (adjust as needed)
#$ -pe omp 16           # Request 16 cores for parallel processing
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/preprocess  # Specify the path for standard output logs
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/preprocess_error  # Specify the path for error logs
#$ -M preshita@bu.edu

# Keep track of information related to the current job
echo "# -------------------------------------------------"
echo "Start date: $(date)"
echo "Job name: $JOB_NAME"
echo "Job ID: $JOB_ID  $SGE_TASK_ID"
echo "Running in directory: $PWD"
echo "# -------------------------------------------------"

module load python3  # Load the Python module (adjust version as needed)

path="/projectnb/paxlab/presh/software/Monopogen" # where Monopogen is downloaded
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${path}/apps
python ${path}/src/Monopogen.py  preProcess --help

python ${path}/src/Monopogen.py preProcess \
  -b /projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/lowseq_bam_list.txt \
  -o /projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/lowseq \
  -a ${path}/apps -t 16
