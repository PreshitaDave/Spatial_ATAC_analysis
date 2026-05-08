#!/bin/bash -l

#$ -P paxlab            # Specify the SCC project name you want to use
#$ -l h_rt=24:00:00     # Specify the hard time limit for the job
#$ -N preprocess        # Give job a name
#$ -l mem_per_core=8G  # Request memory per core (adjust as needed)
#$ -pe omp 8           # Request 8 cores for parallel processing
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs  # Specify the path for standard output logs
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs  # Specify the path for error logs

module load python3

path="/projectnb/paxlab/presh/software/Monopogen" # where Monopogen is downloaded
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${path}/apps
python ${path}/src/Monopogen.py  preProcess --help

python ${path}/src/Monopogen.py preProcess \
  -b /projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq_bam_list.txt \
  -o /projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq \
  -a ${path}/apps -t 8