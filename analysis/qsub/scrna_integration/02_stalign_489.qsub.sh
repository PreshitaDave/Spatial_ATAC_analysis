#!/bin/bash -l
#$ -P paxlab
#$ -N stalign_489
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 4
#$ -l mem_per_core=16G
#$ -l h_rt=03:00:00
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/02_stalign_489_$JOB_ID.log
#$ -m bea
#$ -M preshita@bu.edu

set -euo pipefail
module load miniconda
conda activate /projectnb/paxlab/presh/env/conda_env/mosaicfield_env
python "/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/scrna_seq/integration/stalign/02_stalign_489.py"
