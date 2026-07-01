#!/bin/bash -l
#$ -P paxlab
#$ -N cmp_align
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 2
#$ -l mem_per_core=8G
#$ -l h_rt=01:00:00
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/03_compare_alignments_$JOB_ID.log
#$ -m bea
#$ -M preshita@bu.edu

set -euo pipefail
module load miniconda
conda activate /projectnb/paxlab/presh/env/conda_env/mosaicfield_env
python "/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/scrna_seq/integration/stalign/03_compare_alignments.py"
