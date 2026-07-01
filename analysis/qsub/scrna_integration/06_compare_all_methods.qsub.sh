#!/bin/bash -l
#$ -P paxlab
#$ -N compare_int
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 2
#$ -l mem_per_core=8G
#$ -l h_rt=01:00:00
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/06_compare_all_methods_$JOB_ID.log
#$ -m bea
#$ -M preshita@bu.edu

set -euo pipefail
module load R
Rscript "/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/scrna_seq/integration/comparison/06_compare_all_methods.R"
