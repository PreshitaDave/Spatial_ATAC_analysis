#!/bin/bash -l
#$ -P paxlab
#$ -N compare_methods_489
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 4
#$ -l mem_per_core=8G
#$ -l h_rt=02:00:00
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/05_compare_methods_489_$JOB_ID.log
#$ -m bea
#$ -M preshita@bu.edu

set -euo pipefail
module load R
Rscript "/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/scrna_seq/integration/489/05_compare_methods_489.R"
