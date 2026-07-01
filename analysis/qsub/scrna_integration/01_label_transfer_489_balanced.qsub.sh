#!/bin/bash -l
#$ -P paxlab
#$ -N ltx_489_bal
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=16G
#$ -l h_rt=04:00:00
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/01_label_transfer_489_balanced_$JOB_ID.log
#$ -m bea
#$ -M preshita@bu.edu

set -euo pipefail
module load R
Rscript "/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/scrna_seq/integration/489/01_label_transfer_489.R" balanced
