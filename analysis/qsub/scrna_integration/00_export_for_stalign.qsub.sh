#!/bin/bash -l
#$ -P paxlab
#$ -N exp_stalign
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 4
#$ -l mem_per_core=16G
#$ -l h_rt=01:00:00
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/00_export_for_stalign_$JOB_ID.log
#$ -m bea
#$ -M preshita@bu.edu

set -euo pipefail
module load R
Rscript "/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/scrna_seq/integration/stalign/00_export_for_stalign.R"
