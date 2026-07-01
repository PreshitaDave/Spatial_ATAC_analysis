#!/bin/bash -l
#$ -P paxlab
#$ -N proj_atac_489
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=16G
#$ -l h_rt=06:00:00
#$ -j y
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/04_project_atac_489_$JOB_ID.log
#$ -m bea
#$ -M preshita@bu.edu

set -euo pipefail
module load R
Rscript "/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/scrna_seq/integration/489/04_project_atac_onto_scrna_489.R"
