#!/bin/bash -l

#$ -P paxlab
#$ -N somatic_compare
#$ -pe omp 6
#$ -l h_rt=4:00:00
#$ -l mem_per_core=4G
#$ -o /projectnb/paxlab/presh/qsub_logs_deep/somatic_compare.out
#$ -e /projectnb/paxlab/presh/qsub_logs_deep/somatic_compare.err
#$ -j n

module load R/4.4.0

cd /projectnb/paxlab/presh/projects/spatial_atac/analysis/src

Rscript 9_somatic_snv_comparison.R
