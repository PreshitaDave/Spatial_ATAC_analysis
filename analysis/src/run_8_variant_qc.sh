#!/bin/bash -l
#$ -P paxlab
#$ -N variant_qc
#$ -l h_rt=24:00:00
#$ -l mem_per_core=8G
#$ -pe omp 8
#$ -o /projectnb/paxlab/presh/qsub_logs_deep/variant_qc_$JOB_ID.out
#$ -e /projectnb/paxlab/presh/qsub_logs_deep/variant_qc_$JOB_ID.err
#$ -j n

module load R/4.4.0

cd /projectnb/paxlab/presh/projects/spatial_atac/analysis/src

Rscript 8_variant_qc_comparison.R
