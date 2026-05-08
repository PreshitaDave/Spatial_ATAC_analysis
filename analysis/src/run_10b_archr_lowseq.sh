#!/bin/bash -l
#$ -N archr_lowseq_var
#$ -j y
#$ -o /projectnb/paxlab/presh/qsub_logs_deep/
#$ -pe omp 2
#$ -l mem_per_core=8G
#$ -l h_rt=12:00:00
#$ -P paxlab

module load R/4.4.0
export NSLOTS=${NSLOTS:-6}

cd /projectnb/paxlab/presh/projects/spatial_atac/analysis/src
Rscript 10b_archr_variant_plotting_lowseq.R

