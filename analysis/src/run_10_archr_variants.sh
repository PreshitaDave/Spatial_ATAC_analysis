#!/bin/bash -l

#$ -P paxlab
#$ -N archr_variants
#$ -pe omp 6
#$ -l h_rt=8:00:00
#$ -l mem_per_core=8G
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/archr_variants.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/archr_variants.err
#$ -j n

module load R/4.4.0

cd /projectnb/paxlab/presh/projects/spatial_atac/analysis/src

Rscript 10_archr_variant_plotting.R
