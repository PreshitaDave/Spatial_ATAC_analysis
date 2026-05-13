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

export RESULT_DIR="${RESULT_DIR:-/projectnb/paxlab/presh/projects/spatial_atac/Data/05_results/variant_calling/somatic_comparison/tables}"
export PLOT_DIR="${PLOT_DIR:-/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/comparison/somatic}"

cd /projectnb/paxlab/presh/projects/spatial_atac

Rscript analysis/src/pipeline/archr/10_archr_variant_plotting_deepseq.R
