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

export TABLE_DIR="${TABLE_DIR:-/projectnb/paxlab/presh/projects/spatial_atac/Data/05_results/variant_calling/somatic_comparison/tables}"
export PLOT_DIR="${PLOT_DIR:-/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/comparison/somatic}"

cd /projectnb/paxlab/presh/projects/spatial_atac

Rscript analysis/src/pipeline/somatic/9_somatic_snv_comparison.R
