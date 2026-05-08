#!/bin/bash -l
#$ -P paxlab
#$ -N allelo_deep
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l mem_per_core=8G
#$ -l h_rt=48:00:00
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/allelo_deep.$JOB_ID.out
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/allelo_deep.$JOB_ID.err
#$ -j n

set -euo pipefail

if ! command -v module >/dev/null 2>&1; then
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
fi

module load R
module load samtools

cd /projectnb/paxlab/presh/projects/spatial_atac

# Step 1: ensure deepseq BAMs expose full-cell CB tags for VarTrix.
bash analysis/src/alleloscope/deepseq/retag_deepseq_bams_cb16.sh

# Step 2: clear stale per-chromosome VarTrix outputs from failed all-zero runs.
if [[ -d Data/alleloscope/deepseq/vartrix ]]; then
  find Data/alleloscope/deepseq/vartrix -mindepth 1 -maxdepth 1 -type d -name 'chr*' -exec rm -rf {} +
fi

# Step 3: run Alleloscope prep using retagged BAM directory.
export DEEPSEQ_BAM_DIR="/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq/Bam_cb16"
R --vanilla -q -f analysis/src/alleloscope/deepseq/prepare_alleloscope_deepseq.R
