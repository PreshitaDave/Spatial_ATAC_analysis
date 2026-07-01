#!/bin/bash
#$ -P paxlab
#$ -N build_archr_v2
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 8
#$ -l h_rt=04:00:00
#$ -l mem_per_core=10G
#$ -j y
#$ -o analysis/qsub_logs/build_tissue/build_archr_v2_$JOB_ID.log
# ============================================================================
# QC-corrected ArchR variant project build (v2): filters to no_edge_effect
# barcodes, applies TSS/nFrags QC, removes high-nFrags outliers, then runs
# LSI -> Clusters -> ImputeWeights -> UMAP. Depends on the arrow file
# produced by run_create_arrow_variants_v2.qsub.sh.
#
# Usage: qsub run_build_archr_variant_project_v2.qsub.sh <tissue> <tilesize> <binarize>
# Example: qsub run_build_archr_variant_project_v2.qsub.sh deepseq_488B 5000 FALSE
# ============================================================================

set -eo pipefail

TISSUE="${1:-deepseq_488B}"
TILESIZE="${2:-5000}"
BINARIZE="${3:-FALSE}"

echo "[$(date +'%F %T')] [SETUP] Job ID: $JOB_ID"
echo "[$(date +'%F %T')] [SETUP] Hostname: $(hostname)"
echo "[$(date +'%F %T')] [SETUP] NSLOTS: $NSLOTS"
echo "[$(date +'%F %T')] [SETUP] Params: tissue=$TISSUE tilesize=$TILESIZE binarize=$BINARIZE"

set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

module load R
Rscript --version

mkdir -p /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs/build_tissue

Rscript /projectnb/paxlab/presh/projects/spatial_atac/analysis/src/build_tissue/binsize_test_archr/build_archr_variant_project_v2.R \
  "$TISSUE" "$TILESIZE" "$BINARIZE"

EXIT_CODE=$?
echo "[$(date +'%F %T')] [DONE] build_archr_variant_project_v2.R exited with code: $EXIT_CODE"
exit $EXIT_CODE
