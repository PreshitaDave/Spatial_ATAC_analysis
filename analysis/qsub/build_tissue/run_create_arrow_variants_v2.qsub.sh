#!/bin/bash
#$ -P paxlab
#$ -N create_arrow_v2
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 16
#$ -l h_rt=12:00:00
#$ -l mem_per_core=8G
#$ -j y
#$ -o analysis/qsub_logs/build_tissue/create_arrow_v2_$JOB_ID.log
# ============================================================================
# QC-corrected arrow creation (v2): restricts createArrowFiles() to the
# no_edge_effect barcode whitelist via validBarcodes, with permissive
# minTSS=2/minFrags=100/maxFrags=Inf so the real filtering happens downstream
# in build_archr_variant_project_v2.R (edge effect + high-nFrags removal).
#
# Usage: qsub run_create_arrow_variants_v2.qsub.sh <tissue> <tilesize> <binarize>
# Example: qsub run_create_arrow_variants_v2.qsub.sh deepseq_488B 5000 FALSE
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

export NSLOTS=$NSLOTS

Rscript /projectnb/paxlab/presh/projects/spatial_atac/analysis/src/build_tissue/binsize_test_archr/create_arrow_variants_v2.R \
  "$TISSUE" "$TILESIZE" "$BINARIZE"

EXIT_CODE=$?
echo "[$(date +'%F %T')] [DONE] create_arrow_variants_v2.R exited with code: $EXIT_CODE"
exit $EXIT_CODE
