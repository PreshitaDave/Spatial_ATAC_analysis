#!/bin/bash
# Simple parse test for patched R script
# SGE options are set at submission time
echo "Started on $(hostname) at $(date)"
Rscript -e 'parse(file="/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/pipeline/somatic/9_somatic_snv_comparison.R")'
RC=$?
echo "Rscript parse exit code: $RC"
exit $RC
