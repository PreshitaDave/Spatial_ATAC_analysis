#!/bin/bash
#$ -l h_rt=00:15:00 -pe omp 2 -P paxlab -l mem_per_core=4G -o analysis/qsub_logs/test_params_$JOB_ID.log

set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

module load R
cd /projectnb/paxlab/presh/projects/spatial_atac

Rscript << 'EOFR'
library(numbat)
cat("=== run_numbat PARAMETERS ===\n")
print(formals(run_numbat))
EOFR
