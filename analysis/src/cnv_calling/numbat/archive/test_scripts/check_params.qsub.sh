#!/bin/bash
#$ -N check_params
#$ -l h_rt=00:15:00
#$ -pe omp 2
#$ -P paxlab
#$ -l mem_per_core=4G
#$ -o analysis/qsub_logs/check_params_$JOB_ID.log

set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

module load R

Rscript --vanilla << 'EOFR'
library(numbat)
cat("\n=== NUMBAT run_numbat() function signature ===\n")
print(args(run_numbat))
EOFR
