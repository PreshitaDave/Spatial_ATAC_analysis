#!/bin/bash -l
#$ -P paxlab
#$ -N dbg_exact_nfr
#$ -wd /projectnb/paxlab/presh/projects/spatial_atac
#$ -pe omp 1
#$ -l h_rt=00:20:00
#$ -l mem_per_core=8G
#$ -j n

set -euo pipefail
if [[ -f /etc/profile.d/modules.sh ]]; then
  source /etc/profile.d/modules.sh
fi
module load R

Rscript - <<'RS'
project_root <- "/projectnb/paxlab/presh/projects/spatial_atac"
fragments_path <- file.path(project_root, "Data", "01_inputs", "fragments", "deepseq_488B", "deepseq_488B.fragments.sort.filtered.bed.gz")
run_shell <- function(cmd) {
  cat("[dbg3] cmd=", cmd, "\n", sep="")
  status <- system2("/bin/bash", c("-o", "pipefail", "-c", cmd))
  cat("[dbg3] status=", status, "\n", sep="")
  status
}
tmp_file <- tempfile(pattern = "dbg3_nfrags_", fileext = ".tsv")
reader <- "gzip -dc"
cmd <- sprintf(
  "%s %s | awk 'NF>=4{bc=$4; sub(/-1$/, \"\", bc); n[bc]++} END{for(b in n) print b \"\\t\" n[b]}' > %s",
  reader,
  shQuote(fragments_path),
  shQuote(tmp_file)
)
st <- run_shell(cmd)
cat("[dbg3] tmp_file=", tmp_file, "\n", sep="")
cat("[dbg3] exists=", file.exists(tmp_file), "\n", sep="")
if (file.exists(tmp_file)) {
  finfo <- file.info(tmp_file)
  cat("[dbg3] size=", finfo$size, "\n", sep="")
  if (!is.na(finfo$size) && finfo$size > 0) {
    con <- file(tmp_file, "r")
    lines <- readLines(con, n = 3)
    close(con)
    cat("[dbg3] head:\n")
    cat(paste(lines, collapse="\n"), "\n")
  }
}
quit(status = ifelse(st == 0, 0, 1))
RS
