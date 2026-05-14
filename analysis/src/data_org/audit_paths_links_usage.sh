#!/bin/bash
set -euo pipefail

ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
cd "$ROOT"
OUTDIR="analysis/documentation/audit_reports"
mkdir -p "$OUTDIR"

if command -v rg >/dev/null 2>&1; then
  SEARCH_TOOL="rg"
else
  SEARCH_TOOL="grep"
fi

SCRIPT_GLOBS=("*.R" "*.Rmd" "*.sh" "*.qsub.sh" "*.py")

if [[ "$SEARCH_TOOL" == "rg" ]]; then
  REF_SEARCH_CMD=(rg -n)
  for g in "${SCRIPT_GLOBS[@]}"; do
    REF_SEARCH_CMD+=( -g "$g" )
  done
  REF_SEARCH_CMD+=(analysis .github)
fi

# 1) Symlink audit
{
  echo -e "symlink_path\ttarget\ttarget_exists\thas_lnk_suffix"
  while IFS= read -r lnk; do
    tgt=$(readlink "$lnk" || true)
    if [[ -n "$tgt" ]] && [[ -e "$lnk" ]]; then ex=1; else ex=0; fi
    if [[ "$lnk" == *.lnk ]]; then suf=1; else suf=0; fi
    echo -e "$lnk\t$tgt\t$ex\t$suf"
  done < <(find Data analysis -type l | sort)
} > "$OUTDIR/symlink_audit.tsv"

# 2) Script inventory + likely usage
scripts_file="$OUTDIR/script_inventory.tsv"
echo -e "script_path\tref_count\tlikely_unused" > "$scripts_file"
while IFS= read -r s; do
  base=$(basename "$s")
  if [[ "$SEARCH_TOOL" == "rg" ]]; then
    cnt=$(("${REF_SEARCH_CMD[@]}" --fixed-strings "$base" 2>/dev/null || true) | wc -l | tr -d ' ')
  else
    cnt=$( (grep -R -n -F --include='*.R' --include='*.Rmd' --include='*.sh' --include='*.qsub.sh' --include='*.py' -- "$base" analysis .github 2>/dev/null || true) | wc -l | tr -d ' ')
  fi
  # one hit is usually itself; <=1 treated as likely unused/orphan for manual review
  if [[ "$cnt" -le 1 ]]; then lu=1; else lu=0; fi
  echo -e "$s\t$cnt\t$lu" >> "$scripts_file"
done < <(find analysis -type f \( -name '*.R' -o -name '*.Rmd' -o -name '*.sh' -o -name '*.qsub.sh' -o -name '*.py' \) | sort)

# 3) Absolute path checks inside scripts
abs_file="$OUTDIR/absolute_path_checks.tsv"
echo -e "script_path\tpath_literal\texists" > "$abs_file"
while IFS= read -r s; do
  (grep -Eo '/projectnb[^" )]+' "$s" 2>/dev/null || true) | sort -u | while IFS= read -r p; do
    [[ -e "$p" ]] && ex=1 || ex=0
    echo -e "$s\t$p\t$ex" >> "$abs_file"
  done
done < <(find analysis -type f \( -name '*.R' -o -name '*.Rmd' -o -name '*.sh' -o -name '*.qsub.sh' -o -name '*.py' \) | sort)

# 4) Summary
{
  echo "Audit generated at: $(date)"
  echo "Symlinks total: $(($(wc -l < "$OUTDIR/symlink_audit.tsv")-1))"
  echo "Symlinks missing .lnk: $(awk -F'\t' 'NR>1 && $4==0 {c++} END{print c+0}' "$OUTDIR/symlink_audit.tsv")"
  echo "Broken symlinks: $(awk -F'\t' 'NR>1 && $3==0 {c++} END{print c+0}' "$OUTDIR/symlink_audit.tsv")"
  echo "Scripts inventoried: $(($(wc -l < "$scripts_file")-1))"
  echo "Likely unused scripts: $(awk -F'\t' 'NR>1 && $3==1 {c++} END{print c+0}' "$scripts_file")"
  echo "Absolute paths missing: $(awk -F'\t' 'NR>1 && $3==0 {c++} END{print c+0}' "$abs_file")"
} > "$OUTDIR/summary.txt"

echo "[done] Wrote:"
echo "  $OUTDIR/symlink_audit.tsv"
echo "  $OUTDIR/script_inventory.tsv"
echo "  $OUTDIR/absolute_path_checks.tsv"
echo "  $OUTDIR/summary.txt"
