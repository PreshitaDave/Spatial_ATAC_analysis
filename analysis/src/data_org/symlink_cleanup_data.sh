#!/bin/bash
set -euo pipefail

# Usage:
#   bash analysis/src/data_org/symlink_cleanup_data.sh            # dry-run
#   bash analysis/src/data_org/symlink_cleanup_data.sh --apply    # apply changes
#   bash analysis/src/data_org/symlink_cleanup_data.sh --apply --materialize
#
# Behavior:
# - Scans symlinks under Data/
# - If link path does not end with .lnk, either:
#   (a) rename to <name>.lnk, or
#   (b) replace with copied file when --materialize is set and target is a file

ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
cd "$ROOT"

APPLY=0
MATERIALIZE=0
if [[ "${1:-}" == "--apply" ]]; then
  APPLY=1
  shift || true
fi
if [[ "${1:-}" == "--materialize" ]]; then
  MATERIALIZE=1
fi

REPORT_DIR="analysis/documentation/audit_reports"
mkdir -p "$REPORT_DIR"
REPORT_FILE="$REPORT_DIR/symlink_cleanup_data.$(date +%s).tsv"

echo -e "action\tpath\ttarget\tstatus\tnote" > "$REPORT_FILE"

while IFS= read -r lnk; do
  [[ "$lnk" == *.lnk ]] && continue

  target=$(readlink "$lnk" || true)
  exists=0
  [[ -e "$lnk" ]] && exists=1

  if [[ "$MATERIALIZE" -eq 1 ]]; then
    if [[ "$exists" -eq 1 && -f "$lnk" ]]; then
      if [[ "$APPLY" -eq 1 ]]; then
        tmp_copy="${lnk}.tmpcopy.$$"
        cp -a --dereference "$lnk" "$tmp_copy"
        rm -f "$lnk"
        mv "$tmp_copy" "$lnk"
      fi
      echo -e "materialize\t$lnk\t$target\t$exists\tfile_copied_or_dryrun" >> "$REPORT_FILE"
      continue
    fi
    echo -e "skip\t$lnk\t$target\t$exists\tnot_regular_file_or_missing_target" >> "$REPORT_FILE"
    continue
  fi

  new_path="${lnk}.lnk"
  if [[ -e "$new_path" || -L "$new_path" ]]; then
    echo -e "skip\t$lnk\t$target\t$exists\tdestination_exists" >> "$REPORT_FILE"
    continue
  fi

  if [[ "$APPLY" -eq 1 ]]; then
    mv "$lnk" "$new_path"
  fi
  echo -e "rename_to_lnk\t$lnk\t$target\t$exists\t${new_path}" >> "$REPORT_FILE"
done < <(find Data -type l | sort)

echo "Wrote report: $REPORT_FILE"
awk -F'\t' 'NR>1 {a[$1]++} END {for(k in a) print k"="a[k]}' "$REPORT_FILE" | sort
