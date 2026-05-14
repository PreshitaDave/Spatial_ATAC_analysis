#!/usr/bin/env bash
set -euo pipefail

# Organize Data/01_inputs for spatial ATAC processing.
# - Fragments/BED files go into per-object folders.
# - Barcode caches stay in per-object barcode folders.
# - Symlinks are renamed with .lnk suffix for clarity.
# - Old archive-like clutter is consolidated under Data/01_inputs/archive.

PROJECT_ROOT="${PROJECT_ROOT:-/projectnb/paxlab/presh/projects/spatial_atac}"
INPUT_ROOT="$PROJECT_ROOT/Data/01_inputs"
FRAG_ROOT="$INPUT_ROOT/fragments"
BC_ROOT="$INPUT_ROOT/barcodes/tissue_barcodes"
BAM_ROOT="$INPUT_ROOT/bam"
ARCHIVE_ROOT="$INPUT_ROOT/archive"
TS="$(date +%Y%m%d_%H%M%S)"
RUN_ARCHIVE="$ARCHIVE_ROOT/unused_not_used_$TS"

OBJECTS=(deepseq_488B deepseq_489 lowseq_488B lowseq_489)

log() {
  echo "[$(date +'%F %T')] $*"
}

safe_mv() {
  local src="$1"
  local dst="$2"
  if [[ -e "$src" || -L "$src" ]]; then
    mkdir -p "$(dirname "$dst")"
    if [[ -e "$dst" || -L "$dst" ]]; then
      mkdir -p "$RUN_ARCHIVE/conflicts"
      log "Conflict at destination, archiving source: $src"
      mv -f "$src" "$RUN_ARCHIVE/conflicts/"
    else
      mv -f "$src" "$dst"
      log "Moved: $src -> $dst"
    fi
  fi
}

rename_symlink_with_lnk() {
  local link_path="$1"
  [[ "$(basename "$link_path")" == *.lnk ]] && return 0

  local new_path="${link_path}.lnk"
  if [[ -e "$new_path" || -L "$new_path" ]]; then
    if [[ -L "$new_path" && "$(readlink "$link_path")" == "$(readlink "$new_path")" ]]; then
      rm -f "$link_path"
      log "Removed duplicate symlink: $link_path"
      return 0
    fi
    mkdir -p "$RUN_ARCHIVE/symlink_name_conflicts"
    mv -f "$link_path" "$RUN_ARCHIVE/symlink_name_conflicts/"
    log "Archived conflicting symlink: $link_path"
    return 0
  fi

  mv -f "$link_path" "$new_path"
  log "Renamed symlink: $link_path -> $new_path"
}

mkdir -p "$FRAG_ROOT" "$BC_ROOT" "$BAM_ROOT" "$ARCHIVE_ROOT" "$RUN_ARCHIVE"
for obj in "${OBJECTS[@]}"; do
  mkdir -p "$FRAG_ROOT/$obj" "$BC_ROOT/$obj"
done

# 1) Move known tissue fragment sources from nested old folders into per-object fragment folders.
declare -A SRC_FRAG
SRC_FRAG[deepseq_488B]="$FRAG_ROOT/fragments/deepseq/tissue/deepseq_488B.fragments.tsv.gz"
SRC_FRAG[deepseq_489]="$FRAG_ROOT/fragments/deepseq/tissue/deepseq_489.fragments.tsv.gz"
SRC_FRAG[lowseq_488B]="$FRAG_ROOT/fragments/lowseq/tissue/lowseq_488B.fragments.tsv.gz"
SRC_FRAG[lowseq_489]="$FRAG_ROOT/fragments/lowseq/tissue/lowseq_489.fragments.tsv.gz"

for obj in "${OBJECTS[@]}"; do
  src="${SRC_FRAG[$obj]}"
  dst="$FRAG_ROOT/$obj/$obj.fragments.sort.filtered.bed.gz"
  if [[ -f "$src" && ! -e "$dst" ]]; then
    safe_mv "$src" "$dst"
  fi

  src_tbi="${src}.tbi"
  dst_tbi="${dst}.tbi"
  if [[ -f "$src_tbi" && ! -e "$dst_tbi" ]]; then
    safe_mv "$src_tbi" "$dst_tbi"
  fi

done

# 2) Move BED/fragment coordinate files accidentally placed in barcode folders into matching fragment object folder.
# Do not move nFrags cache tables (they belong in barcode object folders).
for obj in "${OBJECTS[@]}"; do
  if [[ -d "$BC_ROOT/$obj" ]]; then
    while IFS= read -r f; do
      safe_mv "$f" "$FRAG_ROOT/$obj/$(basename "$f")"
    done < <(find "$BC_ROOT/$obj" -maxdepth 1 -type f \( -name '*.bed' -o -name '*.bed.gz' -o -name '*.fragments.tsv.gz' -o -name '*.fragments.tsv.gz.tbi' \))
  fi

done

# 3) Ensure each object has a canonical fragment path in its own folder.
for obj in "${OBJECTS[@]}"; do
  if [[ ! -e "$FRAG_ROOT/$obj/$obj.fragments.sort.filtered.bed.gz" && -e "$FRAG_ROOT/$obj/$obj.fragments.sort.filtered.bed" ]]; then
    # Keep existing extension if source is non-gz.
    log "Canonical .bed.gz missing for $obj, found .bed file at $FRAG_ROOT/$obj/$obj.fragments.sort.filtered.bed"
  fi
done

# 4) Rename all symlinks under Data/01_inputs to include .lnk suffix.
while IFS= read -r lnk; do
  rename_symlink_with_lnk "$lnk"
done < <(find "$INPUT_ROOT" -type l)

# 5) Remove broken symlinks (now with .lnk suffix) by archiving them.
while IFS= read -r broken; do
  mkdir -p "$RUN_ARCHIVE/broken_symlinks"
  mv -f "$broken" "$RUN_ARCHIVE/broken_symlinks/"
  log "Archived broken symlink: $broken"
done < <(find "$INPUT_ROOT" -xtype l)

# 6) Consolidate old archive-like folders under Data/01_inputs/archive.
# Move archived_old_files_* from barcode tree.
for d in "$BC_ROOT"/archived_old_files_*; do
  [[ -e "$d" ]] || continue
  safe_mv "$d" "$RUN_ARCHIVE/$(basename "$d")"
done

# Move accidental analysis folder under fragments to archive.
if [[ -d "$FRAG_ROOT/analysis" ]]; then
  safe_mv "$FRAG_ROOT/analysis" "$RUN_ARCHIVE/fragments_analysis_misplaced"
fi

# Move nested old fragments staging folder to archive after migration.
if [[ -d "$FRAG_ROOT/fragments" ]]; then
  safe_mv "$FRAG_ROOT/fragments" "$RUN_ARCHIVE/fragments_old_nested_source"
fi

# 7) Summary.
log "Organization complete."
log "Run archive: $RUN_ARCHIVE"

log "Fragments per object:"
for obj in "${OBJECTS[@]}"; do
  echo "--- $FRAG_ROOT/$obj"
  ls -la "$FRAG_ROOT/$obj" || true
done

log "Barcode object folders:"
for obj in "${OBJECTS[@]}"; do
  echo "--- $BC_ROOT/$obj"
  ls -la "$BC_ROOT/$obj" || true
done

log "Symlinks currently present (all should end in .lnk):"
find "$INPUT_ROOT" -type l -printf '%p -> %l\n' | sed -n '1,300p'
