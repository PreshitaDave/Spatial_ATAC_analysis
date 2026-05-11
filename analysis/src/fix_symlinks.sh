#!/bin/bash
set -euo pipefail
git checkout -b tidy/fix-analysis-symlinks

find analysis -type l -print0 | while IFS= read -r -d '' link; do
  echo "----------------------------------------"
  echo "Link: $link"
  tgt=$(readlink -f "$link" 2>/dev/null || true)
  echo "Current target: ${tgt:-<BROKEN>}"
  read -p "Enter NEW target path (absolute or repo-relative), or ENTER to skip: " newt
  if [ -z "$newt" ]; then
    echo "Skipped."
    continue
  fi

  # normalize newt to absolute if repo-relative
  if [[ "$newt" != /* ]]; then
    newt="$(pwd)/${newt#./}"
  fi

  mkdir -p "$(dirname "$newt")"

  # move existing target into new location if it exists and user wants to preserve it
  if [ -e "$tgt" ]; then
    echo "Target exists at $tgt"
    read -p "Move existing target to $newt? [y/N]: " mvok
    if [[ "$mvok" =~ ^[Yy]$ ]]; then
      # prefer git mv if tracked
      if git ls-files --error-unmatch "$tgt" >/dev/null 2>&1; then
        git mv "$tgt" "$newt"
      else
        mv "$tgt" "$newt"
      fi
    else
      echo "Not moving existing file."
    fi
  else
    echo "No existing file at current target."
  fi

  # replace symlink with relative symlink to newt
  rm -f "$link"
  rel=$(realpath --relative-to="$(dirname "$link")" "$newt")
  ln -s "$rel" "$link"
  git add "$link"
  echo "Replaced $link -> $rel and staged change."
done

echo "Committing staged changes (if any)..."
git commit -m "chore: resolve symlinks under analysis/" || echo "Nothing to commit"
echo "Pushing branch tidy/fix-analysis-symlinks"
git push -u origin tidy/fix-analysis-symlinks