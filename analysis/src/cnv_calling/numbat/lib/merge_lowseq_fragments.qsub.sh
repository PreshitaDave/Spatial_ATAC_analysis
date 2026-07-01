#!/bin/bash
# Merge lowseq fragment files (488B + 489 → combined)
# Purpose: Create combined fragment file for NUMBAT input preparation

set -eo pipefail

# CRITICAL: Initialize module system before using 'module' command
set +u
for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
  if [[ -f "$profile_file" ]]; then
    . "$profile_file" 2>/dev/null || true
    break
  fi
done
set -u

PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
FRAG_488B="${PROJECT_ROOT}/Data/01_inputs/fragments/lowseq_488B/lowseq_488B.fragments.sort.filtered.bed.gz"
FRAG_489="${PROJECT_ROOT}/Data/01_inputs/fragments/lowseq_489/lowseq_489.fragments.sort.filtered.bed.gz"
OUTPUT_DIR="${PROJECT_ROOT}/Data/01_inputs/fragments/lowseq_combined"
COMBINED_FRAG="${OUTPUT_DIR}/lowseq_combined.fragments.sort.filtered.bed.gz"
COMBINED_INDEX="${COMBINED_FRAG}.tbi"

LOG_FILE="${PROJECT_ROOT}/analysis/qsub_logs/merge_lowseq_fragments_$$.log"
exec > >(tee "${LOG_FILE}")
exec 2>&1

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting lowseq fragment merge..."

# Step 1: Load modules
echo "[STEP 1] Loading modules..."
module load samtools htslib bedtools
which tabix >/dev/null 2>&1 || { echo "[ERROR] tabix not found"; exit 1; }
which bgzip >/dev/null 2>&1 || { echo "[ERROR] bgzip not found"; exit 1; }
echo "✓ Modules loaded"

# Step 2: Verify inputs
echo "[STEP 2] Verifying input files..."
[[ ! -f "$FRAG_488B" ]] && { echo "[ERROR] File not found: $FRAG_488B"; exit 1; }
[[ ! -f "$FRAG_489" ]] && { echo "[ERROR] File not found: $FRAG_489"; exit 1; }
echo "✓ Input files verified"
echo "  488B: $(du -h $FRAG_488B | cut -f1)"
echo "  489: $(du -h $FRAG_489 | cut -f1)"

# Step 3: Create output directory
mkdir -p "$OUTPUT_DIR"

# Step 4: Merge fragments
echo "[STEP 3] Merging fragment files..."
echo "This may take 30-60 minutes..."
zcat "$FRAG_488B" "$FRAG_489" | \
  sort -k1,1V -k2,2n -k3,3n | \
  bgzip -c > "$COMBINED_FRAG"

if [[ ! -f "$COMBINED_FRAG" ]]; then
  echo "[ERROR] Merge failed - output file not created"
  exit 1
fi
echo "✓ Fragments merged: $(du -h $COMBINED_FRAG | cut -f1)"

# Step 5: Create index
echo "[STEP 4] Creating index..."
tabix -p bed "$COMBINED_FRAG"
if [[ ! -f "$COMBINED_INDEX" ]]; then
  echo "[ERROR] Index creation failed"
  exit 1
fi
echo "✓ Index created: $(du -h $COMBINED_INDEX | cut -f1)"

# Step 6: Verify combined file
echo "[STEP 5] Verifying combined file..."
COMBINED_COUNT=$(zcat "$COMBINED_FRAG" | wc -l)
echo "✓ Combined fragment count: $COMBINED_COUNT"

echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Lowseq fragment merge complete!"
echo "Output: $COMBINED_FRAG"
