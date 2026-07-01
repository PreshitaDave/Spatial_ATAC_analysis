#!/bin/bash
set -eo pipefail

echo "=========================================="
echo "NUMBAT Reference Generation - FIXED"
echo "=========================================="
echo ""
echo "Step 1: Verify compute node allocation"
HOSTNAME=$(hostname)
if [[ $HOSTNAME == scc1* ]]; then
  echo "❌ ERROR: Still on login node ($HOSTNAME)"
  echo "❌ Run this script INSIDE a qrsh compute allocation!"
  exit 1
else
  echo "✓ Compute node: $HOSTNAME"
fi

echo ""
echo "Step 2: Running reference generation"
echo "Script: analysis/src/cnv_calling/numbat/generate_all_tissue_references.test.sh"
echo ""

bash analysis/src/cnv_calling/numbat/generate_all_tissue_references.test.sh

echo ""
echo "=========================================="
echo "✓ Reference generation complete!"
echo "=========================================="
