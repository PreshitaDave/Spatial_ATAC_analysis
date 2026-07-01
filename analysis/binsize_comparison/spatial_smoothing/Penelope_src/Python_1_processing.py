"""
Python_1_processing.py — D02475 PanNET spatial ATAC-seq
========================================================

Step 1 of the spatial_epigenome_FFPE pipeline structure, adapted for:
  - Human genome (hg38)
  - 10x Visium-style tissue_positions_list.csv spatial barcodes
  - PanNET (pancreatic neuroendocrine tumor) sample D02475

Replaces the BC_A/BC_B custom barcode decoding used in the FFPE mouse brain
pipeline with a direct read of the tissue_positions_list.csv whitelist.

Output: step1_imported.h5ad  (checkpoint for Python_2)
"""

import snapatac2 as snap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
from resource_logger import ResourceLogger
rl = ResourceLogger("Python_1_processing")

# ============================================================
# Configuration — edit paths here only
# ============================================================

FRAGMENT_FILE = (
    "/projectnb/rd-spat/HOME/pvarela/New_Pipeline/SnapATAC2_Signac_2/D02475_preprocess_out/chromap_output/fragments.tsv.gz"
)
SPATIAL_CSV = (
    "/projectnb/rd-spat/DATA/Lab/NET_project/Spatial_ATACseq/"
    "spatials/D2475/spatial/tissue_positions_list.csv"
)
SAMPLE_ID  = "D02475"
OUTPUT_DIR = "/projectnb/rd-spat/HOME/pvarela/New_Pipeline/SnapATAC2_Signac_2/D02475_out"

# Minimum fragments per barcode to import.
# Low floor — on-tissue spatial mask does the primary filtering.
# 200 keeps low-coverage on-tissue spots that 1000 would discard.
MIN_FRAGS = 200

# ============================================================
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# Read tissue_positions_list.csv
# Columns: barcode, in_tissue, array_row, array_col, pxl_col, pxl_row
# ============================================================

print("Reading spatial positions from tissue_positions_list.csv...")
spatial = pd.read_csv(
    SPATIAL_CSV,
    header=None,
    names=["barcode", "in_tissue", "array_row", "array_col", "pxl_col", "pxl_row"],
)
spatial = spatial[spatial["in_tissue"] == 1].copy()

spatial.index = spatial["barcode"]

in_tissue_barcodes = spatial.index.tolist()
print(f"  In-tissue spots: {len(in_tissue_barcodes)}")

# ============================================================
# Import fragments
# whitelist restricts import to in-tissue barcodes only
# ============================================================

rl.step("spatial CSV read")
print("Importing fragment file...")
data = snap.pp.import_fragments(
    FRAGMENT_FILE,
    chrom_sizes=snap.genome.hg38,
    sorted_by_barcode=False,
    whitelist=in_tissue_barcodes,
    min_num_fragments=MIN_FRAGS,
)
print(f"  Cells after import (min_frags={MIN_FRAGS}): {data.n_obs}")
rl.step("fragment import")

if data.n_obs == 0:
    raise RuntimeError(
        f"No barcodes passed min_num_fragments={MIN_FRAGS}. "
        "Check that FRAGMENT_FILE barcodes match the tissue_positions_list.csv whitelist."
    )

# ============================================================
# Attach spatial coordinates to obs
# ============================================================

print("Attaching spatial coordinates...")
aligned = spatial.reindex(data.obs_names)

data.obs["array_row"] = aligned["array_row"].values
data.obs["array_col"] = aligned["array_col"].values

# pxl_col = x (horizontal), pxl_row = y (vertical)
x_vals = aligned["pxl_col"].values.astype(str)
y_vals = aligned["pxl_row"].values.astype(str)
data.obs["x_coord"] = x_vals
data.obs["y_coord"] = y_vals

n_missing = int(np.isnan(aligned["pxl_col"].values.astype(float)).sum())
if n_missing > 0:
    print(f"  Warning: {n_missing} cells have no spatial coordinates.")

# ============================================================
# Fragment size distribution — QC plot
# ============================================================

print("Generating fragment size distribution plot...")
fig = snap.pl.frag_size_distr(data, show=False)
fig.update_layout(
    xaxis_title="Fragment size (bp)",
    yaxis_title="Count",
    width=400,
    height=450,
    plot_bgcolor="rgba(0,0,0,0)",
    paper_bgcolor="rgba(0,0,0,0)",
    xaxis=dict(
        showline=True, showgrid=False, linewidth=1, linecolor="black",
        ticks="outside", tickwidth=2, tickcolor="black", ticklen=10,
        mirror=True,
        title_font=dict(family="Arial", size=18, color="black"),
        tickfont=dict(family="Arial", size=16, color="black"),
    ),
    yaxis=dict(
        showline=True, showgrid=False, linewidth=1, linecolor="black",
        ticks="outside", tickwidth=2, tickcolor="black", ticklen=10,
        mirror=True,
        title_font=dict(family="Arial", size=18, color="black"),
        tickfont=dict(family="Arial", size=16, color="black"),
    ),
)
fig.write_image(f"{OUTPUT_DIR}/fragment_size_dist.png", scale=4)
print(f"  Saved: {OUTPUT_DIR}/fragment_size_dist.png")
rl.step("fragment size distribution plot")

# ============================================================
# Summary
# ============================================================

print(f"\nSample: {SAMPLE_ID}")
print(f"  Total in-tissue spots (CSV): {len(in_tissue_barcodes)}")
print(f"  Cells imported (>= {MIN_FRAGS} fragments): {data.n_obs}")

n_frag = list(data.obs["n_fragment"])
print(f"  Fragment count — median: {int(np.median(n_frag)):,}  "
      f"min: {int(np.min(n_frag)):,}  max: {int(np.max(n_frag)):,}")

# ============================================================
# Save checkpoint for Python_2_qc_clustering.py
# ============================================================

out_h5ad = f"{OUTPUT_DIR}/step1_imported.h5ad"
data.write(out_h5ad)
print(f"\nCheckpoint saved: {out_h5ad}")
rl.step("checkpoint write")
rl.done()
print("Run Python_2_qc_clustering.py next.")
