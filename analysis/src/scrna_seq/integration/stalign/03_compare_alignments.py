#!/usr/bin/env python3
"""
03_compare_alignments.py
Compare MOSAICField nonlinear alignment vs STAlign RNA-guided re-alignment.
Generates a multi-panel diagnostic figure and summary table for both tissues.

Prerequisites: 02_stalign_{488B,489}.py must have been run.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
import warnings
warnings.filterwarnings("ignore")

BASE_OUT = "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration"
PLOT_DIR = os.path.join(BASE_OUT, "comparison")
os.makedirs(PLOT_DIR, exist_ok=True)

TISSUES = ["488B", "489"]

def load_tissue(tissue):
    stalign_dir = os.path.join(BASE_OUT, "stalign", tissue)
    coords_path = os.path.join(stalign_dir, f"stalign_atac_new_coords_{tissue}.csv")
    meta_path   = os.path.join(stalign_dir, f"atac_metadata_{tissue}.csv")
    summary_path = os.path.join(stalign_dir, f"stalign_summary_{tissue}.csv")

    if not os.path.exists(coords_path):
        print(f"SKIP {tissue}: {coords_path} not found (run 02_stalign_{tissue}.py first)")
        return None, None, None

    coords = pd.read_csv(coords_path, index_col=0)
    meta   = pd.read_csv(meta_path, index_col=0) if os.path.exists(meta_path) else pd.DataFrame()
    summary = pd.read_csv(summary_path) if os.path.exists(summary_path) else pd.DataFrame()
    return coords, meta, summary

# ── Figure: side-by-side MOSAICField vs STAlign per tissue ───────────────────
fig = plt.figure(figsize=(20, 12))
gs = gridspec.GridSpec(len(TISSUES), 4, figure=fig, hspace=0.4, wspace=0.3)
all_summaries = []

for t_idx, tissue in enumerate(TISSUES):
    coords, meta, summary = load_tissue(tissue)
    if coords is None:
        continue

    all_summaries.append({"tissue": tissue,
                           **(summary.iloc[0].to_dict() if len(summary) > 0 else {})})

    disp = coords["displacement_um"].values
    x_orig = coords["x_um_orig"].values
    y_orig = coords["y_um_orig"].values
    x_new  = coords["x_um_stalign"].values
    y_new  = coords["y_um_stalign"].values
    dx = x_new - x_orig
    dy = y_new - y_orig

    # Panel 1: MOSAICField coords colored by predicted cell type
    ax1 = fig.add_subplot(gs[t_idx, 0])
    ct = meta["predicted_type"].values if "predicted_type" in meta.columns and len(meta) == len(coords) else None
    if ct is not None:
        cell_types = pd.Categorical(ct)
        cmap = plt.cm.get_cmap("Set1", len(cell_types.categories))
        for i, ctype in enumerate(cell_types.categories):
            mask = ct == ctype
            ax1.scatter(x_orig[mask], y_orig[mask], c=[cmap(i)], s=0.5, alpha=0.5, label=ctype)
        ax1.legend(markerscale=4, fontsize=5, loc="upper right")
    else:
        ax1.scatter(x_orig, y_orig, s=0.5, alpha=0.4, c="steelblue")
    ax1.set_title(f"{tissue}: MOSAICField coords\n(predicted cell types)", fontsize=9)
    ax1.set_aspect("equal"); ax1.set_xlabel("X (µm)", fontsize=7); ax1.set_ylabel("Y (µm)", fontsize=7)

    # Panel 2: STAlign-refined coords
    ax2 = fig.add_subplot(gs[t_idx, 1])
    if ct is not None:
        for i, ctype in enumerate(cell_types.categories):
            mask = ct == ctype
            ax2.scatter(x_new[mask], y_new[mask], c=[cmap(i)], s=0.5, alpha=0.5, label=ctype)
    else:
        ax2.scatter(x_new, y_new, s=0.5, alpha=0.4, c="darkorange")
    ax2.set_title(f"{tissue}: STAlign-refined coords", fontsize=9)
    ax2.set_aspect("equal"); ax2.set_xlabel("X (µm)", fontsize=7); ax2.set_ylabel("Y (µm)", fontsize=7)

    # Panel 3: displacement vectors (quiver on subsample)
    ax3 = fig.add_subplot(gs[t_idx, 2])
    n = len(x_orig)
    stride = max(1, n // 2000)
    idx = np.arange(0, n, stride)
    sc3 = ax3.scatter(x_orig[idx], y_orig[idx], c=disp[idx], s=0.8,
                      cmap="viridis", norm=Normalize(vmin=0, vmax=np.percentile(disp, 95)))
    ax3.quiver(x_orig[idx], y_orig[idx], dx[idx], dy[idx],
               angles="xy", scale_units="xy", scale=0.5,
               width=0.002, alpha=0.5, color="white")
    plt.colorbar(sc3, ax=ax3, label="displacement (µm)", fraction=0.046)
    ax3.set_title(f"{tissue}: MOSAICField→STAlign displacement\n(arrows ×2 magnified)", fontsize=9)
    ax3.set_aspect("equal"); ax3.set_xlabel("X (µm)", fontsize=7); ax3.set_ylabel("Y (µm)", fontsize=7)

    # Panel 4: displacement distribution
    ax4 = fig.add_subplot(gs[t_idx, 3])
    ax4.hist(disp, bins=80, color="steelblue", alpha=0.8, edgecolor="none")
    ax4.axvline(np.median(disp), color="red", linestyle="--", label=f"Median={np.median(disp):.1f} µm")
    ax4.axvline(np.mean(disp),   color="orange", linestyle="--", label=f"Mean={np.mean(disp):.1f} µm")
    ax4.set_xlabel("Displacement (µm)", fontsize=8)
    ax4.set_ylabel("ATAC spots", fontsize=8)
    ax4.set_title(f"{tissue}: displacement distribution", fontsize=9)
    ax4.legend(fontsize=7)

fig.suptitle("MOSAICField → STAlign alignment comparison", fontsize=13, y=1.01)
out_fig = os.path.join(PLOT_DIR, "stalign_vs_mosaicfield_comparison.pdf")
fig.savefig(out_fig, bbox_inches="tight", dpi=150)
print(f"Saved: {out_fig}")

# ── Summary table ─────────────────────────────────────────────────────────────
if all_summaries:
    summary_df = pd.DataFrame(all_summaries)
    summary_csv = os.path.join(PLOT_DIR, "stalign_alignment_summary.csv")
    summary_df.to_csv(summary_csv, index=False)
    print(f"Saved: {summary_csv}")
    print(summary_df.to_string(index=False))

print("\nAlignment comparison complete.")
