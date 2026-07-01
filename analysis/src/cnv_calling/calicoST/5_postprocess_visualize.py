#!/usr/bin/env python3
"""
5_postprocess_visualize.py

Parse CalicoST outputs and generate visualizations:
  - Spatial scatter plot colored by clone label
  - Spatial scatter plot colored by tumor purity
  - Genome-wide CNV heatmap per clone
  - Phylogenetic tree (if phylogeography output available)

Usage:
    python 5_postprocess_visualize.py <tissue> [--n-clones N] [--spatial-weight W]
    Example: python 5_postprocess_visualize.py lowseq_489 --n-clones 3
"""

import sys
import os
import argparse
import logging
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger()

PROJECT_ROOT = "/projectnb/paxlab/presh/projects/spatial_atac"

CHROM_ORDER = [str(c) for c in range(1, 23)]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("tissue", help="Tissue name (e.g. lowseq_489)")
    p.add_argument("--n-clones", type=int, default=3)
    p.add_argument("--spatial-weight", type=float, default=1.0)
    p.add_argument("--n-clones-purity", type=int, default=5,
                   help="n_clones used in purity run (to locate tumorprop file)")
    p.add_argument("--output-base", default=None,
                   help="Override output base directory (default: Data/04_analysis/cnv/calicoST/<tissue>)")
    return p.parse_args()


def find_result_dir(base_dir, n_clones, spatial_weight):
    """Find CalicoST result subdirectory (tries rectangle0 first, then any)."""
    candidate = Path(base_dir) / f"clone{n_clones}_rectangle0_w{spatial_weight:.1f}"
    if candidate.exists():
        return candidate
    # Try other random seeds
    for d in sorted(Path(base_dir).glob(f"clone{n_clones}_rectangle*_w{spatial_weight:.1f}")):
        if d.is_dir():
            return d
    return None


def load_spatial_coords(tissue):
    inter_dir = f"{PROJECT_ROOT}/Data/04_analysis/cnv/calicoST/{tissue}/intermediate"
    df = pd.read_csv(f"{inter_dir}/spatial_coords.csv")
    return df.set_index("archr_barcode")


def plot_spatial_clones(clone_df, spatial_df, out_path, title):
    clones = sorted(clone_df["clone_label"].unique())
    n = len(clones)
    colors = plt.cm.get_cmap("tab10", max(n, 3)).colors

    merged = clone_df.join(spatial_df[["x_spatial", "y_spatial"]], how="left")
    merged = merged.dropna(subset=["x_spatial", "y_spatial"])

    fig, ax = plt.subplots(figsize=(6, 6))
    for i, cl in enumerate(clones):
        sub = merged[merged["clone_label"] == cl]
        ax.scatter(sub["x_spatial"], sub["y_spatial"],
                   c=[colors[i % len(colors)]], s=3, alpha=0.8, label=f"Clone {cl}")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("X spatial")
    ax.set_ylabel("Y spatial")
    ax.legend(markerscale=3, fontsize=7, loc="upper right")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    logger.info(f"  Saved: {out_path}")


def plot_spatial_purity(tumorprop_df, spatial_df, out_path, title):
    merged = tumorprop_df.join(spatial_df[["x_spatial", "y_spatial"]], how="left")
    merged = merged.dropna(subset=["x_spatial", "y_spatial", "Tumor"])

    fig, ax = plt.subplots(figsize=(6, 6))
    sc = ax.scatter(merged["x_spatial"], merged["y_spatial"],
                    c=merged["Tumor"], cmap="RdYlGn_r", s=3, alpha=0.8,
                    vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax, label="Tumor proportion")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("X spatial")
    ax.set_ylabel("Y spatial")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    logger.info(f"  Saved: {out_path}")


def plot_cnv_heatmap(cnv_df, out_path, title):
    """
    cnv_df: rows = genomic segments, columns = clone IDs, values = copy number.
    """
    if cnv_df is None or cnv_df.empty:
        logger.warning("  No CNV data to plot.")
        return

    # Sort by chromosome then start position
    if "CHR" in cnv_df.columns and "START" in cnv_df.columns:
        cnv_df["CHR_ord"] = pd.Categorical(cnv_df["CHR"].astype(str),
                                           categories=CHROM_ORDER, ordered=True)
        cnv_df = cnv_df.sort_values(["CHR_ord", "START"]).drop(columns="CHR_ord")

    clone_cols = [c for c in cnv_df.columns if c not in ("CHR", "START", "END",
                                                          "ARM", "N_SNPS", "bin_id")]
    if not clone_cols:
        logger.warning("  No clone columns found in CNV data.")
        return

    mat = cnv_df[clone_cols].values.astype(float)

    # Color map centered at 2 (diploid)
    vmin, vmax = 0, min(6, np.nanmax(mat) + 1)
    cmap = plt.cm.get_cmap("RdBu_r", int(vmax - vmin + 1))

    fig, ax = plt.subplots(figsize=(max(6, len(clone_cols) * 1.5), 8))
    im = ax.imshow(mat.T, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                   interpolation="nearest")
    plt.colorbar(im, ax=ax, label="Copy number", fraction=0.03, pad=0.02)
    ax.set_yticks(range(len(clone_cols)))
    ax.set_yticklabels([f"Clone {c}" for c in clone_cols], fontsize=8)
    ax.set_xlabel("Genomic bin")
    ax.set_title(title, fontsize=10)

    # Add chromosome boundary lines
    if "CHR" in cnv_df.columns:
        chrom_vals = cnv_df["CHR"].astype(str).values
        boundaries = np.where(np.diff(chrom_vals != np.roll(chrom_vals, 1)))[0]
        for b in boundaries:
            ax.axvline(b, color="black", linewidth=0.5, alpha=0.5)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"  Saved: {out_path}")


def main():
    args = parse_args()
    tissue = args.tissue

    _base = args.output_base or f"{PROJECT_ROOT}/Data/04_analysis/cnv/calicoST/{tissue}"
    base_dir_cna    = f"{_base}/cna"
    base_dir_purity = f"{_base}/purity"
    plots_dir       = f"{_base}/plots"
    Path(plots_dir).mkdir(parents=True, exist_ok=True)

    logger.info(f"=== Postprocessing CalicoST results for: {tissue} ===")

    # Load spatial coordinates
    try:
        spatial_df = load_spatial_coords(tissue)
    except FileNotFoundError:
        logger.error("spatial_coords.csv not found. Run 1_export_archr_data.R first.")
        sys.exit(1)

    # ----------------------------------------------------------------
    # 1. Clone labels spatial plot (from CNA run)
    # ----------------------------------------------------------------
    cna_dir = find_result_dir(base_dir_cna, args.n_clones, args.spatial_weight)
    if cna_dir:
        clone_file = cna_dir / "clone_labels.tsv"
        if clone_file.exists():
            clone_df = pd.read_csv(clone_file, sep="\t", index_col=0)
            clone_df.columns = ["clone_label"]
            clone_df.index.name = "archr_barcode"
            plot_spatial_clones(
                clone_df, spatial_df,
                out_path=f"{plots_dir}/{tissue}_spatial_clones.pdf",
                title=f"{tissue}: CalicoST clone labels (n={args.n_clones})"
            )
        else:
            logger.warning(f"clone_labels.tsv not found in {cna_dir}")
    else:
        logger.warning(f"CNA result directory not found under {base_dir_cna}")

    # ----------------------------------------------------------------
    # 2. Tumor purity spatial plot (from purity run)
    # ----------------------------------------------------------------
    purity_dir = find_result_dir(base_dir_purity, args.n_clones_purity, args.spatial_weight)
    if purity_dir:
        tp_file = purity_dir / "tumorprop_spots.tsv"
        if tp_file.exists():
            tp_df = pd.read_csv(tp_file, sep="\t", index_col=0)
            tp_df.index.name = "archr_barcode"
            plot_spatial_purity(
                tp_df, spatial_df,
                out_path=f"{plots_dir}/{tissue}_spatial_purity.pdf",
                title=f"{tissue}: Tumor proportion per spot"
            )
        else:
            logger.warning(f"tumorprop_spots.tsv not found in {purity_dir}")
    else:
        logger.warning(f"Purity result directory not found under {base_dir_purity}")

    # ----------------------------------------------------------------
    # 3. CNV heatmap per clone
    # ----------------------------------------------------------------
    if cna_dir:
        # Try multiple ploidy outputs, use whichever exists
        for ploidy in ("", "_diploid", "_triploid", "_tetraploid"):
            cnv_file = cna_dir / f"cnv{ploidy}_seglevel.tsv"
            if cnv_file.exists():
                cnv_df = pd.read_csv(cnv_file, sep="\t")
                plot_cnv_heatmap(
                    cnv_df,
                    out_path=f"{plots_dir}/{tissue}_cnv{ploidy}_heatmap.pdf",
                    title=f"{tissue}: CNV per clone{ploidy.replace('_', ' ')} (segment level)"
                )

    logger.info(f"=== Postprocessing complete. Plots in: {plots_dir} ===")


if __name__ == "__main__":
    main()
