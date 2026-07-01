#!/usr/bin/env python3
"""
6_enhanced_visualize.py

Enhanced CalicoST visualizations for spatial ATAC-seq CNV results.
Uses CalicoST's built-in utils_plotting functions plus custom plots.

Outputs (to Data/04_analysis/cnv/calicoST/<tissue>/plots/my_plots/):
  1. acn_heatmap_chisel.pdf      — allele-specific CN, 16-color chisel palette
  2. total_cn_heatmap.pdf        — total CN, 6-state amp/bamp/neu/bdel/del/loh
  3. rdr_baf_genome.pdf          — per-clone RDR + BAF genome profiles
  4. rdr_baf_2d_scatter.pdf      — 2D RDR vs phased AF scatter per clone
  5. spatial_clone_uncertainty.pdf — clone map with assignment certainty
  6. spatial_normal_tumor.pdf    — normal vs tumor spatial overlay
  7. spatial_tumor_purity.pdf    — estimated tumor purity per spot
  8. clone_summary.pdf           — clone sizes, per-chr CNV burden, mean CN heatmap

Usage:
    cd /projectnb/paxlab/presh/projects/spatial_atac/analysis/src/cnv_calling/calicoST
    /projectnb/paxlab/presh/env/calicost_env/bin/python3 6_enhanced_visualize.py lowseq_489
"""

import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn
from pathlib import Path

CALICOST_SRC = "/projectnb/paxlab/presh/software/CalicoST/src"
sys.path.insert(0, CALICOST_SRC)

from calicost.utils_plotting import (
    plot_acn_from_df, plot_total_cn, plot_rdr_baf_from_df,
    plot_2dscatter_rdrbaf_from_df, get_full_palette,
)
from calicost.utils_hmrf import merge_pseudobulk_by_index

PROJECT_ROOT = "/projectnb/paxlab/presh/projects/spatial_atac"
CHROM_ORDER = [str(c) for c in range(1, 23)]


def cn_to_category(a, b, neutral_total=2):
    """Map (major, minor) CN to 6-state category, relative to neutral_total."""
    total = int(a) + int(b)
    minor = min(int(a), int(b))
    if total == 0:
        return "del"
    if minor == 0 and total <= neutral_total:
        return "loh"
    if total < neutral_total:
        return "bdel" if total == neutral_total - 1 else "del"
    if total == neutral_total:
        return "neu"
    # total > neutral_total → amplification
    return "bamp" if int(a) == int(b) else "amp"


def main():
    tissue = sys.argv[1] if len(sys.argv) > 1 else "lowseq_489"

    BASE      = f"{PROJECT_ROOT}/Data/04_analysis/cnv/calicoST/{tissue}"
    # run2: results live under run2/ subfolder
    import glob as _glob
    _cna_dirs = sorted(_glob.glob(f"{BASE}/run2/cna/clone*_rectangle0_w*.0"))
    if not _cna_dirs:
        raise FileNotFoundError(
            f"No CNA output found in {BASE}/run2/cna/ — has run2 completed?\n"
            f"Expected dirs matching: {BASE}/run2/cna/clone*_rectangle0_w*.0"
        )
    CNA_DIR   = _cna_dirs[-1]   # pick most recent clone* dir
    PARSED    = f"{BASE}/parsed_inputs"
    OUT_DIR   = f"{BASE}/run2/plots/my_plots"
    Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

    # ----------------------------------------------------------------
    # Load data
    # ----------------------------------------------------------------
    seg          = pd.read_csv(f"{CNA_DIR}/cnv_seglevel.tsv", sep="\t")
    clone_labels = pd.read_csv(f"{CNA_DIR}/clone_labels.tsv", sep="\t")
    posterior    = np.load(f"{CNA_DIR}/posterior_clone_probability.npy")  # (n_spots, n_clones)
    meta         = pd.read_csv(f"{PARSED}/table_meta.csv.gz", sep="\t")

    # normal_candidate_barcodes.txt stores integer row indices (not barcode strings)
    _normal_txt = f"{CNA_DIR}/normal_candidate_barcodes.txt"
    if os.path.exists(_normal_txt):
        with open(_normal_txt) as f:
            normal_indices = set(int(line.strip()) for line in f if line.strip())
        normal_barcodes = set(clone_labels.iloc[sorted(normal_indices)]["BARCODES"].tolist())
        print(f"  Loaded {len(normal_barcodes)} normal barcodes from normal_candidate_barcodes.txt")
    else:
        # Fall back to the normalidx_file used as run2 purity input (actual barcode strings)
        _normalidx = f"{PROJECT_ROOT}/analysis/binsize_comparison/normal_barcodes/{tissue}_normal_barcodes.csv"
        if os.path.exists(_normalidx):
            with open(_normalidx) as f:
                normal_barcodes = set(line.strip() for line in f if line.strip())
            print(f"  normal_candidate_barcodes.txt not found — using normalidx_file: {len(normal_barcodes)} barcodes")
        else:
            normal_barcodes = set()
            print(f"  WARNING: no normal barcode source found — is_normal will be empty")

    bd = dict(np.load(f"{CNA_DIR}/binned_data.npz",           allow_pickle=True))
    rd = dict(np.load(f"{CNA_DIR}/rdrbaf_final_nstates7_smp.npz", allow_pickle=True))

    single_X           = bd["single_X"]            # (n_bins, 2, n_spots)
    single_base_nb_mean = bd["single_base_nb_mean"] # (n_bins, n_spots)
    single_total_bb_RD  = bd["single_total_bb_RD"]  # (n_bins, n_spots)
    assignment          = rd["new_assignment"]       # (n_spots,)

    n_spots      = single_X.shape[2]
    n_sub_clones = len(np.unique(assignment))       # 4 RDR sub-clones
    sub_clone_ids = sorted(np.unique(assignment).tolist())
    chr_col       = seg["CHR"].values

    # BAF-level clones come from cnv_seglevel.tsv columns
    # e.g. ['clone0 A', 'clone0 B', 'clone1 A', 'clone1 B']
    baf_clone_ids = sorted(set(
        int(col.split()[0].replace("clone", ""))
        for col in seg.columns if col.startswith("clone")
    ))
    n_baf_clones = len(baf_clone_ids)

    # Map RDR sub-clones → BAF clone (by offset order)
    # sub-clones are numbered sequentially across BAF clones
    sub_per_baf = n_sub_clones // n_baf_clones      # typically 2
    sub_to_baf  = {sc: sc // sub_per_baf for sc in sub_clone_ids}

    print(f"Tissue: {tissue}")
    print(f"  Spots: {n_spots}, Bins: {single_X.shape[0]}")
    print(f"  BAF clones: {n_baf_clones}  ({baf_clone_ids})")
    print(f"  RDR sub-clones: {n_sub_clones}  (sub→BAF: {sub_to_baf})")
    for sc in sub_clone_ids:
        print(f"    Sub-clone {sc}: {np.sum(assignment == sc)} spots → BAF clone {sub_to_baf[sc]}")

    tab10 = plt.get_cmap("tab10")
    clone_colors = [tab10(c) for c in range(n_sub_clones)]

    # ----------------------------------------------------------------
    # Load tetraploid seglevel — has all 4 sub-clone columns with real
    # CN variation; use for coloring RDR/BAF plots and total-CN heatmap.
    # Default diploid seglevel is all (1,1) for this dataset.
    # ----------------------------------------------------------------
    seg_tet = pd.read_csv(f"{CNA_DIR}/cnv_tetraploid_seglevel.tsv", sep="\t")
    tet_clone_ids = sorted(set(
        int(col.split()[0].replace("clone", ""))
        for col in seg_tet.columns if col.startswith("clone")
    ))
    # Neutral total for tetraploid ploidy = 4 (i.e. (2,2))
    TET_NEUTRAL = 4

    # ----------------------------------------------------------------
    # Build per-sub-clone pseudobulk RDR/BAF dataframe
    # Color by tetraploid CN since diploid calls are all (1,1)
    # ----------------------------------------------------------------
    clone_index = [np.where(assignment == sc)[0] for sc in sub_clone_ids]
    X_pb, base_pb, tot_bb_pb = merge_pseudobulk_by_index(
        single_X, single_base_nb_mean, single_total_bb_RD, clone_index
    )  # X_pb: (n_bins, 2, n_sub_clones)

    df_rdr = pd.DataFrame({"CHR": chr_col})
    for idx, sc in enumerate(sub_clone_ids):
        rdr = X_pb[:, 0, idx] / np.maximum(base_pb[:, idx], 1e-8)
        raw_baf = np.where(tot_bb_pb[:, idx] > 0,
                           X_pb[:, 1, idx] / tot_bb_pb[:, idx], 0.5)
        baf = np.minimum(raw_baf, 1.0 - raw_baf)   # fold to [0, 0.5]
        df_rdr[f"clone{sc} RD"]  = rdr
        df_rdr[f"clone{sc} BAF"] = baf
        # Use tetraploid A/B for color coding (sub-clone IDs match)
        if sc in tet_clone_ids:
            df_rdr[f"clone{sc} A"] = seg_tet[f"clone{sc} A"].values
            df_rdr[f"clone{sc} B"] = seg_tet[f"clone{sc} B"].values
        else:
            baf_c = sub_to_baf[sc]
            df_rdr[f"clone{sc} A"] = seg[f"clone{baf_c} A"].values
            df_rdr[f"clone{sc} B"] = seg[f"clone{baf_c} B"].values

    # ----------------------------------------------------------------
    # Build total-CN categorical dataframe for plot_total_cn
    # Use tetraploid seglevel (4 sub-clones, real variation)
    # ----------------------------------------------------------------
    df_total = pd.DataFrame({"CHR": chr_col})
    for sc in tet_clone_ids:
        df_total[f"clone {sc}"] = [
            cn_to_category(seg_tet[f"clone{sc} A"].iloc[i],
                           seg_tet[f"clone{sc} B"].iloc[i],
                           neutral_total=TET_NEUTRAL)
            for i in range(len(seg_tet))
        ]

    n_clones = n_sub_clones   # alias for downstream spatial/summary sections

    # ----------------------------------------------------------------
    # Merge spatial metadata with clone labels
    # ----------------------------------------------------------------
    merged = meta.merge(clone_labels, on="BARCODES")
    merged["is_normal"] = merged["BARCODES"].isin(normal_barcodes)

    certainty   = posterior.max(axis=1)                    # (n_spots,)
    spot_order  = clone_labels["BARCODES"].tolist()        # matches assignment order
    meta_indexed = meta.set_index("BARCODES")

    # Map spot index → row in merged (to align posterior with spatial coords)
    bc_to_postidx = {bc: i for i, bc in enumerate(spot_order)}
    merged["post_idx"]   = merged["BARCODES"].map(bc_to_postidx)
    merged["certainty"]  = merged["post_idx"].map(lambda i: certainty[i] if pd.notna(i) else np.nan)

    # ================================================================
    # Plot 1: ACN heatmap (chisel palette)
    # Use tetraploid if diploid has only one unique state (all 1,1)
    # ================================================================
    print("Plotting 1/8: acn_heatmap_chisel.pdf ...")
    unique_pairs = set()
    for baf_c in baf_clone_ids:
        unique_pairs.update(zip(seg[f"clone{baf_c} A"], seg[f"clone{baf_c} B"]))

    if len(unique_pairs) < 2:
        seg_tet = pd.read_csv(f"{CNA_DIR}/cnv_tetraploid_seglevel.tsv", sep="\t")
        unique_pairs_tet = set()
        for baf_c in baf_clone_ids:
            unique_pairs_tet.update(zip(seg_tet[f"clone{baf_c} A"], seg_tet[f"clone{baf_c} B"]))
        seg_acn = seg_tet
        acn_label = "Allele-Specific Copy Number (tetraploid ploidy assumption)"
        print(f"  All bins are diploid (1,1) under default ploidy — using tetraploid. Unique states: {sorted(unique_pairs_tet)}")
    else:
        seg_acn = seg
        acn_label = "Allele-Specific Copy Number"

    fig, ax = plt.subplots(1, 1, figsize=(22, 2.5 + n_baf_clones * 0.8))
    try:
        plot_acn_from_df(seg_acn, ax, add_legend=True, add_arrow=True)
        ax.set_title(f"{tissue}: {acn_label}", fontsize=13, pad=18)
        fig.tight_layout()
        fig.savefig(f"{OUT_DIR}/acn_heatmap_chisel.pdf", bbox_inches="tight")
        print("  Saved: acn_heatmap_chisel.pdf")
    except Exception as e:
        print(f"  WARNING: plot_acn_from_df failed ({e}) — saving note")
        ax.text(0.5, 0.5, f"All bins called as diploid (1,1)\nNo CNA detected under default ploidy",
                ha="center", va="center", transform=ax.transAxes, fontsize=14)
        ax.set_title(f"{tissue}: {acn_label}", fontsize=13)
        fig.savefig(f"{OUT_DIR}/acn_heatmap_chisel.pdf", bbox_inches="tight")
    plt.close(fig)

    # ================================================================
    # Plot 2: Total CN heatmap (6-state, tetraploid ploidy)
    # ================================================================
    print("Plotting 2/8: total_cn_heatmap.pdf ...")
    fig, ax = plt.subplots(1, 1, figsize=(22, 2.5 + n_sub_clones * 0.8))
    plot_total_cn(df_total, ax, palette_mode=6, add_legend=True, legend_position="lower center")
    ax.set_title(
        f"{tissue}: Total Copy Number (tetraploid ploidy, neutral = 4 copies)\n"
        f"neu=(2,2), bamp=(3,3)+, amp=unbalanced high, bdel=3 total, del≤2 total",
        fontsize=11, pad=18,
    )
    fig.tight_layout()
    fig.savefig(f"{OUT_DIR}/total_cn_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)

    # ================================================================
    # Plot 2b: HMM state assignments per bin per sub-clone
    # (7 HMM states correspond to different RDR/BAF combinations)
    # ================================================================
    print("Plotting 2b: hmm_state_heatmap.pdf ...")
    pred_cnv = rd["pred_cnv"]   # (n_bins, n_sub_clones)
    n_states = int(pred_cnv.max()) + 1

    fig, axes = plt.subplots(n_sub_clones, 1, figsize=(22, 1.5 * n_sub_clones + 1), sharex=True)
    if n_sub_clones == 1:
        axes = [axes]
    state_cmap = plt.get_cmap("tab10", n_states)
    for idx, sc in enumerate(sub_clone_ids):
        n_in_clone = np.sum(assignment == sc)
        im = axes[idx].imshow(
            pred_cnv[:, idx].reshape(1, -1),
            aspect="auto", cmap=state_cmap, vmin=0, vmax=n_states - 1,
            interpolation="none",
        )
        axes[idx].set_yticks([0])
        axes[idx].set_yticklabels([f"Sub-clone {sc}\n(n={n_in_clone})"], fontsize=9)
        axes[idx].set_xticks([])

    # Chromosome boundaries
    for chrom in np.unique(chr_col):
        idx_last = np.where(chr_col == chrom)[0][-1]
        for ax in axes:
            ax.axvline(idx_last, color="white", linewidth=0.5)

    # Chromosome labels on last panel
    for chrom in np.unique(chr_col):
        idxs = np.where(chr_col == chrom)[0]
        mid = int(np.median(idxs))
        axes[-1].text(mid, 1.05, str(chrom), ha="center", va="bottom",
                      fontsize=7, transform=axes[-1].get_xaxis_transform())

    plt.colorbar(im, ax=axes, label=f"HMM state (0–{n_states-1})", shrink=0.5, pad=0.01)
    axes[0].set_title(
        f"{tissue}: HMM state assignment per bin per sub-clone\n"
        f"(states reflect RDR+BAF combinations; integer CN calls may be diploid if signal is weak)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(f"{OUT_DIR}/hmm_state_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)

    # ================================================================
    # Plot 3: RDR + BAF genome profiles per clone
    # ================================================================
    print("Plotting 3/8: rdr_baf_genome.pdf ...")
    fig, axes = plot_rdr_baf_from_df(
        df_rdr, rdr_ylim=3, baf_ylim=0.55,
        baf_yticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5],
        pointsize=5, linewidth=0, add_legend=True,
    )
    fig.suptitle(f"{tissue}: Read Depth Ratio & Phased B-Allele Frequency per Clone",
                 y=1.01, fontsize=12)
    fig.savefig(f"{OUT_DIR}/rdr_baf_genome.pdf", bbox_inches="tight")
    plt.close(fig)

    # ================================================================
    # Plot 4: 2D RDR vs phased AF scatter (one panel per sub-clone)
    # ================================================================
    print("Plotting 4/8: rdr_baf_2d_scatter.pdf ...")
    fig, axes = plt.subplots(1, n_sub_clones,
                              figsize=(3.8 * n_sub_clones, 4),
                              dpi=150, facecolor="white")
    if n_sub_clones == 1:
        axes = [axes]
    for idx, sc in enumerate(sub_clone_ids):
        add_leg = (idx == n_sub_clones - 1)
        plot_2dscatter_rdrbaf_from_df(
            df_rdr, axes[idx], cid=str(sc),
            baf_xlim=0.55, rdr_ylim=3, pointsize=12,
            linewidth=0, add_legend=add_leg,
        )
        axes[idx].set_title(f"Sub-clone {sc}\n(n={np.sum(assignment == sc)})", fontsize=10)
    fig.suptitle(f"{tissue}: RDR vs Phased AF per Sub-Clone", y=1.03, fontsize=12)
    fig.tight_layout()
    fig.savefig(f"{OUT_DIR}/rdr_baf_2d_scatter.pdf", bbox_inches="tight")
    plt.close(fig)

    # ================================================================
    # Plot 5: Spatial clone map with assignment certainty
    # ================================================================
    print("Plotting 5/8: spatial_clone_uncertainty.pdf ...")
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))

    for c in range(n_clones):
        mask = merged["clone_label"] == c
        cert_vals = merged.loc[mask, "certainty"].fillna(0.5).values
        axes[0].scatter(
            merged.loc[mask, "X"], merged.loc[mask, "Y"],
            c=[clone_colors[c]], s=3,
            alpha=cert_vals * 0.85 + 0.1,
            label=f"Clone {c} (n={mask.sum()})", rasterized=True,
        )
    axes[0].set_title("Clone assignment\n(opacity = certainty)", fontsize=11)
    axes[0].set_aspect("equal")
    axes[0].legend(markerscale=3, loc="upper right", fontsize=8)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("X"); axes[0].set_ylabel("Y")

    sc = axes[1].scatter(
        merged["X"], merged["Y"],
        c=merged["certainty"].fillna(0.5),
        s=3, cmap="RdYlGn", vmin=0, vmax=1, rasterized=True,
    )
    plt.colorbar(sc, ax=axes[1], label="Assignment certainty")
    axes[1].set_title("Clone assignment certainty\n(green = high, red = low)", fontsize=11)
    axes[1].set_aspect("equal")
    axes[1].invert_yaxis()
    axes[1].set_xlabel("X"); axes[1].set_ylabel("Y")

    fig.suptitle(f"{tissue}: Spatial Clone Map with Assignment Certainty", fontsize=13)
    fig.tight_layout()
    fig.savefig(f"{OUT_DIR}/spatial_clone_uncertainty.pdf", bbox_inches="tight")
    plt.close(fig)

    # ================================================================
    # Plot 6: Spatial normal vs tumor
    # ================================================================
    print("Plotting 6/8: spatial_normal_tumor.pdf ...")
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))

    for c in range(n_clones):
        mask = (merged["clone_label"] == c) & (~merged["is_normal"])
        axes[0].scatter(
            merged.loc[mask, "X"], merged.loc[mask, "Y"],
            c=[clone_colors[c]], s=3, alpha=0.7,
            label=f"Clone {c}", rasterized=True,
        )
    norm_mask = merged["is_normal"]
    axes[0].scatter(
        merged.loc[norm_mask, "X"], merged.loc[norm_mask, "Y"],
        c="black", s=6, marker="x",
        label=f"Normal (n={norm_mask.sum()})", rasterized=True,
    )
    axes[0].set_title("Clone map with normal spots (×)", fontsize=11)
    axes[0].set_aspect("equal")
    axes[0].legend(markerscale=2, loc="upper right", fontsize=8)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("X"); axes[0].set_ylabel("Y")

    axes[1].scatter(
        merged.loc[~norm_mask, "X"], merged.loc[~norm_mask, "Y"],
        c="tomato", s=3, alpha=0.6, label="Tumor", rasterized=True,
    )
    axes[1].scatter(
        merged.loc[norm_mask, "X"], merged.loc[norm_mask, "Y"],
        c="steelblue", s=4, label=f"Normal (n={norm_mask.sum()})", rasterized=True,
    )
    axes[1].set_title(
        f"Tumor vs Normal\n({norm_mask.sum()} normal / {len(merged)} total)", fontsize=11
    )
    axes[1].set_aspect("equal")
    axes[1].legend(markerscale=3, loc="upper right", fontsize=8)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("X"); axes[1].set_ylabel("Y")

    fig.suptitle(f"{tissue}: Tumor vs Normal Spatial Distribution", fontsize=13)
    fig.tight_layout()
    fig.savefig(f"{OUT_DIR}/spatial_normal_tumor.pdf", bbox_inches="tight")
    plt.close(fig)

    # ================================================================
    # Plot 7: Estimated tumor purity per spot
    # Approach: identify the near-normal clone (highest fraction of
    # normal_candidate_barcodes) → purity = 1 − P(normal clone).
    # ================================================================
    print("Plotting 7/8: spatial_tumor_purity.pdf ...")

    # Find the normal clone: whichever clone has the most normal candidate spots
    clone_normal_counts = np.array([
        np.sum((clone_labels["clone_label"] == c) &
               clone_labels["BARCODES"].isin(normal_barcodes))
        for c in range(n_clones)
    ])
    normal_clone = int(np.argmax(clone_normal_counts))
    print(f"  Normal-enriched clone: {normal_clone} "
          f"({clone_normal_counts[normal_clone]} / {norm_mask.sum()} normal spots)")

    # Per-spot tumor purity = 1 − posterior probability of the normal clone
    purity_est = 1.0 - posterior[:, normal_clone]   # (n_spots,)
    merged["purity"] = merged["post_idx"].map(
        lambda i: purity_est[i] if pd.notna(i) else np.nan
    )

    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))

    sc1 = axes[0].scatter(
        merged["X"], merged["Y"],
        c=merged["purity"].fillna(0),
        s=3, cmap="RdYlBu_r", vmin=0, vmax=1, rasterized=True,
    )
    plt.colorbar(sc1, ax=axes[0], label="Estimated tumor purity")
    axes[0].set_title("Estimated tumor purity\n(1 − P[normal clone])", fontsize=11)
    axes[0].set_aspect("equal")
    axes[0].invert_yaxis()
    axes[0].set_xlabel("X"); axes[0].set_ylabel("Y")

    # Histogram of purity per clone (show counts; posterior is near-binary
    # so density=True would explode to infinity for delta-like distributions)
    for c in range(n_clones):
        mask = merged["clone_label"] == c
        vals = merged.loc[mask, "purity"].dropna().values
        axes[1].hist(vals, bins=np.linspace(0, 1, 21), alpha=0.6,
                     color=clone_colors[c], label=f"Clone {c} (n={mask.sum()})")
    axes[1].axvline(0.5, color="black", linestyle="--", linewidth=1, label="0.5 threshold")
    axes[1].set_xlabel("Estimated tumor purity  (1 − P[normal clone])")
    axes[1].set_ylabel("Number of spots")
    axes[1].set_xlim([0, 1])
    axes[1].set_title("Purity distribution per clone", fontsize=11)
    axes[1].legend(fontsize=9)

    fig.suptitle(
        f"{tissue}: Estimated Tumor Purity\n"
        f"(normal clone = {normal_clone}, based on {clone_normal_counts[normal_clone]} normal spots)",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(f"{OUT_DIR}/spatial_tumor_purity.pdf", bbox_inches="tight")
    plt.close(fig)

    # ================================================================
    # Plot 8: Clone summary
    # ================================================================
    print("Plotting 8/8: clone_summary.pdf ...")
    chr_list = [c for c in CHROM_ORDER if int(c) in np.unique(chr_col)]

    fig, axes = plt.subplots(3, 1, figsize=(16, 13))

    # --- Top: clone sizes ---
    clone_sizes = [int(np.sum(assignment == c)) for c in range(n_clones)]
    bars = axes[0].bar(range(n_clones), clone_sizes,
                       color=clone_colors[:n_clones], edgecolor="white", width=0.6)
    axes[0].set_xticks(range(n_clones))
    axes[0].set_xticklabels([f"Clone {c}" for c in range(n_clones)], fontsize=11)
    axes[0].set_ylabel("Number of spots", fontsize=11)
    axes[0].set_title(f"{tissue}: Clone sizes (total {n_spots} spots)", fontsize=12)
    for bar, v in zip(bars, clone_sizes):
        axes[0].text(bar.get_x() + bar.get_width() / 2, v + 5, str(v),
                     ha="center", va="bottom", fontsize=10)

    # --- Middle: fraction of bins with CN ≠ tetraploid neutral (4) per chromosome ---
    # Use tetraploid seglevel which has all 4 sub-clone columns with real variation
    tet_chr = seg_tet["CHR"].values
    x = np.arange(len(chr_list))
    width = 0.8 / n_sub_clones
    for idx, sc in enumerate(sub_clone_ids):
        fracs = []
        for chrom in chr_list:
            mask = tet_chr == int(chrom)
            total_cn = seg_tet.loc[mask, f"clone{sc} A"].values + seg_tet.loc[mask, f"clone{sc} B"].values
            fracs.append(np.mean(total_cn != TET_NEUTRAL))
        axes[1].bar(x + idx * width - 0.4 + width / 2, fracs, width,
                    color=clone_colors[idx], label=f"Clone {sc}", alpha=0.85, edgecolor="white")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(chr_list, fontsize=8)
    axes[1].set_xlabel("Chromosome", fontsize=11)
    axes[1].set_ylabel(f"Fraction of bins with CN ≠ {TET_NEUTRAL}", fontsize=11)
    axes[1].set_title("Per-chromosome CNV burden per sub-clone (tetraploid neutral)", fontsize=12)
    axes[1].legend(fontsize=9, ncol=n_sub_clones)
    axes[1].set_ylim([0, 1.05])

    # --- Bottom: mean total CN heatmap (sub-clones × chromosomes, tetraploid) ---
    mean_cn = np.zeros((n_sub_clones, len(chr_list)))
    for idx, sc in enumerate(sub_clone_ids):
        for ci, chrom in enumerate(chr_list):
            mask = tet_chr == int(chrom)
            mean_cn[idx, ci] = np.mean(
                seg_tet.loc[mask, f"clone{sc} A"].values + seg_tet.loc[mask, f"clone{sc} B"].values
            )

    im = axes[2].imshow(mean_cn, aspect="auto", cmap="RdBu_r",
                        vmin=TET_NEUTRAL - 2, vmax=TET_NEUTRAL + 2)
    axes[2].set_xticks(range(len(chr_list)))
    axes[2].set_xticklabels(chr_list, fontsize=8)
    axes[2].set_yticks(range(n_sub_clones))
    axes[2].set_yticklabels([f"Clone {sc}" for sc in sub_clone_ids], fontsize=10)
    axes[2].set_xlabel("Chromosome", fontsize=11)
    axes[2].set_title(f"Mean total copy number per clone per chromosome  ({TET_NEUTRAL} = tetraploid neutral)", fontsize=12)
    plt.colorbar(im, ax=axes[2], label="Mean total CN", shrink=0.8)

    # Add text annotations
    for idx in range(n_sub_clones):
        for chi in range(len(chr_list)):
            axes[2].text(chi, idx, f"{mean_cn[idx, chi]:.1f}",
                         ha="center", va="center", fontsize=6, color="black")

    fig.tight_layout()
    fig.savefig(f"{OUT_DIR}/clone_summary.pdf", bbox_inches="tight")
    plt.close(fig)

    print(f"\nAll 8 plots saved to: {OUT_DIR}/")


if __name__ == "__main__":
    main()
