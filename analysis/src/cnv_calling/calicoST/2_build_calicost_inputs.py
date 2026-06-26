#!/usr/bin/env python3
"""
2_build_calicost_inputs.py

Convert spatial ATAC (ArchR TileMatrix export) + numbat allele counts into
the 6 CalicoST parsed_inputs checkpoint files. When these files exist,
CalicoST's run_parse_n_load() (parse_input.py:228) skips parse_visium()
entirely, allowing us to bypass the Visium/spaceranger requirement.

Phasing: raw ref/alt allele counts are written WITHOUT haplotype assignment
(B = alt allele AD, TOT = total DP). CalicoST performs its own initial phase
estimation via pseudobulk spatial partitioning internally.

Usage:
    python 2_build_calicost_inputs.py <tissue> [--snps-per-bin N]
    Example: python 2_build_calicost_inputs.py lowseq_489 --snps-per-bin 200

Outputs (written to Data/04_analysis/cnv/calicoST/<tissue>/parsed_inputs/):
    table_bininfo.csv.gz
    table_rdrbaf.csv.gz
    table_meta.csv.gz
    adjacency_mat.npz
    smooth_mat.npz
    exp_counts.pkl
"""

import sys
import os
import argparse
import logging
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger()

PROJECT_ROOT = "/projectnb/paxlab/presh/projects/spatial_atac"
CALICOST_DIR = "/projectnb/paxlab/presh/software/CalicoST"
CALICOST_SRC = f"{CALICOST_DIR}/src"

# Add CalicoST source to path for importing its utilities
sys.path.insert(0, CALICOST_SRC)

CHROM_ORDER = [str(c) for c in range(1, 23)]

# ============================================================================
# Argument parsing
# ============================================================================

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("tissue", help="Tissue name (e.g. lowseq_489)")
    p.add_argument("--snps-per-bin", type=int, default=200,
                   help="Target number of SNPs per genomic bin (default: 200)")
    p.add_argument("--min-bin-size-bp", type=int, default=200_000,
                   help="Minimum bin size in base pairs (default: 200000)")
    p.add_argument("--nu", type=float, default=1.0,
                   help="Phase switch nu parameter for CalicoST (default: 1.0)")
    p.add_argument("--logphase-shift", type=float, default=-2.0,
                   help="Log phase shift parameter for CalicoST (default: -2.0)")
    p.add_argument("--maxspots-pooling", type=int, default=7,
                   help="Max spots pooled in HMRF smooth_mat (default: 7)")
    return p.parse_args()


# ============================================================================
# Helper: load CalicoST utility functions
# ============================================================================

def load_calicost_utils():
    try:
        from calicost.utils_IO import get_position_cM_table, compute_phase_switch_probability_position
        from calicost.utils_hmrf import multislice_adjacency
        return get_position_cM_table, compute_phase_switch_probability_position, multislice_adjacency
    except ImportError as e:
        logger.error(f"Cannot import CalicoST utilities: {e}")
        logger.error(f"Ensure CALICOST_SRC={CALICOST_SRC} is correct and env is activated.")
        sys.exit(1)


# ============================================================================
# Step 1: Load ArchR tile matrix and spatial coords (from script 1 outputs)
# ============================================================================

def load_archr_exports(intermediate_dir):
    logger.info("Loading ArchR tile matrix exports...")

    tile_mat = mmread(f"{intermediate_dir}/archr_tilematrix.mtx").tocsc()  # tiles × cells
    tile_ranges = pd.read_csv(f"{intermediate_dir}/tile_ranges.csv")
    barcodes_df = pd.read_csv(f"{intermediate_dir}/barcodes.csv")
    spatial_df = pd.read_csv(f"{intermediate_dir}/spatial_coords.csv")

    barcodes = barcodes_df["barcode"].tolist()
    logger.info(f"  Tile matrix: {tile_mat.shape[0]} tiles x {tile_mat.shape[1]} cells")
    logger.info(f"  Spatial coords: {len(spatial_df)} rows, "
                f"{spatial_df['x_spatial'].notna().sum()} with spatial data")

    return tile_mat, tile_ranges, barcodes, spatial_df


# ============================================================================
# Step 2: Load and filter numbat allele counts
# ============================================================================

def load_allele_counts(allele_file, archr_bare_barcodes):
    logger.info(f"Loading numbat allele counts from: {allele_file}")
    df = pd.read_csv(allele_file, sep="\t", compression="gzip")
    logger.info(f"  Raw rows: {len(df)}, unique cells: {df['cell'].nunique()}, "
                f"unique SNPs: {df['snp_id'].nunique()}")

    # Filter to cells present in ArchR project
    archr_set = set(archr_bare_barcodes)
    df = df[df["cell"].isin(archr_set)].copy()
    logger.info(f"  After filtering to ArchR barcodes: {len(df)} rows, "
                f"{df['cell'].nunique()} cells, {df['snp_id'].nunique()} SNPs")

    # Filter to autosomes only (chr 1-22)
    df["CHROM"] = df["CHROM"].astype(str)
    df = df[df["CHROM"].isin(CHROM_ORDER)].copy()
    logger.info(f"  After autosome filter: {len(df)} rows")

    return df


# ============================================================================
# Step 3: Load reference annotation files
# ============================================================================

def load_hgtable(hgtable_file):
    logger.info(f"Loading gene annotation: {hgtable_file}")
    hgt = pd.read_csv(hgtable_file, sep="\t", header=0)
    # Columns: (index), name2, chrom, cdsStart, cdsEnd
    hgt.columns = hgt.columns.str.strip()
    # Normalize chromosome names (remove chr prefix if present)
    hgt["chrom"] = hgt["chrom"].astype(str).str.replace("^chr", "", regex=True)
    hgt = hgt[hgt["chrom"].isin(CHROM_ORDER)].copy()
    return hgt


def load_filter_lists(ig_gene_file, hla_bed_file):
    ig_genes = set()
    if os.path.exists(ig_gene_file):
        with open(ig_gene_file) as f:
            ig_genes = {line.strip() for line in f if line.strip()}

    hla_regions = []
    if os.path.exists(hla_bed_file):
        hla_df = pd.read_csv(hla_bed_file, sep="\t", header=None,
                              names=["chrom", "start", "end"])
        hla_df["chrom"] = hla_df["chrom"].astype(str).str.replace("^chr", "", regex=True)
        hla_regions = hla_df.to_dict("records")

    logger.info(f"  IG gene filter: {len(ig_genes)} genes")
    logger.info(f"  HLA region filter: {len(hla_regions)} regions")
    return ig_genes, hla_regions


def is_in_hla(chrom, start, end, hla_regions):
    for r in hla_regions:
        if str(chrom) == str(r["chrom"]) and start < r["end"] and end > r["start"]:
            return True
    return False


# ============================================================================
# Step 4: Genomic binning of SNPs
# ============================================================================

def build_snp_bins(allele_df, snps_per_bin, min_bin_size_bp):
    """
    Group SNPs into genomic bins. Returns DataFrame with columns:
    bin_id, chr, start, end, snp_ids (list)
    """
    logger.info(f"Building genomic bins (target {snps_per_bin} SNPs/bin, "
                f"min {min_bin_size_bp}bp)...")

    # Unique SNP positions
    snp_pos = allele_df[["snp_id", "CHROM", "POS"]].drop_duplicates("snp_id")
    snp_pos["CHROM"] = pd.Categorical(snp_pos["CHROM"], categories=CHROM_ORDER, ordered=True)
    snp_pos = snp_pos.sort_values(["CHROM", "POS"]).reset_index(drop=True)

    bins = []
    bin_id = 0
    for chrom in CHROM_ORDER:
        chrom_snps = snp_pos[snp_pos["CHROM"] == chrom].reset_index(drop=True)
        if len(chrom_snps) == 0:
            continue
        i = 0
        while i < len(chrom_snps):
            j = min(i + snps_per_bin, len(chrom_snps))
            chunk = chrom_snps.iloc[i:j]
            bin_start = int(chunk["POS"].min())
            bin_end   = int(chunk["POS"].max())
            # Enforce minimum bin size
            if (bin_end - bin_start) < min_bin_size_bp and j < len(chrom_snps):
                # Extend until min size is met
                while j < len(chrom_snps) and (chrom_snps.iloc[j]["POS"] - bin_start) < min_bin_size_bp:
                    j += 1
                chunk = chrom_snps.iloc[i:j]
                bin_end = int(chunk["POS"].max())
            bins.append({
                "bin_id": bin_id,
                "CHR": chrom,
                "START": bin_start,
                "END": bin_end,
                "snp_ids": chunk["snp_id"].tolist(),
                "N_SNPS": len(chunk)
            })
            bin_id += 1
            i = j

    bin_df = pd.DataFrame(bins)
    logger.info(f"  Created {len(bin_df)} genomic bins across {bin_df['CHR'].nunique()} chromosomes")
    return bin_df


# ============================================================================
# Step 5: Compute per-bin per-cell allele counts (unphased: B=alt, TOT=total)
# ============================================================================

def compute_baf_per_bin(allele_df, bin_df, barcodes):
    """
    For each (bin, cell): sum AD (alt) → B, sum DP (total) → TOT.
    Returns two dense arrays: tot_mat (n_bins × n_cells), b_mat (n_bins × n_cells).
    """
    logger.info("Computing per-bin allele counts (unphased: B=ALT, TOT=DP)...")

    # Build SNP → bin_id mapping
    snp_to_bin = {}
    for _, row in bin_df.iterrows():
        for s in row["snp_ids"]:
            snp_to_bin[s] = row["bin_id"]

    allele_df = allele_df.copy()
    allele_df["bin_id"] = allele_df["snp_id"].map(snp_to_bin)
    allele_df = allele_df.dropna(subset=["bin_id"])
    allele_df["bin_id"] = allele_df["bin_id"].astype(int)

    # Barcode index
    barcode_to_idx = {b: i for i, b in enumerate(barcodes)}
    allele_df["cell_idx"] = allele_df["cell"].map(barcode_to_idx)
    allele_df = allele_df.dropna(subset=["cell_idx"])
    allele_df["cell_idx"] = allele_df["cell_idx"].astype(int)

    n_bins  = len(bin_df)
    n_cells = len(barcodes)
    tot_mat = np.zeros((n_bins, n_cells), dtype=np.float32)
    b_mat   = np.zeros((n_bins, n_cells), dtype=np.float32)

    for (bin_id, cell_idx), grp in allele_df.groupby(["bin_id", "cell_idx"]):
        tot_mat[bin_id, cell_idx] = grp["DP"].sum()
        b_mat[bin_id, cell_idx]   = grp["AD"].sum()

    # Fraction of bins with any coverage
    covered = np.sum(np.sum(tot_mat, axis=1) > 0)
    logger.info(f"  {covered} / {n_bins} bins have at least one allele count")
    return tot_mat, b_mat


# ============================================================================
# Step 6: Aggregate ATAC tile counts to genomic bins (RDR proxy)
# ============================================================================

def compute_rdr_per_bin(tile_mat, tile_ranges, bin_df):
    """
    Sum ATAC tile counts within each genomic bin's coordinates.
    tile_mat: scipy sparse (tiles × cells), CSC format for efficient column ops.
    Returns dense array exp_mat (n_bins × n_cells).
    """
    logger.info("Aggregating ATAC tile counts per genomic bin...")

    tile_chrom = tile_ranges["chr"].astype(str).str.replace("^chr", "", regex=True).values
    tile_start = tile_ranges["start"].values
    tile_end   = tile_ranges["end"].values

    n_bins  = len(bin_df)
    n_cells = tile_mat.shape[1]
    exp_mat = np.zeros((n_bins, n_cells), dtype=np.float32)

    tile_mat_csr = tile_mat.tocsr()  # row slicing is fast in CSR

    for _, brow in bin_df.iterrows():
        bid   = int(brow["bin_id"])
        chrom = str(brow["CHR"])
        bstart = int(brow["START"])
        bend   = int(brow["END"])

        # Tiles overlapping this bin
        mask = ((tile_chrom == chrom) &
                (tile_start < bend) &
                (tile_end   > bstart))
        tile_idx = np.where(mask)[0]

        if len(tile_idx) == 0:
            continue

        exp_mat[bid, :] = np.asarray(tile_mat_csr[tile_idx, :].sum(axis=0)).flatten()

    covered = np.sum(np.sum(exp_mat, axis=1) > 0)
    logger.info(f"  {covered} / {n_bins} bins have ATAC coverage")
    return exp_mat


# ============================================================================
# Step 7: Compute LOG_PHASE_TRANSITION using CalicoST utilities
# ============================================================================

def compute_log_phase_transition(bin_df, get_position_cM_table,
                                 compute_phase_switch_probability_position,
                                 geneticmap_file, nu, logphase_shift):
    """
    Compute log phase switch probability per bin boundary.
    Uses CalicoST's own get_position_cM_table() and
    compute_phase_switch_probability_position() to match internal behavior.
    Returns array of length n_bins.
    """
    logger.info("Computing LOG_PHASE_TRANSITION from genetic map...")

    log_trans = np.full(len(bin_df), np.log(0.5))  # default: maximum phase switch

    for chrom in bin_df["CHR"].unique():
        chrom_bins = bin_df[bin_df["CHR"] == chrom].sort_values("START").reset_index()
        if len(chrom_bins) < 2:
            continue

        # Build pairs of (bin_start, bin_end) positions for CalicoST's function
        # CalicoST expects sorted_chr_pos as array of [start, end] pairs interleaved
        positions = []
        for _, br in chrom_bins.iterrows():
            positions.append((int(chrom), int(br["START"])))
            positions.append((int(chrom), int(br["END"])))

        sorted_chr_pos = np.array(positions)  # shape (2*n_bins, 2): (chr, pos)
        try:
            position_cM = get_position_cM_table(sorted_chr_pos, geneticmap_file)
            phase_switch_prob = compute_phase_switch_probability_position(
                position_cM, sorted_chr_pos, nu)
            # CalicoST takes every other element (bin boundaries, not midpoints)
            log_trans_chrom = np.minimum(
                np.log(0.5),
                np.log(phase_switch_prob[1::2]) - logphase_shift
            )
            # Map back to original bin_df indices
            orig_idx = chrom_bins["index"].values
            log_trans[orig_idx] = log_trans_chrom[:len(orig_idx)]
        except Exception as e:
            logger.warning(f"  Phase transition failed for chr{chrom}: {e}. Using log(0.5).")

    logger.info(f"  LOG_PHASE_TRANSITION range: [{log_trans.min():.3f}, {log_trans.max():.3f}]")
    return log_trans


# ============================================================================
# Step 8: Annotate bins with genes (from hgTables), filter IG/HLA
# ============================================================================

def annotate_bins_with_genes(bin_df, hgt, ig_genes, hla_regions):
    logger.info("Annotating bins with genes (filtering IG/HLA)...")

    gene_col = []
    for _, brow in bin_df.iterrows():
        chrom  = str(brow["CHR"])
        bstart = int(brow["START"])
        bend   = int(brow["END"])

        if is_in_hla(chrom, bstart, bend, hla_regions):
            gene_col.append("")
            continue

        overlap = hgt[
            (hgt["chrom"] == chrom) &
            (hgt["cdsStart"] < bend) &
            (hgt["cdsEnd"]   > bstart)
        ]
        gene_names = [g for g in overlap["name2"].unique()
                      if isinstance(g, str) and g not in ig_genes and g.strip()]
        gene_col.append(" ".join(gene_names))

    bin_df = bin_df.copy()
    bin_df["INCLUDED_GENES"] = gene_col
    bin_df["ARM"] = "."
    bin_df["INCLUDED_SNP_IDS"] = bin_df["snp_ids"].apply(lambda x: " ".join(x))
    return bin_df


# ============================================================================
# Step 9: Build adjacency and smooth matrices using CalicoST's own function
# ============================================================================

def build_adjacency(spatial_df, barcodes, tot_mat, exp_mat, multislice_adjacency_fn,
                    maxspots_pooling, construct_adjacency_w):
    """
    Build hexagonal-grid adjacency and KNN smooth matrices.
    Uses array_row/array_col for hexagonal grid construction.
    """
    import anndata

    logger.info("Building adjacency and smooth matrices...")

    # Align spatial coords to barcodes order
    # spatial_df has archr_barcode column which matches barcodes list
    spatial_indexed = spatial_df.set_index("archr_barcode")
    coords_xy = np.zeros((len(barcodes), 2), dtype=float)
    for i, bc in enumerate(barcodes):
        if bc in spatial_indexed.index:
            row = spatial_indexed.loc[bc]
            coords_xy[i, 0] = row["x_spatial"] if pd.notna(row["x_spatial"]) else 0.0
            coords_xy[i, 1] = row["y_spatial"] if pd.notna(row["y_spatial"]) else 0.0

    # Minimal AnnData for multislice_adjacency
    n_cells = len(barcodes)
    n_bins  = exp_mat.shape[0]

    # exp_counts: sparse DataFrame (n_cells × n_bins), used for KNN (not hexagon)
    # For hexagon method, this is not used in adjacency construction but is required
    # as a checkpoint file, so we create a minimal version
    exp_counts_df = pd.DataFrame.sparse.from_spmatrix(
        sp.csc_matrix(exp_mat.T.astype(np.float32)),
        index=barcodes,
        columns=[f"bin_{i}" for i in range(n_bins)]
    )

    # single_total_bb_RD: (n_bins × n_cells) — allele depth per bin per cell
    single_total_bb_RD = sp.csc_matrix(tot_mat)

    sample_ids   = np.zeros(n_cells, dtype=int)
    sample_list  = ["lowseq_489"]

    try:
        adjacency_mat, smooth_mat = multislice_adjacency_fn(
            sample_ids=sample_ids,
            sample_list=sample_list,
            coords=coords_xy,
            single_total_bb_RD=single_total_bb_RD,
            exp_counts=exp_counts_df,
            across_slice_adjacency_mat=None,
            construct_adjacency_method="hexagon",
            maxspots_pooling=maxspots_pooling,
            construct_adjacency_w=construct_adjacency_w
        )
        logger.info(f"  Adjacency matrix shape: {adjacency_mat.shape}")
        logger.info(f"  Smooth matrix shape: {smooth_mat.shape}")
    except Exception as e:
        logger.warning(f"  multislice_adjacency failed ({e}); building simple kNN adjacency...")
        from sklearn.neighbors import NearestNeighbors
        k = min(6, n_cells - 1)
        nn = NearestNeighbors(n_neighbors=k + 1, metric="euclidean")
        nn.fit(coords_xy)
        dists, indices = nn.kneighbors(coords_xy)
        # Build sparse adjacency from kNN (excluding self)
        rows, cols = [], []
        for i in range(n_cells):
            for j in indices[i, 1:]:  # skip self (index 0)
                rows.append(i)
                cols.append(j)
        adjacency_mat = sp.csr_matrix(
            (np.ones(len(rows)), (rows, cols)), shape=(n_cells, n_cells)
        )
        # Symmetric smooth matrix: row-normalize
        smooth_mat = adjacency_mat.copy().astype(float)
        row_sums = np.asarray(smooth_mat.sum(axis=1)).flatten()
        row_sums[row_sums == 0] = 1
        smooth_mat = sp.diags(1.0 / row_sums) @ smooth_mat
        smooth_mat = sp.csr_matrix(smooth_mat)

    return adjacency_mat, smooth_mat, exp_counts_df


# ============================================================================
# Main
# ============================================================================

def main():
    args = parse_args()
    tissue = args.tissue

    # Paths
    inter_dir       = f"{PROJECT_ROOT}/Data/04_analysis/cnv/calicoST/{tissue}/intermediate"
    out_dir         = f"{PROJECT_ROOT}/Data/04_analysis/cnv/calicoST/{tissue}/parsed_inputs"
    allele_file     = (f"{PROJECT_ROOT}/Data/04_analysis/cnv/numbat/inputs/"
                       f"{tissue}/pileup_phase/{tissue}_allele_counts.tsv.gz")
    geneticmap_file = f"{CALICOST_DIR}/GRCh38_resources/genetic_map_GRCh38_merged.tab.gz"
    hgtable_file    = f"{CALICOST_DIR}/GRCh38_resources/hgTables_hg38_gencode.txt"
    ig_gene_file    = f"{CALICOST_DIR}/GRCh38_resources/ig_gene_list.txt"
    hla_bed_file    = f"{CALICOST_DIR}/GRCh38_resources/HLA_regions.bed"

    Path(out_dir).mkdir(parents=True, exist_ok=True)

    logger.info(f"=== Building CalicoST parsed_inputs for: {tissue} ===")

    # Load CalicoST utility functions
    get_cM, compute_phase_sw, multislice_adj = load_calicost_utils()

    # --- Load inputs ---
    tile_mat, tile_ranges, barcodes, spatial_df = load_archr_exports(inter_dir)
    bare_barcodes = spatial_df["bare_barcode"].tolist()  # same order as barcodes

    allele_df = load_allele_counts(allele_file, bare_barcodes)
    hgt       = load_hgtable(hgtable_file)
    ig_genes, hla_regions = load_filter_lists(ig_gene_file, hla_bed_file)

    # --- Genomic bins ---
    bin_df = build_snp_bins(allele_df, args.snps_per_bin, args.min_bin_size_bp)

    # --- Allele counts per bin (unphased: B = alt ALT, TOT = DP) ---
    tot_mat, b_mat = compute_baf_per_bin(allele_df, bin_df, bare_barcodes)

    # --- ATAC RDR per bin ---
    exp_mat = compute_rdr_per_bin(tile_mat, tile_ranges, bin_df)

    # --- LOG_PHASE_TRANSITION ---
    log_trans = compute_log_phase_transition(
        bin_df, get_cM, compute_phase_sw,
        geneticmap_file, args.nu, args.logphase_shift
    )

    # --- Annotate bins ---
    bin_df = annotate_bins_with_genes(bin_df, hgt, ig_genes, hla_regions)

    # --- NORMAL_COUNT: mean ATAC count per bin (diploid baseline, no normal spots) ---
    normal_count = exp_mat.mean(axis=1)  # (n_bins,)

    # ============================================================
    # Build table_bininfo
    # ============================================================
    table_bininfo = pd.DataFrame({
        "bin_id":              bin_df["bin_id"].values,
        "CHR":                 bin_df["CHR"].values,
        "ARM":                 bin_df["ARM"].values,
        "START":               bin_df["START"].values,
        "END":                 bin_df["END"].values,
        "LOG_PHASE_TRANSITION": log_trans,
        "INCLUDED_GENES":      bin_df["INCLUDED_GENES"].values,
        "INCLUDED_SNP_IDS":    bin_df["INCLUDED_SNP_IDS"].values,
        "NORMAL_COUNT":        normal_count,
        "N_SNPS":              bin_df["N_SNPS"].values,
    })

    # ============================================================
    # Build table_rdrbaf (long format: n_bins * n_cells rows)
    # ============================================================
    logger.info("Building table_rdrbaf (long format)...")
    n_bins  = len(bin_df)
    n_cells = len(barcodes)

    # CalicoST expects column-major order ("F" order) when reshaping back:
    # table_rdrbaf["EXP"].values.reshape((n_bins, n_cells), order="F")
    # This means: all bins for cell 0, then all bins for cell 1, ...
    bc_repeat  = np.tile(barcodes, n_bins)      # cell 0 repeated n_bins times, etc... NO
    # Correct F-order: for each cell, all bins
    bc_col = np.repeat(barcodes, n_bins)         # cell 0 repeated n_bins, then cell 1 ...
    # Wait — F order reshape means: element [bin_i, cell_j] = flat[bin_i + n_bins*cell_j]
    # So flat array in F order is: bin0cell0, bin1cell0, ..., binNcell0, bin0cell1, ...
    # Which means: repeat each barcode n_bins times (wrong) OR tile barcodes n_bins times?
    # F order: data[i,j] -> flat[i + n_rows*j]. So column j has indices [j*n_rows:(j+1)*n_rows]
    # → for each cell (column j), we have n_bins rows. So we should repeat each barcode n_bins times.
    bc_col  = np.repeat(barcodes, n_bins)        # cell 0 × n_bins, cell 1 × n_bins, ...
    exp_col = exp_mat.T.flatten(order="C")       # same as exp_mat.flatten(order="F")
    tot_col = tot_mat.T.flatten(order="C")
    b_col   = b_mat.T.flatten(order="C")

    table_rdrbaf = pd.DataFrame({
        "BARCODES": bc_col,
        "EXP":      exp_col.astype(np.float32),
        "TOT":      tot_col.astype(np.float32),
        "B":        b_col.astype(np.float32),
    })

    # ============================================================
    # Build table_meta
    # ============================================================
    logger.info("Building table_meta...")
    spatial_indexed = spatial_df.set_index("archr_barcode")

    x_vals = np.array([
        spatial_indexed.loc[bc, "x_spatial"] if bc in spatial_indexed.index
        else np.nan for bc in barcodes
    ])
    y_vals = np.array([
        spatial_indexed.loc[bc, "y_spatial"] if bc in spatial_indexed.index
        else np.nan for bc in barcodes
    ])

    table_meta = pd.DataFrame({
        "BARCODES": barcodes,
        "SAMPLE":   tissue,
        "X":        x_vals,
        "Y":        y_vals,
    })

    # ============================================================
    # Build adjacency + smooth matrices + exp_counts
    # ============================================================
    adjacency_mat, smooth_mat, exp_counts_df = build_adjacency(
        spatial_df, barcodes, tot_mat, exp_mat, multislice_adj,
        args.maxspots_pooling, construct_adjacency_w=1.0
    )

    # ============================================================
    # Write checkpoint files
    # ============================================================
    logger.info(f"Writing parsed_inputs to: {out_dir}")

    table_bininfo.to_csv(f"{out_dir}/table_bininfo.csv.gz", index=False, sep="\t")
    table_rdrbaf.to_csv(f"{out_dir}/table_rdrbaf.csv.gz",   index=False, sep="\t")
    table_meta.to_csv(f"{out_dir}/table_meta.csv.gz",       index=False, sep="\t")
    sp.save_npz(f"{out_dir}/adjacency_mat.npz", sp.csr_matrix(adjacency_mat))
    sp.save_npz(f"{out_dir}/smooth_mat.npz",    sp.csr_matrix(smooth_mat))
    exp_counts_df.to_pickle(f"{out_dir}/exp_counts.pkl")

    # Also save gene_snp_info (needed by load_tables_to_matrices if called)
    # Build minimal df_gene_snp from bin_df
    df_gene_snp_rows = []
    for _, brow in bin_df.iterrows():
        for sid in brow["snp_ids"]:
            df_gene_snp_rows.append({
                "bin_id": brow["bin_id"],
                "CHR":    brow["CHR"],
                "START":  brow["START"],
                "END":    brow["END"],
                "snp_id": sid,
                "gene":   None,
                "block_id": None,
            })
    df_gene_snp = pd.DataFrame(df_gene_snp_rows)
    df_gene_snp.to_csv(f"{out_dir}/gene_snp_info.csv.gz", index=False, sep="\t")

    logger.info("=== parsed_inputs complete ===")
    logger.info(f"  table_bininfo : {len(table_bininfo)} bins")
    logger.info(f"  table_rdrbaf  : {len(table_rdrbaf)} rows ({n_bins} bins x {n_cells} cells)")
    logger.info(f"  table_meta    : {len(table_meta)} spots")
    logger.info(f"  adjacency_mat : {adjacency_mat.shape}")
    logger.info(f"  smooth_mat    : {smooth_mat.shape}")
    logger.info(f"  exp_counts    : {exp_counts_df.shape}")


if __name__ == "__main__":
    main()
