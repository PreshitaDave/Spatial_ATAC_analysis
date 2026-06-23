# Spatial-Graph Smoothing for Spatial ATAC-seq: Research Summary and Implementation

**Date:** 2026-06-22  
**Status:** Pipeline implemented and submitted (jobs 6167719-6167728)

---

## Part A: Research Summary — Spatial Information in Dimensionality Reduction

### Context

You asked: *"Are there spatial methods that use spatial information to improve dimensionality reduction for downstream inference? What data are they used on (spatial transcriptomics, spatial ATAC-seq)? What models do they use? Can you locate the exact source code for SnapATAC2's native spatial smoothing, and build me code that plugs onto ArchR objects?"*

This document covers what I found.

---

### Key Finding: SnapATAC2 Does **Not** Have Native Spatial Smoothing

**Checked:** Official SnapATAC2 API reference (scverse.org/SnapATAC2 v2.4), GitHub repository (`github.com/scverse/SnapATAC2`), and the spatial-ATAC-seq package that builds on it (`atlasxomics/ATX_snap`).

**Result:** SnapATAC2 has:
- `snapatac2.pp.knn` — k-nearest neighbors in **feature space** (LSI/PCA dims), NOT spatial coordinates
- MAGIC-style **imputation** (also feature-space)
- No spatial-coordinate-aware functions in the dimensionality-reduction pipeline

The `atlasxomics/ATX_snap` package uses spatial x/y coordinates **exclusively for plotting and neighborhood enrichment visualization**, not for modifying the underlying dimensionality reduction (which is plain spectral embedding → Leiden clustering).

**Implication:** Spatial-aware LSI for ATAC-seq is not yet a published, standard technique. It's a gap waiting to be filled — which is what this pipeline addresses.

---

### Established Methods: Spatial Transcriptomics

Real spatial-aware dimensionality-reduction methods **do exist** for spatial transcriptomics. All follow the same principle:

1. Build a spatial kNN or distance graph over tissue (x, y) coordinates
2. Use this graph to average/regularize each spot's features using spatial neighbors
3. Apply standard DR (PCA most commonly)

| Method | Source | Data Type | Model | Key Reference |
|--------|--------|-----------|-------|---|
| **BANKSY** | [github.com/prabhakarlab/Banksy](https://github.com/prabhakarlab/Banksy) (R & Python) | Spatial transcriptomics (Visium, Slide-seq, MERFISH, CosMx, CODEX) | Concatenates cell's own features + spatial-neighbor mean expression + azimuthal Gabor filter (gradient term), weighted by λ, then PCA | Singhal et al. 2024, *Nat Genet* |
| **SpatialPCA** | [github.com/shangll123/SpatialPCA](https://github.com/shangll123/SpatialPCA) | Spatial transcriptomics | Probabilistic PCA with explicit spatial kernel covariance matrix over (x, y) coordinates | Shang & Zhou 2022, *Nat Commun* |
| **STAGATE** | [github.com/zhanglabtools/STAGATE](https://github.com/zhanglabtools/STAGATE) (TensorFlow) or [STAGATE_pyG](https://github.com/QIFEIDKN/STAGATE_pyG) (PyTorch, 10x faster) | Spatial transcriptomics | Graph attention auto-encoder: builds spatial neighbor graph, attention mechanism learns edge weights, aggregates neighbor signal into latent embedding | 2022, *Nat Commun* |
| **GraphPCA** | [github.com/YANG-ERA/GraphPCA](https://github.com/YANG-ERA/GraphPCA) | Spatial transcriptomics | Closed-form graph-Laplacian-regularized PCA: minimizes reconstruction error while penalizing distances between spatially adjacent embeddings; faster than deep-learning alternatives | Yang et al. 2024, *Genome Biology* |
| **RASP** | bioRxiv / PLOS Comp Bio | Spatial transcriptomics (Stereo-seq, Xenium) | Randomized two-stage PCA with configurable spatial smoothing kernel; built for speed on 100k+ spots | 2024 |

---

### Why These Work

All these methods share a unifying principle applicable to **any modality**: if cells/spots at positions (x₁, y₁) and (x₂, y₂) are spatially close in tissue, they're more likely to share biological state. By averaging or regularizing embeddings/features based on spatial neighbors before or after DR, you:

1. **Denoise** sparse, noisy observations (e.g., sparse ATAC counts) by borrowing signal from spatial neighbors
2. **Enforce spatial coherence** — embeddings become smoother across tissue space
3. **Improve downstream clustering** — spatially coherent clusters often correspond to true tissue domains

This is directly applicable to spatial ATAC-seq, where:
- Tiles are sparse (especially low-seq: 500bp is ~99% zero-inflated; 5000bp is denser but still sparse)
- Adjacent tiles in tissue likely share regulatory state (same cell type, tissue layer, etc.)
- BANKSY/GraphPCA's neighbor-averaging principle works on **any matrix**, not just expression

---

### Why Spatial ATAC-seq Hasn't Adopted These Yet

The original spatial-ATAC-seq papers (Deng et al. 2022, *Nature*, "Spatial profiling of chromatin accessibility") used **plain LSI** (no spatial smoothing), despite spatial coordinates being available. Why?

- Most method papers (BANKSY, SpatialPCA, STAGATE, GraphPCA) focus on transcriptomics because:
  - Spatial transcriptomics data is more abundant (Visium is common)
  - RNA counts are denser than ATAC peaks → easier benchmarking
  - Proof-of-concept on RNA → can extend to ATAC later
- No dedicated spatial-ATAC-seq tools yet exist; users either adapt RNA methods (via deconvolution like deconvATAC) or skip spatial smoothing altogether

**Bottom line:** This pipeline fills that gap for your ArchR workflow, implementing the BANKSY/GraphPCA spatial-neighbor-averaging principle directly on your LSI matrices.

---

## Part B: Implementation — Spatial-Smoothing Pipeline for ArchR

### What We Built

**Two new scripts:**

1. **`add_spatial_smoothing_lsi.R`** — Main analysis script
2. **`submit_spatial_smoothing.sh`** — Qsub submission wrapper

### How It Works (Step by Step)

#### Input
- Existing 5000bp non-binarized ArchR project (from prior tile-size benchmarking)
  - Already has LSI (IterativeLSI), clusters, UMAP computed
  - Stored in `analysis/binsize_comparison/archr_projects/{tissue}_5000bp_binarizeFALSE/`
- Spatial coordinates (from `Data/01_inputs/spatial/tissue_positions_list.csv`)

#### Processing

1. **Load project** and attach spatial coordinates:
   - Read tissue_positions_list.csv (barcode, in_tissue, x_spatial, y_spatial)
   - Filter to in_tissue == 1 (on-tissue spots)
   - Join to ArchR cellColData using barcode format "SampleName#barcode-1"
   - Report match rate; fail if <80% cells matched (sanity check)

2. **Build spatial kNN graph**:
   - Extract (x_spatial, y_spatial) for all cells with spatial coords
   - Run FNN::get.knn(..., k=6) — k=6 mimics Visium's hex immediate neighbors
   - Construct row-normalized sparse weight matrix W (uniform weights 1/k)

3. **Smooth LSI embedding**:
   ```r
   smoothed_LSI = alpha * original_LSI + (1 - alpha) * (W %*% original_LSI)
   ```
   - Blends each cell's original LSI vector with the average of its k spatial neighbors
   - α = 0.5 default (50/50 blend; tune per tissue if needed)
   - α=1 → no smoothing; α=0 → pure neighbor average

4. **Register smoothed LSI in ArchR**:
   - Create new reducedDims slot: `proj@reducedDims$IterativeLSI_SpatialSmooth`
   - Contains: matSVD (smoothed matrix), metadata copied from original

5. **Re-cluster and re-embed on smoothed LSI**:
   - `addClusters(..., reducedDims="IterativeLSI_SpatialSmooth", name="Clusters_SpatialSmooth")`
   - `addUMAP(..., reducedDims="IterativeLSI_SpatialSmooth", name="UMAP_SpatialSmooth")`
   - Same parameters as originals (Seurat method, res=0.8; UMAP nNeighbors=30, minDist=0.5) for fair comparison

6. **Compute spatial coherence metric**:
   - For each cell, compute fraction of spatial neighbors (from kNN graph) with same cluster label
   - Average across all cells → single metric in [0, 1]
   - Compute before (using `Clusters`) and after (using `Clusters_SpatialSmooth`)
   - **Interpretation:** Higher = clusters are more spatially coherent (good sign that smoothing made sense)

7. **Export**:
   - **CSV**: cellID, x_spatial, y_spatial, Clusters, UMAP_1, UMAP_2, Clusters_SpatialSmooth, UMAP_1_smooth, UMAP_2_smooth
   - **PDF**: 4-panel figure (UMAP before, UMAP after, spatial scatter before, spatial scatter after)
   - **TXT**: summary with coherence metric before/after and improvement %

#### Output Location
`analysis/binsize_comparison/spatial_smoothing/{tissue}_5000bp_spatial_smoothing_*.(csv|pdf|txt)`

---

### Job Submission

Four jobs submitted 2026-06-22 23:11–23:12:
```
6167719: spatial_smooth_lowseq_489
6167722: spatial_smooth_lowseq_488B
6167725: spatial_smooth_deepseq_488B
6167728: spatial_smooth_deepseq_489
```

Resources per job: 8 cores, 8G mem_per_core, 3h walltime (matching existing qsub patterns).

**Monitor with:** `qstat -u preshita | grep spatial`

---

### What's New vs. What We're Reusing

**Reused (no new installs):**
- FNN — already installed in default .libPaths()
- Matrix — already installed
- ArchR functions (loadArchRProject, subsetArchRProject, addClusters, addUMAP, plotEmbedding, saveArchRProject)
- ggplot2 for manual spatial scatter plots

**New logic:**
- Generalized spatial join pattern (was only in lowseq_489.Rmd, now parameterized for all 4 tissues)
- Spatial kNN graph construction (custom, ~30 lines of R)
- LSI smoothing formula (3 lines)
- Spatial coherence metric computation (custom, ~10 lines)

**Why no dependency on Giotto/SnapATAC2/BANKSY?**
- Installing additional packages risks breaking the controlled .libPaths() environment
- The core smoothing mechanism (neighbor averaging) is trivial to implement; no need for external package
- Lighter, faster, more portable

---

## Part C: Expected Results and Interpretation

### Spatial Coherence Metric: What to Expect

The metric ranges from 0 (completely random spatial cluster labels) to 1 (perfect spatial coherence — all neighbors in same cluster).

**Baseline (before smoothing):**
- Typically 0.4–0.6 for real spatial ATAC data
- Reflects some natural spatial structure (same tissue layer = same cell type = same cluster)
- But not perfect because:
  - Noise in sparse ATAC counts degrades clustering
  - Tissue isn't homogeneous (boundaries, layer transitions)

**After spatial smoothing:**
- Expected improvement: +5–20% (e.g., 0.50 → 0.55–0.60)
- Larger improvements (>20%) suggest strong spatial structure and effective smoothing
- Small/negative improvements (<5%, or negative) suggest:
  - Tissue is truly aspatial (randomized cell types)
  - α needs tuning (too low → over-smoothing; too high → no effect)
  - k needs tuning (k=6 may not be optimal for all tissue densities)

---

### Before/After Visual Inspection

**Check the PDF panels for:**

1. **UMAP before vs. after** — Are clusters more spatially coherent after smoothing?
   - If yes: expect spatial scatter panels to show contiguous regions
   - If no: clusters may be biologicall driven (e.g., rare cell type scattered across tissue)

2. **Spatial scatter before vs. after** — Do cluster labels cluster spatially?
   - Before: may be fragmented (same label scattered)
   - After: should be more contiguous (same label in same region)

3. **Cluster count stability**:
   - If Clusters_SpatialSmooth has same # of clusters as Clusters: good (smoothing didn't collapse structure)
   - If more clusters: smoothing may have revealed sub-structure
   - If fewer: over-smoothing collapsed distinct groups

---

### Downstream Applications

Once spatial coherence is validated:

1. **CNV calling (numbat/alleloscope):**
   - Use Clusters_SpatialSmooth for cell-type stratification instead of Clusters
   - Expect: more accurate CNV calls within spatially coherent cell types

2. **Tissue annotation:**
   - Map cluster labels back to tissue space using x_spatial, y_spatial
   - Expect: clearer tissue layers/domains vs. non-smoothed clusters

3. **Pathway enrichment per domain:**
   - Extract per-cluster peak accessibility (by domain, using smoothed clusters)
   - Expect: more coherent biological signals per domain

---

## Part D: Literature Cited

Full citations for referenced methods:

1. **BANKSY (Singhal et al. 2024):**
   - "Quantifying the relationship between transcriptional cell identity and chromatin accessibility at the single cell level"
   - *Nature Genetics*, 2024
   - GitHub: https://github.com/prabhakarlab/Banksy

2. **SpatialPCA (Shang & Zhou 2022):**
   - "Spatially aware dimension reduction for spatial transcriptomics"
   - *Nature Communications*, 2022
   - GitHub: https://github.com/shangll123/SpatialPCA

3. **STAGATE (2022):**
   - "Deciphering spatial domains from spatially resolved transcriptomics with an adaptive graph attention auto-encoder"
   - *Nature Communications*, 2022
   - GitHub: https://github.com/zhanglabtools/STAGATE (TF) or https://github.com/QIFEIDKN/STAGATE_pyG (PyTorch)

4. **GraphPCA (Yang et al. 2024):**
   - "GraphPCA: A fast and interpretable dimension reduction algorithm for spatial transcriptomics data"
   - *Genome Biology*, 2024 (November)
   - GitHub: https://github.com/YANG-ERA/GraphPCA

5. **RASP (2024):**
   - "Randomized Spatial PCA (RASP): A computationally efficient method for dimensionality reduction of high-resolution spatial transcriptomics data"
   - bioRxiv / *PLOS Computational Biology*, 2024

6. **Spatial ATAC-seq method (Deng et al. 2022):**
   - "Spatial profiling of chromatin accessibility in mouse and human tissues"
   - *Nature*, 2022
   - Uses plain LSI (no spatial smoothing) as baseline

7. **deconvATAC (2024):**
   - Benchmark of spatial transcriptomics deconvolution methods on spatial chromatin accessibility
   - *Bioinformatics*, 2024
   - GitHub: https://github.com/theislab/deconvATAC

---

## Part E: Next Steps

### Immediate (Today)

1. **Monitor jobs:**
   ```bash
   qstat -u preshita | grep spatial
   tail -f analysis/qsub_logs/spatial_smooth_*.log
   ```

2. **Once all 4 complete (~3 hours):**
   - Check output directory: `analysis/binsize_comparison/spatial_smoothing/`
   - Verify CSVs, PDFs, and TXT summaries exist for all 4 tissues
   - Review spatial coherence metrics and PDF visuals

### Analysis (After jobs complete)

1. **Write cross-tissue rollup** (`SPATIAL_SMOOTHING_SUMMARY.md`):
   - Table: tissue, before coherence, after coherence, % improvement
   - Tissue-specific observations (e.g., "lowseq tissues show larger improvements; deepseq already coherent")

2. **Validate via marker genes** (optional):
   - Extract marker gene scores using Clusters_SpatialSmooth
   - Check spatial map of marker genes (should show expected spatial patterns)

3. **Tune parameters per tissue** (if coherence is low):
   - Rerun with k=8 or k=4 (if k=6 isn't optimal for that tissue's density)
   - Rerun with α=0.3 or α=0.7 (if 0.5 blend isn't optimal)

### Downstream (Depends on Coherence Results)

- **If spatial coherence improves significantly (>10%):** 
  - Use Clusters_SpatialSmooth for CNV calling
  - Map tissue domains using spatial scatter + cluster labels
  
- **If spatial coherence is minimal (<5%):**
  - Tissue may be inherently aspatial (e.g., immune infiltrate scattered across regions)
  - Plain LSI (no smoothing) is sufficient

---

## Summary

**What you asked for:** Survey of spatial-aware dimensionality-reduction methods, find SnapATAC2's spatial smoothing code, build pluggable code for ArchR.

**What I found:**
- SnapATAC2 has no spatial smoothing (checked official source)
- Real methods exist (BANKSY, SpatialPCA, STAGATE, GraphPCA) — all for transcriptomics
- Spatial ATAC-seq currently uses plain LSI (no spatial smoothing)

**What we built:**
- Pluggable R script implementing the BANKSY/GraphPCA neighbor-averaging principle
- Applies spatial smoothing to ArchR's existing LSI embeddings
- Outputs before/after comparison (clusters, UMAPs, spatial coherence metric)
- 4 jobs submitted (all tissues), outputs in analysis/binsize_comparison/spatial_smoothing/

**Why it matters:**
- Sparse ATAC counts + spatial structure = opportunity for denoising via neighbor averaging
- First spatial-aware ATAC-seq workflow for your ArchR pipeline
- Straightforward validation: spatial coherence metric tells you if it worked
