# UMAP Clustering Interpretation: lowseq_489 Binsize Comparison

**Generated:** 2026-06-22  
**Tissue:** lowseq_489 (pilot tissue, N=4,178 cells)  
**Variants Tested:** 500bp vs 5000bp × binarized vs not-binarized  

---

## Summary

All four tile-size/binarization combinations **consistently identify 5 major clusters (C1-C5)**, with minor differences in cell distribution and embedding geometry. Both tile sizes are viable, with trade-offs favoring **5000bp non-binarized** for CNV calling applications.

---

## Cluster Structure Stability

### Cluster Sizes Across Variants

| Metric | 500bp-no-bin | 500bp-bin | 5000bp-no-bin | 5000bp-bin |
|--------|--------------|-----------|---------------|-----------|
| **C1** | 33.0% (1377) | 33.9% (1415) | 22.5% (941) | 22.9% (956) |
| **C2** | 17.0% (712) | 16.3% (682) | 25.3% (1059) | 24.4% (1019) |
| **C3** | 10.0% (418) | 10.5% (439) | 16.2% (675) | 17.0% (712) |
| **C4** | 18.3% (766) | 18.5% (771) | 13.3% (554) | 13.4% (558) |
| **C5** | 21.7% (905) | 20.8% (871) | 22.7% (949) | 22.3% (933) |

**Key observation:** 
- **500bp variants** show skewed distribution: C1 dominates (33%), C3 underrepresented (10%)
- **5000bp variants** show more balanced distribution: all clusters in 13-27% range
- This suggests **5000bp may recover a more balanced representation** of cellular diversity

---

## Embedding Geometry

### UMAP Axis Spread (indicator of cluster separation)

| Variant | X-axis range | Y-axis range | Total spread |
|---------|--------------|--------------|--------------|
| 500bp-no-bin | 14.11 | 5.67 | 19.78 |
| 500bp-bin | 13.87 | 5.85 | 19.72 |
| 5000bp-no-bin | 13.51 | **7.79** | 21.30 |
| 5000bp-bin | 13.18 | **8.12** | 21.30 |

**Interpretation:**
- 500bp variants are **slightly "collapsed"** on Y-axis (5.67-5.85), suggesting tighter clustering
- **5000bp variants expand on Y-axis** (7.79-8.12), creating more dispersed embeddings
- This could reflect **better separation of cell types at higher resolution** (larger tiles have more signal)

---

## Binarization Effect

Comparing `binarized=TRUE` vs `FALSE` within each tile size:

### 500bp
- Cluster sizes **nearly identical** (max diff: 1.7% in C2, <2% elsewhere)
- Embedding geometry **nearly identical** (13.87 vs 14.11, 5.85 vs 5.67)
- **Conclusion:** Binarization has **minimal impact at 500bp**

### 5000bp
- Cluster sizes **nearly identical** (max diff: 1.8% in C3, <2% elsewhere)
- Embedding geometry **nearly identical** (13.18 vs 13.51, 8.12 vs 7.79)
- **Conclusion:** Binarization has **minimal impact at 5000bp**

**Overall:** Binarization (presence/absence vs. counts) is **not a significant factor** for clustering or embedding. This suggests the discretized signal from binarization is sufficient, and you can choose based on downstream analysis needs (e.g., CNV calling may benefit from counts).

---

## Recommendations for Downstream Analysis

### For CNV Calling (numbat/alleloscope)
**Recommendation: 5000bp, not binarized**
- ✅ More balanced cluster representation (better cell-type recovery)
- ✅ Higher count dynamic range (better for statistical modeling)
- ✅ Better Y-axis separation (more distinct embeddings)
- ✅ Matches sparsity comparison preference (highest density: 0.103, density/sparsity optimized)
- ✅ Larger tiles = more signal = more power for diploid CNV detection

### For Rapid Exploration
**Recommendation: Either 500bp or 5000bp**
- ✅ Binarization doesn't meaningfully change clustering
- ✅ Choose based on **computational speed** (500bp: ~6M tiles; 5000bp: ~600k tiles)
- ✅ For interactive exploration, 500bp may be faster

### For Validation Across Both Modalities
**Recommendation: Test marker genes on both tile sizes**
- The PDF `archr_umap_cluster_comparison.pdf` (Page 1) shows visual cluster separation
- Check Pages 2+ for marker gene expression consistency:
  - Do immune markers (CD3D, CD19, CD14) segregate to expected clusters?
  - Do hormone receptors (ESR1, PGR) concentrate in specific clusters?
  - Is there strong batch effects or noise?

---

## Next Steps

1. **Expand to all tissues** (underway, 18 jobs submitted):
   - lowseq_488B, deepseq_488B, deepseq_489
   - Same 4 variants each
   
2. **Marker gene validation**:
   - Visual inspection of `archr_umap_cluster_comparison.pdf` (pages 2+)
   - Confirm expected tissue-specific and lineage-specific markers

3. **CNV calling benchmarking**:
   - Once expanded tissues complete, run numbat/alleloscope with 5000bp (recommended) and 500bp
   - Compare diploid CNV sensitivity and specificity

4. **Final tile-size decision**:
   - Based on marker gene patterns + CNV calling performance → lock tile size for production

---

## Technical Notes

- All UMAPs computed with: LSI (2 iterations, 30 dims) → Seurat clustering (res=0.8) → UMAP (n_neighbors=30, min_dist=0.5)
- Marker genes extracted: CD19, MS4A1, TERT, FOXA1, SOX17, HOXD9, KLRC1, GNLY, TPSAB1, CD34, GATA1, PAX5, MME, CD14, MPO, CD3D, CD8A, ESR1, ERBB2, PGR (20 genes)
- Gene scores computed via ArchR (10kb window around gene)
