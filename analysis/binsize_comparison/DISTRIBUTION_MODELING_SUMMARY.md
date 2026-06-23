# Distribution Modeling Summary: All Tissues

**Completed:** 2026-06-22 22:03 UTC  
**Scope:** 4 tissues × 2 tile sizes = 8 analyses  
**Result:** All completed successfully, consistent NB findings

---

## Key Findings

### Universal Negative Binomial Distribution

All tissues show **99.9% overdispersion** across both 500bp and 5000bp tile sizes, validating NB as the correct statistical model for spatial ATAC tile counts.

| Tissue | Sequencing | Informative Tiles | Overdispersed | α (dispersion) |
|--------|------------|--------------------|---------------|----------------|
| lowseq_489 | low | 79.4% | **99.9%** | 133.74 |
| lowseq_488B | low | 90.9% | **99.9%** | 127.29 |
| deepseq_489 | deep | 90.7% | **99.9%** | 72.79 |
| deepseq_488B | deep | 91.9% | **99.9%** | 59.80 |

### Sequencing Depth Effect on Dispersion

**Critical observation:** Deeper sequencing produces **lower dispersion parameters**.

- **Low-seq tissues:** α ≈ 130 (high biological/technical variance relative to mean)
- **Deep-seq tissues:** α ≈ 60-73 (lower relative variance at higher counts)

**Interpretation:**
- In low-seq (sparse), count variance is dominated by **sampling noise** and **sparse signal**
- In deep-seq (dense), variance is smaller relative to higher mean counts
- This is **expected behavior** for count data (Poisson-like for dense, NB for sparse)

**Implication for CNV calling:**
- Low-seq has **more noise** but also more **detectability** of heterozygous differences (higher variance)
- Deep-seq has **less noise** but requires higher CNV signal magnitude to stand out
- Both are valid; choice depends on **sensitivity vs. background noise** tradeoff

---

## Tile-Size Consistency

Both 500bp and 5000bp tile sizes show the **same distributional pattern** within each tissue:
- 500bp: ~79-92% informative (more sparse, more zero-inflated)
- 5000bp: ~91% informative (more counts, denser)
- Both remain **>99% overdispersed** regardless of tile size

**Recommendation:** Tile size choice should be driven by **sparsity/count-range requirements** for downstream analysis, not distributional concerns (both are NB).

---

## Statistical Validation

**Chi-squared test for Poisson null (D=1.0):**
- Observed dispersion (D) across all tiles: far exceeds 1.0
- Fraction D > 1.2 (significantly overdispersed): **99.9%**
- Fraction 0.8 ≤ D ≤ 1.2 (Poisson-consistent): **0%**

This rules out:
- ❌ Poisson distribution (too much variance)
- ❌ Binomial/multinomial (different generative process)

And validates:
- ✅ Negative Binomial (accounts for overdispersion)
- ✅ Zero-inflated NB (potential, but simple NB sufficient)

---

## Downstream Analysis Implications

### For Differential Accessibility (DESeq2-style)
- Use **NB-based test** (NOT Poisson)
- α estimates can be shared-across-regions (prior) or estimated per-region
- Recommended: Use lowseq tissue α (~130) as prior for stability on shallow samples

### For CNV Calling (numbat/alleloscope)
- NB likelihood is appropriate
- Low-seq (higher α) → higher variance detection, more sensitive to heterozygosity
- Deep-seq (lower α) → lower variance, better at finding strong imbalances
- Consider **tissue depth** when setting significance thresholds

### Model Choice
All models (Poisson, NB, ZIP, ZINB) were tested on lowseq_489 pilot:
- **Winner: Negative Binomial** (simplest adequate model, AIC/BIC competitive with zero-inflated variants)
- Zero-inflation not strong enough to justify complexity
- Standard NB sufficient for both sparse (lowseq) and dense (deepseq)

---

## Recommendations

1. **Use NB distribution** for all downstream statistical models
2. **Tile-size choice** (500bp vs 5000bp) based on analysis requirements:
   - 500bp: higher spatial resolution, sparser counts, more informative tiles needed for stable estimates
   - 5000bp: richer count distributions, more robust inference on smaller datasets
3. **Stratify by sequencing depth** in CNV analysis:
   - Low-seq: more sensitive to subtle CNVs (high noise floor)
   - Deep-seq: detect high-magnitude CNVs with precision (lower noise)
4. **Share dispersion estimates** across tissues when possible (α ≈ 100-130 for low-seq, 60-75 for deep-seq)

---

## Files Generated

- lowseq_489_distribution_modeling_report.txt
- lowseq_489_500bp_meanvar_plot.pdf
- lowseq_489_5000bp_meanvar_plot.pdf
- lowseq_489_5000bp_tile_model_summary.csv
- lowseq_488B_distribution_modeling_report.txt
- lowseq_488B_5000bp_meanvar_plot.pdf
- deepseq_488B_distribution_modeling_report.txt
- deepseq_488B_500bp_meanvar_plot.pdf
- deepseq_488B_5000bp_meanvar_plot.pdf
- deepseq_489_distribution_modeling_report.txt
- deepseq_489_500bp_meanvar_plot.pdf
- deepseq_489_5000bp_meanvar_plot.pdf

All in: `analysis/binsize_comparison/distribution_modeling/`
