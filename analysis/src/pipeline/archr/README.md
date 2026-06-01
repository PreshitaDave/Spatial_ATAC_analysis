# archr

Purpose: scripts in this folder for the Spatial ATAC workflow.

## Pipeline Scripts (in order)

### 0_create_archr_qc_cluster.R ⭐ NEW
**Purpose**: Create ArchR objects from fragments, apply QC filtering, remove doublets, cluster, and generate comprehensive PDF reports

**Input**:
- Fragment files: `Data/01_inputs/fragments/{object}/{object}.fragments.sort.filtered.bed.gz`
- Barcode files: `Data/01_inputs/barcodes/tissue_barcodes/{object}/{object}.no_edge_effect.barcodes.tsv`

**Output**:
- ArchR objects: `Data/01_outputs/archR_objects/{object}/`
- PDF reports: `analysis/plots/cnv_analysis/archR_qc_{object}.pdf`

**Features**:
- Arrow file creation
- QC filtering (TSS, nFrags)
- Doublet removal
- LSI dimensionality reduction
- K-means clustering
- UMAP embedding
- Comprehensive PDF with all plots

**Run via**: `qsub analysis/qsub/pipeline/run_archR_qc_cluster.qsub.sh`
**Documentation**: `README_ArchR_Pipeline.md`

### Other Scripts (Legacy/Reference)
- 10_archr_variant_plotting_deepseq.R: Plotting somatic variants on ArchR objects
- 10b_archr_variant_plotting_lowseq.R: Lowseq variant plotting
- 2_scATAC_giotto_obj_creation.Rmd: Giotto object creation (reference)
- 3_ArchR_ATAC_analysis.Rmd: Basic ArchR analysis template (reference)

## Workflow

1. **Fragment files created** (from `build_tissue/` pipeline)
2. **ArchR objects created** (0_create_archr_qc_cluster.R) ← YOU ARE HERE
3. **Alleloscope CNV analysis** (uses clusters from step 2)
4. **Variant analysis** (uses ArchR objects from step 2)

## Output Organization

```
Data/01_outputs/archR_objects/
├── deepseq_488B/
│   ├── arrows/ (Arrow files)
│   ├── *_archR_project_final/ (Final ArchR project)
│   └── *_archR_final.rds (Quick-load RDS)
...

analysis/plots/cnv_analysis/
├── archR_qc_deepseq_488B.pdf
├── archR_qc_deepseq_489.pdf
├── archR_qc_lowseq_488B.pdf
├── archR_qc_lowseq_489.pdf
└── archR_processing_summary.tsv
```

## Notes
- Paths should be validated before qsub submission.
- Prefer absolute paths in qsub scripts.
- ⭐ Use 0_create_archr_qc_cluster.R for standard preprocessing
- All scripts use `hg38` genome (configure in script if needed)
- Support both 4 individual objects (tissue×depth) and combined analyses
