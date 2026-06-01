# ArchR QC, Clustering & Analysis Pipeline

## Overview

This pipeline creates ArchR objects from filtered fragment files, applies comprehensive QC filtering, removes doublets, performs clustering and UMAP analysis, and saves all results in organized PDFs.

## Pipeline Steps

1. **Arrow File Creation**: Converts BED fragment files to Arrow format (ArchR's efficient storage format)
2. **ArchR Project Creation**: Creates ArchR projects from Arrow files
3. **Barcode Filtering**: Filters to no-edge-effect barcodes only
4. **QC Filtering**: Removes cells failing TSS enrichment and fragment count thresholds
5. **Dimensionality Reduction**: Performs iterative LSI on tile matrix
6. **Clustering**: K-means clustering via Seurat integration
7. **UMAP Embedding**: Creates 2D UMAP embedding for visualization
8. **Doublet Filtering** (optional): Removes doublet cells if doublet scores available
9. **PDF Reporting**: Generates comprehensive QC and clustering plots

## Input Requirements

### Fragment Files (Required)
- Location: `Data/01_inputs/fragments/{object}/{object}.fragments.sort.filtered.bed.gz`
- Format: Compressed BED files with indexed `.tbi` files
- Objects: `deepseq_488B`, `deepseq_489`, `lowseq_488B`, `lowseq_489`

### Barcode Files (Required)
- Location: `Data/01_inputs/barcodes/tissue_barcodes/{object}/{object}.no_edge_effect.barcodes.tsv`
- Format: One barcode per line (with or without -1 suffix)
- Auto-normalizes to handle barcode format variations

## Output Structure

```
Data/01_outputs/archR_objects/
├── deepseq_488B/
│   ├── arrows/
│   │   └── Deepseq_488B.arrow (Arrow file)
│   ├── deepseq_488B_archR_project/ (Intermediate project)
│   ├── deepseq_488B_archR_project_final/ (Final saved project)
│   └── deepseq_488B_archR_final.rds (Quick-load RDS)
├── deepseq_489/
├── lowseq_488B/
└── lowseq_489/

analysis/plots/cnv_analysis/
├── archR_qc_deepseq_488B.pdf
├── archR_qc_deepseq_489.pdf
├── archR_qc_lowseq_488B.pdf
├── archR_qc_lowseq_489.pdf
└── archR_processing_summary.tsv
```

## Configuration Parameters

Edit in script header:
- `min_tss`: Minimum TSS enrichment score (default: 4)
- `min_frags`: Minimum fragment count (default: 1000)
- `max_frags`: Maximum fragment count (default: Inf)
- `doublet_cutoff`: Doublet score threshold (default: 2)
- `threads`: CPU threads (auto-detected from NSLOTS)

## Running the Pipeline

### Option 1: Via qsub (Recommended for full runs)
```bash
qsub analysis/qsub/pipeline/run_archR_qc_cluster.qsub.sh
```

### Option 2: Direct execution (for testing)
```bash
# On compute node (not login node!)
export PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
export NSLOTS=8
Rscript analysis/src/pipeline/archr/0_create_archr_qc_cluster.R
```

### Option 3: Single object (for debugging)
```bash
# Modify script to process only one object, then run
Rscript analysis/src/pipeline/archr/0_create_archr_qc_cluster.R
```

## Monitoring Execution

### Check job status
```bash
qstat -j <JOB_ID>
```

### View logs in real-time
```bash
tail -f analysis/qsub_logs/build_tissue/archR_qc_cluster_<JOB_ID>.log
```

### After completion
```bash
# View summary
cat analysis/plots/cnv_analysis/archR_processing_summary.tsv

# List all outputs
ls -lh Data/01_outputs/archR_objects/*/

# View PDF reports
ls -lh analysis/plots/cnv_analysis/archR_qc_*.pdf
```

## PDF Report Contents

Each PDF includes:

### Page 1: Title & Metadata
- Pipeline version
- Run timestamp
- Object name

### QC Metrics Pages:
- **TSS vs nFrags scatter plot**: Shows QC filtering boundaries
- **Fragment count histogram**: Distribution of sequencing depth
- **TSS enrichment histogram**: Distribution of TSS scores
- **QC summary bar chart**: Cells passing each threshold

### Clustering Pages:
- **UMAP by clusters**: Cell clusters colored by ID
- **UMAP by TSS**: TSS enrichment values on UMAP
- **UMAP by nFrags**: Fragment counts on UMAP
- *Additional pages as features are computed*

## Troubleshooting

### Issue: "No matching cells found!"
**Cause**: Barcode format mismatch between fragments and barcode file
**Solution**: 
- Check barcode files exist and are not empty
- Verify fragment files have expected BED format
- Check if barcodes need -1 suffix normalization

### Issue: "Missing fragments" or "Missing barcodes"
**Cause**: Input files not in expected locations
**Solution**:
- Run verification step: `ls -lh Data/01_inputs/fragments/*/`
- Run verification step: `ls -lh Data/01_inputs/barcodes/tissue_barcodes/*/`
- Check file permissions (should be readable)

### Issue: Job timeout (h_rt exceeded)
**Cause**: Processing too slow or insufficient resources
**Solution**:
- Increase h_rt in qsub script (currently 24:00:00)
- Process fewer objects at a time
- Increase thread count (pe omp parameter)

### Issue: "Memory exceeded"
**Cause**: Insufficient RAM for large projects
**Solution**:
- Increase mem_per_core in qsub (currently 8G)
- Reduce max_frags parameter to subset cells
- Process objects sequentially instead of together

### Issue: Arrow creation fails
**Cause**: Fragment file format issue or corrupted data
**Solution**:
- Verify fragments are gzipped: `file <filename>`
- Check for corrupted gzip: `gzip -t <filename>`
- Try creating indexed BED: `samtools index <filename>`

## Output File Descriptions

### RDS Objects
- **`{object}_archR_final.rds`**: Compressed ArchR object (fastest to load)
- Load with: `proj <- readRDS("path/to/file.rds")`

### Arrow Files
- **`{sample}.arrow`**: ArchR's native binary format
- Efficient storage, intermediate format

### ArchR Project Directories
- **`_archR_project_final/`**: Full ArchR project directory structure
- Includes Arrow files, LSI, clusters, UMAP, etc.

## Using Results for Downstream Analysis

### Load in R for further analysis
```r
library(ArchR)
addArchRGenome("hg38")

# Quick load from RDS
proj <- readRDS("Data/01_outputs/archR_objects/deepseq_488B/deepseq_488B_archR_final.rds")

# Access key data
cell_metadata <- getCellColData(proj)
umap <- getEmbedding(proj, "UMAP")
clusters <- proj$Clusters
```

### Extract cell barcodes by cluster
```r
for (cluster in unique(proj$Clusters)) {
  cells <- getCellNames(proj)[proj$Clusters == cluster]
  barcodes <- sub("-1$", "", sub(".*#", "", cells))
  writeLines(barcodes, sprintf("cluster_%s_barcodes.txt", cluster))
}
```

### Create pseudo-bulk from clusters
```r
# Get gene score matrix by cluster
gene_scores <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")

# Group by cluster and sum
cluster_means <- lapply(unique(proj$Clusters), function(c) {
  cells <- which(proj$Clusters == c)
  rowMeans(gene_scores[, cells, drop = FALSE])
})
```

## Next Steps (Common Workflows)

1. **Alleloscope (CNV analysis)**:
   - Use barcodes from clusters for phased analysis
   - Use ArchR objects to subset to specific cell types

2. **Somatic variant analysis**:
   - Use cluster information to stratify variants
   - Color variants by cell type on UMAP

3. **Trajectory analysis**:
   - Add pseudotime: `addTrajectory()` in ArchR
   - Correlate variants with pseudotime

4. **Integration with other modalities**:
   - Extract UMAP embeddings for plotting RNA data
   - Use cluster assignments for co-analysis

## Key Outputs for Alleloscope Integration

For running Alleloscope on cluster-specific cells:

1. Extract barcodes from clusters:
```bash
# From within R after loading project
for c in unique(proj$Clusters):
  barcodes_cluster_c <- extract_barcode_by_cluster(proj, c)
  write to: Data/01_outputs/archR_barcodes_cluster_{c}.tsv
```

2. Use with Alleloscope input:
```bash
# Alleloscope can now use these filtered barcode lists
alleloscope --barcodes Data/01_outputs/archR_barcodes_cluster_1.tsv \
            --fragments Data/01_inputs/fragments/{object}/ \
            --output Data/01_outputs/alleloscope_cluster_1/
```

## Performance Notes

- **Per-object runtime**: 30-60 minutes (8 cores, 24 GB RAM)
- **Total runtime for 4 objects**: ~2-4 hours (parallel via NSLOTS)
- **Disk space**: ~10-20 GB per object (Arrow + project)
- **Memory usage**: 8-16 GB per object during processing

## References

- ArchR Documentation: https://www.archrproject.com/
- ArchR Preprocessing: https://www.archrproject.com/articles/AP_Archr_Installation_and_Basic_Usage.html
- Seurat Clustering: https://satijalab.org/seurat/articles/pbmc3k_tutorial

## Contact

For issues or questions, check:
1. This README
2. ArchR official documentation
3. Script logs in `analysis/qsub_logs/`
4. Input file verification in `Data/01_inputs/`
