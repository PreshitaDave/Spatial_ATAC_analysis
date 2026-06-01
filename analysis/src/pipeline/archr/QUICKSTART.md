# ⚡ ArchR QC Pipeline - Quick Start Guide

## What This Does

Creates professional ArchR objects from your fragment files with:
- ✅ Automatic Arrow file creation
- ✅ QC filtering (TSS enrichment, fragment counts)
- ✅ Doublet removal
- ✅ Clustering (K-means via Seurat)
- ✅ UMAP visualization
- ✅ **Complete PDF reports for each object** (all plots in one file)

## Quick Start (5 minutes)

### Step 1: Verify Input Files Exist
```bash
cd /projectnb/paxlab/presh/projects/spatial_atac

# Check fragments
ls -lh Data/01_inputs/fragments/*/

# Check barcodes (no_edge_effect)
ls -lh Data/01_inputs/barcodes/tissue_barcodes/*/
```

### Step 2: Submit Job to Cluster
```bash
# Option A: Via qsub (recommended)
qsub analysis/qsub/pipeline/run_archR_qc_cluster.qsub.sh

# This will output: Job ID like 1234567
```

### Step 3: Monitor Progress
```bash
# Check job status (replace JOB_ID with actual number)
qstat -j JOB_ID

# Watch logs in real-time
tail -f analysis/qsub_logs/build_tissue/archR_qc_cluster_JOB_ID.log
```

### Step 4: Review Results
```bash
# When complete, check outputs:
ls -lh Data/01_outputs/archR_objects/*/

# View PDF reports:
ls -lh analysis/plots/cnv_analysis/archR_qc_*.pdf

# See summary:
cat analysis/plots/cnv_analysis/archR_processing_summary.tsv
```

## What Gets Generated

### PDFs (Ready for presentations!)
```
analysis/plots/cnv_analysis/
├── archR_qc_deepseq_488B.pdf     ← All plots for 488B
├── archR_qc_deepseq_489.pdf      ← All plots for 489
├── archR_qc_lowseq_488B.pdf      ← All plots for 488B
├── archR_qc_lowseq_489.pdf       ← All plots for 489
└── archR_processing_summary.tsv   ← Statistics table
```

Each PDF includes:
- Title page with timestamp
- TSS vs Fragments scatter plot
- Fragment count histogram
- TSS enrichment histogram
- QC summary bar chart
- UMAP by clusters
- UMAP by TSS
- UMAP by fragment counts

### ArchR Objects (Ready for downstream analysis)
```
Data/01_outputs/archR_objects/
├── deepseq_488B/
│   ├── deepseq_488B_archR_final.rds       ← Load this in R!
│   ├── deepseq_488B_archR_project_final/  ← Full ArchR project
│   └── arrows/                             ← Arrow files
...
```

## Using Results in R

### Quick Load
```r
library(ArchR)
addArchRGenome("hg38")

# Load the preprocessed object
proj <- readRDS("Data/01_outputs/archR_objects/deepseq_488B/deepseq_488B_archR_final.rds")

# Explore
print(proj)
table(proj$Clusters)
head(getCellColData(proj))
```

### Extract Barcodes by Cluster (for Alleloscope)
```r
library(ArchR)
proj <- readRDS("Data/01_outputs/archR_objects/deepseq_488B/deepseq_488B_archR_final.rds")

# Get barcodes for cluster 1
cluster_1_cells <- getCellNames(proj)[proj$Clusters == 1]
cluster_1_barcodes <- gsub("-1$", "", gsub(".*#", "", cluster_1_cells))

# Save for use with Alleloscope
writeLines(cluster_1_barcodes, "cluster_1_barcodes.txt")
```

## Troubleshooting

### Issue: Job didn't start
```bash
# Check if qsub errors exist:
cat analysis/qsub_logs/build_tissue/archR_qc_cluster_*.log | grep ERROR

# Most common: Missing module. Re-run with:
module load R
Rscript analysis/src/pipeline/archr/0_create_archr_qc_cluster.R
```

### Issue: "No matching cells found"
```bash
# Check barcode format:
head -n 5 Data/01_inputs/barcodes/tissue_barcodes/*/

# Check fragment barcodes (should match):
zcat Data/01_inputs/fragments/deepseq_488B/deepseq_488B.fragments.sort.filtered.bed.gz | head | awk '{print $4}'
```

### Issue: Job ran out of time/memory
```bash
# Check how long it ran:
qacct -j JOB_ID | grep -E "start_time|end_time"

# Increase resources in qsub script and rerun:
# Edit: analysis/qsub/pipeline/run_archR_qc_cluster.qsub.sh
# Change: -l h_rt=24:00:00 to -l h_rt=48:00:00 (or higher)
#         -pe omp 8 to -pe omp 16 (more cores)
#         mem_per_core=8G to mem_per_core=16G (more RAM)
```

## Advanced: Customize Parameters

Edit these in the script header (line ~25):
```r
min_tss <- 4              # Minimum TSS enrichment
min_frags <- 1000         # Minimum fragments
max_frags <- Inf          # Maximum fragments
doublet_cutoff <- 2       # Doublet score threshold
```

Then run:
```bash
export PROJECT_ROOT="/projectnb/paxlab/presh/projects/spatial_atac"
export NSLOTS=8
Rscript analysis/src/pipeline/archr/0_create_archr_qc_cluster.R
```

## Next Steps (After Pipeline Succeeds)

### 1. Run Alleloscope
```bash
# Use barcodes from clusters
qsub analysis/qsub/alleloscope/run_alleloscope_deepseq.qsub.sh
```

### 2. Somatic Variant Analysis
```bash
# Use ArchR objects for variant stratification
Rscript analysis/src/pipeline/somatic/9_somatic_snv_comparison.R
```

### 3. Explore Clusters
```r
# In R:
proj <- readRDS("Data/01_outputs/archR_objects/deepseq_488B/deepseq_488B_archR_final.rds")

# Get gene scores for cluster markers
library(ComplexHeatmap)
gene_scores <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
# ... perform DE analysis per cluster
```

## File Organization

All outputs follow your organized structure:
- **Code**: `analysis/src/pipeline/archr/`
- **Qsub**: `analysis/qsub/pipeline/`
- **Data**: `Data/01_outputs/archR_objects/`
- **Plots**: `analysis/plots/cnv_analysis/`
- **Logs**: `analysis/qsub_logs/build_tissue/`

## Expected Runtime

- **Single object**: 30-60 minutes
- **All 4 objects**: 60-120 minutes (parallel)
- **With 8 cores**: ~90 minutes

## Key Features of This Pipeline

✅ **Scalable**: Processes all 4 objects in parallel
✅ **Reproducible**: All parameters documented
✅ **Professional**: Beautiful PDF outputs ready for presentations
✅ **Complete**: QC + Clustering + Visualization all in one
✅ **Documented**: Full README with troubleshooting
✅ **Organized**: Outputs in logical folder structure

## Questions?

1. Check: `analysis/src/pipeline/archr/README_ArchR_Pipeline.md`
2. Check logs: `analysis/qsub_logs/build_tissue/archR_qc_cluster_*.log`
3. Check script: `analysis/src/pipeline/archr/0_create_archr_qc_cluster.R`

---

**Status**: ✅ Ready to run!  
**Prerequisite**: Fragment files created successfully (from build_tissue pipeline)  
**Next Pipeline**: Alleloscope analysis (uses ArchR clusters & barcodes)
