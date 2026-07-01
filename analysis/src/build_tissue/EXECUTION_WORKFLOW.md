# Spatial ATAC Edge Effect Filtering & ArchR QC Workflow

**Last Updated**: May 18, 2026  
**Project Root**: `/projectnb/paxlab/presh/projects/spatial_atac`

---

## File Organization

### Core Scripts in This Directory
```
analysis/src/build_tissue/
├── 0_archr_qc_filter_and_plot.R ....................... ⭐ MAIN QC: Filter doublets + generate plots (POST doublet-scoring)
├── 0_archr_plots_qc.R ................................ Reference: Full pipeline (includes addDoubletScores)
├── build_tissue_barcode_lists.R ....................... Initial barcode list generation
├── build_tissue_barcodes_edge_nfrags_plots.R .......... ⭐ MAIN: Edge detection + nFrags filtering
├── build_tissue_files.sh ............................... Data organization wrapper
├── organize_inputs_barcodes_bam_fragments.sh .......... Input data reorganization
├── run_archr_plots_qc.qsub.sh ......................... Job submission for QC filtering & plotting
├── run_build_tissue_barcodes_edge_nfrags_plots.qsub.sh  Job submission for edge filtering
├── run_build_tissue_files.qsub.sh ..................... Job submission for data prep
└── EXECUTION_WORKFLOW.md (this file) .................. Workflow documentation
```

---

## Execution Order

### **Phase 1: Edge Effect Detection & Barcode Filtering**

#### Step 1: Identify Spatial Edges & Fragment Outliers
**Script**: `build_tissue_barcodes_edge_nfrags_plots.R`  
**Submission**: `run_build_tissue_barcodes_edge_nfrags_plots.qsub.sh`

**What it does**:
- Reads spatial coordinates from tissue positions file
- Counts fragments per barcode from fragment BED files
- Identifies edge-effect cells (spatial periphery + high fragments)
- Generates barcode lists (all/edge/no_edge)
- Creates before/after spatial plots and nFrags histograms

**Input**:
- `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/fragments/{deepseq,lowseq}_{488B,489}/...bed.gz`
- Spatial coordinates file (tissue_positions_list.csv)

**Output**:
- Barcode files: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/{sample}/*.barcodes.tsv`
- Plots: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/edge_effect/{sample}/*.png`
- Thresholds: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/edge_effect_nfrags_thresholds.tsv`

**Run Command**:
```bash
cd /projectnb/paxlab/presh/projects/spatial_atac
qsub analysis/src/build_tissue/run_build_tissue_barcodes_edge_nfrags_plots.qsub.sh
```

**Runtime**: ~5-10 minutes  
**Status**: ✅ Complete (last run: 2026-05-15)

---

### **Phase 2: ArchR Object Creation** 
*(Prerequisites: Phase 1 must complete)*

#### Step 2: Create Arrow Files & Initial ArchR Projects
**Location**: `analysis/src/pipeline/archr/` or via ArchR creation scripts  
**Input**: Fragment files from `Data/01_inputs/fragments/`

**What it does**:
- Converts filtered fragment BED files to ArchR Arrow format
- Creates initial ArchR project objects

**Output**:
- Arrow files: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/arrow/{sample}.arrow`
- Projects: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_outputs/archR_objects/{sample}/`

**Status**: ✅ Complete

---

### **Phase 3: Quality Control, Doublet Filtering & Plotting**
*(Prerequisites: Phase 2 must complete; assumes doublet scores already computed)*

#### Step 3: Filter Doublets & Generate QC Plots
**Primary Script**: `0_archr_qc_filter_and_plot.R` ⭐ **USE THIS ONE**  
**Job Submission**: `run_archr_plots_qc.qsub.sh`  
**Location**: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/build_tissue/`

**What it does** (Post-doublet-scoring workflow):
1. Applies basic QC filters (TSS ≥ 3, nFrags ≥ 1000)
2. Uses ArchR's `filterDoublets()` function to filter based on DoubletEnrichment metric
   - DoubletEnrichment = (doublets near cell) / (expected from random distribution)
   - Uses standard parameters: `filterRatio = 1.0`, `cutEnrich = 1`, `cutScore = 2`
3. Extracts final cell barcodes and saves in ArchR object (`proj$FinalCellBarcodes`)
4. Saves barcode lists to files: `{sample}_final_cell_barcodes.txt`
5. Generates publication-ready QC plots (4 panels: TSS dist., nFrags dist., scatter, UMAP)
6. Creates summary statistics TSV with cell counts at each stage

**Why This Script**:
- ✅ Doublet scores already computed (from ArchR creation pipeline)
- ✅ This just applies filtering threshold and generates plots
- ✅ Faster (~5-10 min) than full pipeline
- ✅ Uses ArchR's proper DoubletEnrichment metric (not manual thresholding)
- ✅ Saves filtered barcodes for downstream analysis

**Input**:
- RDS files: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_outputs/archR_objects/{sample}/{sample}_archR_final.rds`
  - Must have doublet scores pre-computed by `addDoubletScores()`

**Output**:
- **Updated ArchR RDS objects** (in-place):
  - `Data/01_outputs/archR_objects/{sample}/{sample}_archR_final.rds` ✅ Updated with `FinalCellBarcodes` field
- **Barcode lists** (final cells only):
  - `Data/01_outputs/archR_objects/{sample}/{sample}_final_cell_barcodes.txt`
- **QC PDFs** (4 files, one per sample):
  - `analysis/plots/archr_obj/archR_qc_{sample}.pdf` (TSS, nFrags, scatter, UMAP)
- **Summary statistics TSV**:
  - `analysis/plots/archr_obj/archR_processing_summary.tsv` (cell counts at each stage)

**Run Command** (interactive):
```bash
cd /projectnb/paxlab/presh/projects/spatial_atac
Rscript analysis/src/build_tissue/0_archr_qc_filter_and_plot.R
```

**Or submit as job** (recommended):
```bash
qsub analysis/src/build_tissue/run_archr_plots_qc.qsub.sh
```

**Runtime**: ~5-10 minutes  
**Status**: 🔄 In Progress (Job ID: 5711413, submitted 2026-05-18)

---

**Reference Scripts** (for documentation only):
- `0_archr_plots_qc.R` - Original script (includes full pipeline with addDoubletScores); kept for reference

---

### **Phase 4: Edge-Effect Filtered ArchR Projects (Optional)**
*(Prerequisites: Phase 1 and 2 must complete)*

#### Step 4: Create ArchR Objects with Edge-Effect Filtering
**Script**: `save_archr_tissue_no_edge.R` (in `analysis/src/data_org/`)  
**Submission**: `run_save_archr_tissue_no_edge.qsub.sh`

**What it does**:
- Loads full ArchR projects
- Subsets to only keep barcodes in `*.no_edge_effect.barcodes.tsv`
- Saves filtered projects to new directory

**Input**:
- Barcode files from Phase 1
- ArchR projects from Phase 2

**Output**:
- Filtered projects: `/projectnb/paxlab/presh/projects/spatial_atac/Data/archr_tissue_no_edge/{dataset}_{tissue}_no_edge/`

**Status**: ✅ Complete

---

## Summary Statistics

### Cell Count Progression

| Tissue | Start | After Edge Filter | After nFrags Filter | Final (ArchR QC) |
|--------|-------|-------------------|---------------------|------------------|
| deepseq_488B | 11,645 | 11,467 (-1.5%) | - | 10,038 (-13.8%) |
| deepseq_489 | 4,671 | 4,622 (-1.0%) | - | 3,978 (-14.8%) |
| lowseq_488B | 11,645 | 11,612 (-0.3%) | - | 10,061 (-13.6%) |
| lowseq_489 | 4,671 | 4,622 (-1.0%) | - | 3,966 (-15.0%) |

**Loss Sources**:
1. Edge-effect filtering: ~1.0% removed (spatial periphery + high fragments)
2. ArchR QC filtering: ~13-15% removed (TSS enrichment < 3, nFrags < 1000)

---

## Key Input/Output Paths

### Fragment Files (Input to Phase 1)
```
Data/01_inputs/fragments/deepseq_488B/deepseq_488B.fragments.sort.filtered.bed.gz
Data/01_inputs/fragments/deepseq_489/deepseq_489.fragments.sort.filtered.bed.gz
Data/01_inputs/fragments/lowseq_488B/lowseq_488B.fragments.sort.filtered.bed.gz
Data/01_inputs/fragments/lowseq_489/lowseq_489.fragments.sort.filtered.bed.gz
```

### Barcode Files (Output from Phase 1 → Input to Phase 4)
```
Data/01_inputs/barcodes/tissue_barcodes/deepseq_488B/
  ├── deepseq_488B.barcodes.tsv (all: 11,645)
  ├── deepseq_488B.edge_effect.barcodes.tsv (removed: 178)
  └── deepseq_488B.no_edge_effect.barcodes.tsv (kept: 11,467)
[Same structure for deepseq_489, lowseq_488B, lowseq_489]
```

### ArchR Projects (Output from Phase 2 → Input to Phase 3)
```
Data/01_outputs/archR_objects/deepseq_488B/deepseq_488B_archR_final.rds
Data/01_outputs/archR_objects/deepseq_489/deepseq_489_archR_final.rds
Data/01_outputs/archR_objects/lowseq_488B/lowseq_488B_archR_final.rds
Data/01_outputs/archR_objects/lowseq_489/lowseq_489_archR_final.rds
```

### QC Outputs (Phase 3)
```
analysis/plots/archr_obj/
  ├── archR_qc_deepseq_488B.pdf
  ├── archR_qc_deepseq_489.pdf
  ├── archR_qc_lowseq_488B.pdf
  ├── archR_qc_lowseq_489.pdf
  └── qc_statistics_detailed.csv

analysis/plots/edge_effect/
  ├── deepseq_488B/{before, after, histogram}.png
  ├── deepseq_489/{before, after, histogram}.png
  ├── lowseq_488B/{before, after, histogram}.png
  └── lowseq_489/{before, after, histogram}.png
```

---

## Troubleshooting

### Phase 1 Fails: Edge Detection Script
- **Check**: Fragment files exist and are readable
  ```bash
  ls -lh Data/01_inputs/fragments/*/
  ```
- **Check**: Spatial coordinates file exists
  ```bash
  ls -lh Data/tissue_positions_list.csv* Data/01_inputs/spatial/tissue_positions_list.csv*
  ```
- **Check**: Output directory is writable
  ```bash
  mkdir -p Data/01_inputs/barcodes/tissue_barcodes && touch Data/01_inputs/barcodes/tissue_barcodes/.test
  ```

### Phase 3 Fails: QC Script
- **Check**: RDS files exist from Phase 2
  ```bash
  ls -lh Data/01_outputs/archR_objects/*/
  ```
- **Check**: R packages installed (ArchR, ggplot2, gridExtra)
  ```bash
  Rscript -e "library(ArchR); library(ggplot2); library(gridExtra); cat('OK\n')"
  ```

---

## Notes for Future Runs

- Phase 1 is **deterministic** (same input → same output)
- Phase 1 parameters can be adjusted via environment variables (see `build_tissue_barcodes_edge_nfrags_plots.R` for defaults)
- Phase 3 is **recommended to run after every ArchR project update** to validate QC metrics
- **Do NOT delete** barcode files from Phase 1 - they are used for Phase 4 and downstream analysis
