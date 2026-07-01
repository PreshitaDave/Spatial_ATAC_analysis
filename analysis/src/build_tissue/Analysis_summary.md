# Spatial ATAC Analysis: Edge Effect Filtering & ArchR QC Summary

**Generated**: May 18, 2026  
**Project Root**: `/projectnb/paxlab/presh/projects/spatial_atac`

## Key Analysis Scripts

### Stage 0: Create Input BAM & Fragment Files
- `/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/build_tissue/build_tissue_files.sh` ⭐ **Merger script** - Merges per-sample BAMs into tissue-level BAMs, creates tissue fragment files (bgzipped + indexed)
- `/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/build_tissue/run_build_tissue_files.qsub.sh` - Job submission wrapper for build_tissue_files.sh

**Output**:
- Tissue BAM files: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/bam/{dataset}_{tissue}.bam`
- Tissue Fragment files: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/fragments/{dataset}_{tissue}/*.bed.gz`

### Stage 0b: Organize Input Data
- `/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/build_tissue/organize_inputs_barcodes_bam_fragments.sh` - Organizes Data/01_inputs for processing (moves fragments to per-object folders, renames symlinks with .lnk suffix, archives old files)
- `/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/data_org/prepare_lowseq_alleloscope_tissue_from_existing.R` - Rebuilds barcode lists and fragment caches

### Stage 1: Edge Effect & Barcode Processing
- `/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/build_tissue/build_tissue_barcodes_edge_nfrags_plots.R` ⭐ **Main edge-effect detection** - Identifies edge cells, calculates nFrags thresholds, generates barcode lists
- `/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/build_tissue/build_tissue_barcode_lists.R` - Creates initial barcode lists from spatial coordinates
- `/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/build_tissue/run_build_tissue_barcodes_edge_nfrags_plots.qsub.sh` - Job submission wrapper

### Stage 2: ArchR Object Creation & QC
- `/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/build_tissue/0_archr_plots_qc.R` ⭐ **Main QC script** - Generates publication-ready QC plots for all samples
- `/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/data_org/save_archr_tissue_no_edge.R` - Saves edge-effect filtered ArchR projects

---

## Execution Order Summary

```
Stage 0: Create BAM & Fragments
         ↓
Stage 0b: Organize & Cache Data
         ↓
Stage 1: Edge Detection & Barcode Filtering ⭐ MAIN
         ↓
Stage 2: ArchR Object Creation
         ↓
Stage 3: Generate QC Plots & Statistics ⭐ MAIN
         ↓
Stage 4 (Optional): Create Edge-Filtered ArchR Projects
```

**Recommended execution order for full pipeline**:
1. Ensure external fragment files exist in `/projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D1942/`
2. Submit Stage 0: `qsub analysis/src/build_tissue/run_build_tissue_files.qsub.sh`
3. Submit Stage 1: `qsub analysis/src/build_tissue/run_build_tissue_barcodes_edge_nfrags_plots.qsub.sh`
4. (Automatically) Stage 2 occurs via external ArchR creation workflow
5. Run Stage 3: `Rscript analysis/src/build_tissue/0_archr_plots_qc.R`
6. (Optional) Run Stage 4: `qsub analysis/src/data_org/run_save_archr_tissue_no_edge.qsub.sh`

---

## Input Files & Their Creation

### Source Data (External)
These files come from the Atlasxomics sequencing pipeline:
- **deepseq**: `/projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D1942/D1942_deepseq/D01942_NG06549_ATv008_1/fragments.sort.bed.gz`
- **lowseq**: `/projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D1942/D1942_lowseq/fragments.tsv.gz`

### Stage 0 Output: Tissue-Level BAM & Fragment Files
These are created by `build_tissue_files.sh` by merging per-sample files:

**BAM Files** (merged from per-sample BAMs by chromosome):
```
Data/01_inputs/bam/
├── deepseq_488B.bam
├── deepseq_489.bam
├── lowseq_488B.bam
└── lowseq_489.bam
```

**Fragment Files** (bgzipped + tabix indexed):
```
Data/01_inputs/fragments/
├── deepseq_488B/deepseq_488B.fragments.sort.filtered.bed.gz
├── deepseq_488B/deepseq_488B.fragments.sort.filtered.bed.gz.tbi
├── deepseq_489/deepseq_489.fragments.sort.filtered.bed.gz
├── deepseq_489/deepseq_489.fragments.sort.filtered.bed.gz.tbi
├── lowseq_488B/lowseq_488B.fragments.sort.filtered.bed.gz
├── lowseq_488B/lowseq_488B.fragments.sort.filtered.bed.gz.tbi
├── lowseq_489/lowseq_489.fragments.sort.filtered.bed.gz
└── lowseq_489/lowseq_489.fragments.sort.filtered.bed.gz.tbi
```

**Spatial Coordinates** (used for edge detection):
```
Data/tissue_positions_list.csv (or Data/01_inputs/spatial/tissue_positions_list.csv)
```

### Stage 1 Output: Barcode Files (Barcode Filtering)

### Deepseq 488B
- All barcodes: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/deepseq_488B/deepseq_488B.barcodes.tsv`
- Edge effect: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/deepseq_488B/deepseq_488B.edge_effect.barcodes.tsv`
- No edge effect: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/deepseq_488B/deepseq_488B.no_edge_effect.barcodes.tsv`
- Fragment counts: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/deepseq_488B/deepseq_488B_nFrags_from_fragments.tsv.gz`

### Deepseq 489
- All barcodes: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/deepseq_489/deepseq_489.barcodes.tsv`
- Edge effect: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/deepseq_489/deepseq_489.edge_effect.barcodes.tsv`
- No edge effect: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/deepseq_489/deepseq_489.no_edge_effect.barcodes.tsv`
- Fragment counts: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/deepseq_489/deepseq_489_nFrags_from_fragments.tsv.gz`

### Lowseq 488B
- All barcodes: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/lowseq_488B/lowseq_488B.barcodes.tsv`
- Edge effect: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/lowseq_488B/lowseq_488B.edge_effect.barcodes.tsv`
- No edge effect: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/lowseq_488B/lowseq_488B.no_edge_effect.barcodes.tsv`
- Fragment counts: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/lowseq_488B/lowseq_488B_nFrags_from_fragments.tsv.gz`

### Lowseq 489
- All barcodes: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/lowseq_489/lowseq_489.barcodes.tsv`
- Edge effect: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/lowseq_489/lowseq_489.edge_effect.barcodes.tsv`
- No edge effect: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/lowseq_489/lowseq_489.no_edge_effect.barcodes.tsv`
- Fragment counts: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/lowseq_489/lowseq_489_nFrags_from_fragments.tsv.gz`

### Stage 2 Output: ArchR Objects & Projects (Input to Stage 3)

These are created by external ArchR creation workflow and used as input for QC plotting:

**Arrow Files**:
```
Data/01_inputs/arrow/
├── deepseq_488B.arrow
├── deepseq_489.arrow
├── lowseq_488B.arrow
└── lowseq_489.arrow
```

**ArchR Project Objects**:
```
Data/01_outputs/archR_objects/
├── deepseq_488B/
│   ├── deepseq_488B_archR_final.rds ⭐ (used by Stage 3)
│   ├── deepseq_488B_archR_project/
│   └── deepseq_488B_archR_project_final/
├── deepseq_489/
│   ├── deepseq_489_archR_final.rds ⭐ (used by Stage 3)
│   ├── deepseq_489_archR_project/
│   └── deepseq_489_archR_project_final/
├── lowseq_488B/
│   ├── lowseq_488B_archR_final.rds ⭐ (used by Stage 3)
│   ├── lowseq_488B_archR_project/
│   └── lowseq_488B_archR_project_final/
└── lowseq_489/
    ├── lowseq_489_archR_final.rds ⭐ (used by Stage 3)
    ├── lowseq_489_archR_project/
    └── lowseq_489_archR_project_final/
```

---

### Stage 3 Output: QC Plots & Statistics (Final Deliverables)

**ArchR QC PDFs** (publication-ready, 4 plots per file):
- Deepseq 488B: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/archr_obj/archR_qc_deepseq_488B.pdf`
- Deepseq 489: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/archr_obj/archR_qc_deepseq_489.pdf`
- Lowseq 488B: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/archr_obj/archR_qc_lowseq_488B.pdf`
- Lowseq 489: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/archr_obj/archR_qc_lowseq_489.pdf`

**QC Statistics CSV** (cell counts & metrics):
- `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/archr_obj/qc_statistics_detailed.csv`

---

### Stage 4 Output: Edge-Filtered ArchR Projects (Optional)

**Edge-Filtered Project Directories** (subsets of Stage 2 projects keeping only no_edge_effect barcodes):
```
Data/archr_tissue_no_edge/
├── deepseq_488B_no_edge/
├── deepseq_489_no_edge/
├── lowseq_488B_no_edge/
└── lowseq_489_no_edge/
```

These are used for downstream CNV calling and multi-omic analysis where edge-effect cells should be excluded.

## Quality Control Plots

### ArchR QC PDFs (Publication-Ready)
- Deepseq 488B: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/archr_obj/archR_qc_deepseq_488B.pdf`
- Deepseq 489: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/archr_obj/archR_qc_deepseq_489.pdf`
- Lowseq 488B: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/archr_obj/archR_qc_lowseq_488B.pdf`
- Lowseq 489: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/archr_obj/archR_qc_lowseq_489.pdf`

### Edge Effect Filtering Plots

#### Deepseq 488B
- Before filter: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/edge_effect/deepseq_488B/deepseq_488B_before_edge_filter.png`
- After filter: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/edge_effect/deepseq_488B/deepseq_488B_after_edge_filter.png`
- nFrags histogram: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/edge_effect/deepseq_488B/deepseq_488B_nFrags_hist_cutoff.png`

#### Deepseq 489
- Before filter: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/edge_effect/deepseq_489/deepseq_489_before_edge_filter.png`
- After filter: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/edge_effect/deepseq_489/deepseq_489_after_edge_filter.png`
- nFrags histogram: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/edge_effect/deepseq_489/deepseq_489_nFrags_hist_cutoff.png`

#### Lowseq 488B
- Before filter: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/edge_effect/lowseq_488B/lowseq_488B_before_edge_filter.png`
- After filter: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/edge_effect/lowseq_488B/lowseq_488B_after_edge_filter.png`
- nFrags histogram: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/edge_effect/lowseq_488B/lowseq_488B_nFrags_hist_cutoff.png`

#### Lowseq 489
- Before filter: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/edge_effect/lowseq_489/lowseq_489_before_edge_filter.png`
- After filter: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/edge_effect/lowseq_489/lowseq_489_after_edge_filter.png`
- nFrags histogram: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/edge_effect/lowseq_489/lowseq_489_nFrags_hist_cutoff.png`

## Edge Effect Thresholds

File: `/projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/barcodes/tissue_barcodes/edge_effect_nfrags_thresholds.tsv`

| Dataset | Tissue | Total Cells | Edge Cells | Kept Cells | nFrags Threshold |
|---------|--------|-------------|------------|-----------|------------------|
| deepseq | 488B   | 11,645      | 178        | 11,467    | 240,000          |
| deepseq | 489    | 4,671       | 49         | 4,622     | 130,000          |
| lowseq  | 488B   | 11,645      | 33         | 11,612    | 70,000           |
| lowseq  | 489    | 4,671       | 49         | 4,622     | 50,000           |

## Cell Count Summary

File: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/plots/archr_obj/qc_statistics_detailed.csv`

**QC Filtering Pipeline**:
- **Basic QC**: TSS enrichment ≥ 3 AND nFrags ≥ 1000 (removes low-quality cells)
- **Doublet Detection**: ArchR's doublet scoring algorithm (flags potential doublets)
- **Final Filter**: Basic QC passed AND NOT flagged as doublet

| Tissue | Initial Cells | After Barcode Filter | After Basic QC | # Doublets Detected | Doublet Rate % | Final Cells | Overall Retention % |
|--------|---------------|----------------------|-----------------|----------------------|-----------------|-------------|-------------------|
| deepseq_488B | 11,382 | 11,359 | 11,319 | 1,281 | 11.3% | 10,038 | 88.2% |
| deepseq_489 | 4,586 | 4,537 | 4,150 | 172 | 4.1% | 3,978 | 86.7% |
| lowseq_488B | 11,382 | 11,349 | 11,348 | 1,287 | 11.3% | 10,061 | 88.4% |
| lowseq_489 | 4,211 | 4,162 | 4,137 | 171 | 4.1% | 3,966 | 94.2% |

**Important Notes**:
- ✅ **YES, QC filter INCLUDES doublet detection** - Cells flagged as doublets by ArchR's `addDoubletScores()` are removed in the final filtering step
- **Doublet Rate** = (Doublets Detected / After Basic QC) × 100:
  - deepseq_488B: 11.3% doublet rate
  - deepseq_489: 4.1% doublet rate  
  - lowseq_488B: 11.3% doublet rate
  - lowseq_489: 4.1% doublet rate
- Higher doublet rates in 488B samples (~11%) suggest either higher cell density or more cell aggregation compared to 489 samples (~4%)

## Processing Pipeline

### Stage 0: Prepare Input Files (Create BAM & Fragments)
**Script**: `build_tissue_files.sh`  
**Job Submission**: `run_build_tissue_files.qsub.sh`

**What it does**:
1. Reads external fragment files from Atlasxomics output
2. Merges BAM files from all samples by tissue region (488B vs 489)
3. Creates tissue-level fragment BED files (compressed + indexed with tabix)
4. Stores BAMs in `Data/01_inputs/bam/` and fragments in `Data/01_inputs/fragments/{object}/`

**Input**: External Atlasxomics data files

**Output**: 
- BAM files: `Data/01_inputs/bam/{dataset}_{tissue}.bam`
- Fragment files: `Data/01_inputs/fragments/{dataset}_{tissue}/{dataset}_{tissue}.fragments.sort.filtered.bed.gz`

**Run Command**:
```bash
cd /projectnb/paxlab/presh/projects/spatial_atac
qsub analysis/src/build_tissue/run_build_tissue_files.qsub.sh
```

**Runtime**: ~20-30 minutes  
**Status**: ✅ Complete (files exist in Data/01_inputs/)

---

### Stage 0b: Organize & Cache Barcode Data
**Script**: `organize_inputs_barcodes_bam_fragments.sh`

**What it does**:
1. Organizes Data/01_inputs structure (moves fragments into per-object folders)
2. Renames symlinks with `.lnk` suffix for clarity
3. Archives old/unused files to `Data/01_inputs/archive/`

**Run Command** *(optional, mainly for cleanup)*:
```bash
bash analysis/src/build_tissue/organize_inputs_barcodes_bam_fragments.sh
```

---

### Stage 1: Spatial Filtering
**Script**: `build_tissue_barcodes_edge_nfrags_plots.R`  
**Job Submission**: `run_build_tissue_barcodes_edge_nfrags_plots.qsub.sh`

**What it does**:
1. Reads spatial coordinates and identifies edge-effect cells
2. Counts fragments per barcode from BED files
3. Calculates nFrags threshold (statistical outlier detection)
4. Marks cells as edge-effect if both: (a) spatially peripheral, (b) high fragments
5. Generates barcode lists (all/edge/no_edge)
6. Creates before/after spatial plots

**Input**:
- Fragment files: `Data/01_inputs/fragments/{dataset}_{tissue}/...bed.gz`
- Spatial file: `Data/tissue_positions_list.csv`

**Output**:
- Barcode files: `Data/01_inputs/barcodes/tissue_barcodes/{sample}/*.barcodes.tsv`
- Plots: `analysis/plots/edge_effect/{sample}/*.png`
- Thresholds: `Data/01_inputs/barcodes/tissue_barcodes/edge_effect_nfrags_thresholds.tsv`

**Run Command**:
```bash
cd /projectnb/paxlab/presh/projects/spatial_atac
qsub analysis/src/build_tissue/run_build_tissue_barcodes_edge_nfrags_plots.qsub.sh
```

**Runtime**: ~5-10 minutes  
**Status**: ✅ Complete (last run: 2026-05-15)

---

### Stage 2: ArchR Object Creation
*(Prerequisites: Stage 0 must complete)*

**Location**: `analysis/src/pipeline/archr/` or external ArchR creation scripts

**What it does**:
1. Creates Arrow files from tissue fragment BED files
2. Initializes ArchR project objects (one per sample)
3. Applies basic QC filters (TSS enrichment, fragment count)
4. Saves final RDS files for each sample

**Input**:
- Fragment files: `Data/01_inputs/fragments/{dataset}_{tissue}/*.bed.gz`
- Optional: Barcode files from Stage 1 for filtering

**Output**:
- Arrow files: `Data/01_inputs/arrow/{sample}.arrow`
- ArchR Projects (RDS): `Data/01_outputs/archR_objects/{sample}/{sample}_archR_final.rds`
- ArchR Project directories: `Data/01_outputs/archR_objects/{sample}/{sample}_archR_project_final/`

**Status**: ✅ Complete

---

### Stage 3: Quality Control, Doublet Filtering & Visualization ⭐ RECOMMENDED
**Location**: `/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/build_tissue/`

**Primary Script**: `0_archr_qc_filter_and_plot.R` ⭐ **USE THIS ONE** (Post-doublet-scoring filtering only)  
**Job Submission**: `run_archr_plots_qc.qsub.sh`

**What it does** (assumes doublet scores already computed):
1. **Applies basic QC filters** (TSS enrichment ≥ 3 AND nFrags ≥ 1000) to identify high-quality cells
2. **Filters doublets** using ArchR's `filterDoublets()` function with:
   - DoubletEnrichment metric (ratio of simulated doublets nearby to expected random distribution)
   - `filterRatio = 1.0`, `cutEnrich = 1`, `cutScore = 2` (standard ArchR parameters)
3. **Extracts final cell barcodes** and stores in ArchR object (`proj$FinalCellBarcodes`)
4. **Saves barcode lists** per sample: `{sample}_final_cell_barcodes.txt`
5. **Generates publication-ready QC plots**:
   - TSS enrichment distribution (post-filter)
   - Fragment count distribution (log10 scale, post-filter)
   - TSS vs Fragment scatter plot (post-filter with thresholds)
   - UMAP clustering (colored by cluster)
6. **Generates summary statistics** with cell counts at each filtering stage

**Doublet Filtering Details** (ArchR Pipeline):
- Uses simulated doublets projected into UMAP embedding
- Computes DoubletEnrichment for each cell: enrichment score = (doublets nearby) / (expected random)
- Filters cells with high DoubletEnrichment using `filterRatio` constraint
- More conservative than manual thresholding; avoids over-filtering

**Input**:
- ArchR RDS files: `Data/01_outputs/archR_objects/{sample}/{sample}_archR_final.rds` (must have doublet scores already computed)
  - deepseq_488B_archR_final.rds
  - deepseq_489_archR_final.rds
  - lowseq_488B_archR_final.rds
  - lowseq_489_archR_final.rds

**Output**:
- **Updated ArchR RDS objects** (same path, with `FinalCellBarcodes` field added):
  - `Data/01_outputs/archR_objects/{sample}/{sample}_archR_final.rds` ✅ Updated in-place
- **Barcode lists** (one per sample):
  - `Data/01_outputs/archR_objects/{sample}/{sample}_final_cell_barcodes.txt`
  - Contains final cell barcodes ready for downstream analysis
- **QC PDFs** (4 files, 1 per sample):
  - `analysis/plots/archr_obj/archR_qc_{sample}.pdf`
  - Shows TSS, nFrags, scatter plot, and UMAP for final filtered cells
- **Summary statistics TSV**:
  - `analysis/plots/archr_obj/archR_processing_summary.tsv`
  - Columns: tissue, initial_cells, after_barcode_filter, after_basic_qc, doublets_detected, final_cells, clusters_found, tss_enrichment_mean, tss_enrichment_sd

**Run Command** (interactive):
```bash
cd /projectnb/paxlab/presh/projects/spatial_atac
Rscript analysis/src/build_tissue/0_archr_qc_filter_and_plot.R
```

**Submit as job** (recommended):
```bash
qsub analysis/src/build_tissue/run_archr_plots_qc.qsub.sh
```

**Runtime**: ~5-10 minutes (post-doublet-scoring only)  
**Status**: 🔄 In Progress (Job ID: 5711413, submitted 2026-05-18)

---

**Alternative Script** (for reference only):
- `0_archr_plots_qc.R` - Original full pipeline (includes addDoubletScores); kept for reference only

---

### Stage 4: Create Edge-Filtered ArchR Projects (Optional)
**Script**: `save_archr_tissue_no_edge.R` (in `analysis/src/data_org/`)  
**Submission**: `run_save_archr_tissue_no_edge.qsub.sh`

**What it does**:
1. **Loads** full ArchR projects (all cells) from Stage 2
2. **Reads** edge-effect barcode lists from Stage 1
3. **Subsets** each project to keep ONLY cells in `*.no_edge_effect.barcodes.tsv`
4. **Saves** filtered projects to new directory structure
5. **Useful for**: Downstream analyses that strictly require non-edge cells (CNV calling, etc.)

**Input**:
- ArchR projects (RDS): `Data/01_outputs/archR_objects/{sample}/{sample}_archR_final.rds`
- Barcode lists: `Data/01_inputs/barcodes/tissue_barcodes/{sample}/{sample}.no_edge_effect.barcodes.tsv`

**Output**:
- Edge-filtered ArchR projects: `Data/archr_tissue_no_edge/{dataset}_{tissue}_no_edge/`
  - deepseq_488B_no_edge/
  - deepseq_489_no_edge/
  - lowseq_488B_no_edge/
  - lowseq_489_no_edge/

**Run Command**:
```bash
cd /projectnb/paxlab/presh/projects/spatial_atac
qsub analysis/src/data_org/run_save_archr_tissue_no_edge.qsub.sh
```

**Runtime**: ~10-20 minutes  
**Status**: ✅ Complete

---

## Quick Reference: File Locations by Stage