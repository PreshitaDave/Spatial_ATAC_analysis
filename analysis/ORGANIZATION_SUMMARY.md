# Analysis Folder Organization (2026-05-13)

## Overview
The analysis folder contains scripts, logs, and plots organized by pipeline stage and analysis type.

## Directory Structure

```
analysis/
├── qsub_logs/                              # Job logs organized by category
│   ├── build_tissue/ (9 files)
│   │   ├── build_tissue.e5612448
│   │   ├── build_tissue.e5612462
│   │   ├── build_tissue.e5612496 (latest)
│   │   ├── build_tissue.o5612448
│   │   ├── build_tissue.o5612462
│   │   ├── build_tissue.o5612496 (latest - SUCCESS)
│   │   ├── build_tissue_barcodes_20260513_210310.log
│   │   ├── build_tissue_barcodes_run_20260513_210054.log
│   │   └── build_tissue_files.20260505_232519.log
│   │
│   ├── diagnostics/ (14 files)
│   │   ├── edge_check.* (4 files)
│   │   ├── edge_nfrags.* (2 files)
│   │   ├── edge_nfrags_unsup.* (2 files)
│   │   ├── edge_regen.* (2 files)
│   │   ├── edge_unsup.* (2 files)
│   │   ├── somatic_compare.err
│   │   ├── somatic_deep.* (2 files)
│   │   ├── somatic.* (3 files: .e4113298, .e4134239, .e4210953, .o4113298, .o4134239, .o4210953)
│   │   ├── variant_qc_*.err (3 files)
│   │   ├── archr_lowseq_var.o4328556
│   │   ├── archr_variants.err
│   │   ├── germline_deepseq (archived)
│   │   ├── germline_deepseq_error (archived)
│   │   ├── preprocess (archived)
│   │   ├── preprocess_error (archived)
│   │   ├── pyscistopic.* (multiple)
│   │   ├── trim.* (multiple)
│   │   ├── download.* (multiple)
│   │   ├── map.* (multiple)
│   │   └── bigwig (archived)
│   │
│   ├── archived/ (12 files)
│   │   ├── dbg_exact_nfr.*
│   │   ├── dbg_nfrags_awk.*
│   │   ├── dbg_pipefail.*
│   │   ├── diagnose_deepseq.*
│   │   ├── test_build_tissue_edge.log
│   │   └── [other debug files]
│   │
│   └── README.md (job log organization guide)
│
├── src/                                    # Analysis scripts organized by pipeline
│   ├── build_tissue/                       # Edge-effect filtering & barcode generation
│   │   ├── build_tissue_barcode_lists.R
│   │   ├── build_tissue_barcodes_edge_nfrags_plots.R ⭐ (main pipeline)
│   │   └── README.md
│   │
│   ├── cnv_calling/                        # CNV and copy number analysis
│   │   ├── alleloscope/                    # Alleloscope for haplotype phasing
│   │   │   └── [alleloscope scripts]
│   │   ├── numbat/                         # NUMBAT for CNV calling
│   │   │   └── [numbat scripts]
│   │   └── README.md
│   │
│   ├── data_org/                           # Data organization utilities
│   │   └── [data organization scripts]
│   │
│   ├── initial_preprocess/                 # Raw data preprocessing
│   │   └── [preprocessing scripts]
│   │
│   ├── pipeline/                           # Main analysis pipelines
│   │   ├── archr/                          # ArchR ATAC-seq analysis
│   │   │   ├── 10_archr_variant_plotting_deepseq.R
│   │   │   ├── 10b_archr_variant_plotting_lowseq.R
│   │   │   ├── 2_scATAC_giotto_obj_creation.Rmd
│   │   │   ├── 3_ArchR_ATAC_analysis.Rmd
│   │   │   └── README.md
│   │   ├── other/                         # Miscellaneous scripts
│   │   ├── plotting/                      # Visualization utilities
│   │   ├── somatic/                       # Somatic variant analysis
│   │   │   ├── 9_somatic_snv_comparison.R
│   │   │   ├── 8_deepseq_variant_analysis.Rmd
│   │   │   ├── 8_variant_qc_comparison.R
│   │   │   └── README.md
│   │   └── variant_calling/               # Variant calling pipelines
│   │
│   ├── pycistopic/                         # Topic modeling for scATAC
│   │   └── [jupyter notebooks]
│   │
│   ├── variant_calling/                    # External variant calling tools
│   │   └── monopgen/
│   │       ├── deepseq/
│   │       ├── lowseq/
│   │       │   ├── 7_variant_calling_atac.Rmd
│   │       │   └── [other scripts]
│   │       └── README.md
│   │
│   └── [Top-level scripts]
│       ├── 1_Lib_install.Rmd
│       ├── 2_scATAC_giotto_obj_creation.Rmd
│       ├── 3_ArchR_ATAC_analysis.Rmd
│       ├── 4_CNV_calling_ATAC.Rmd
│       ├── 5_Xenium_giotto_analysis.Rmd
│       ├── 6_compare_atac_rna.Rmd
│       ├── 7_geseca_analysis_atac.Rmd
│       ├── 8_deepseq_variant_analysis.Rmd
│       └── [etc...]
│
└── plots/                                  # Analysis outputs organized by pipeline
    ├── edge_effect/                        # 🆕 Edge-effect filtering visualizations
    │   ├── deepseq_488B/
    │   │   ├── deepseq_488B_before_edge_filter.png (459K)
    │   │   ├── deepseq_488B_after_edge_filter.png (482K)
    │   │   └── deepseq_488B_nFrags_hist_cutoff.png (38K)
    │   ├── deepseq_489/
    │   ├── lowseq_488B/
    │   └── lowseq_489/
    │
    ├── cnv_analysis/                       # 🆕 CNV calling results (Alleloscope, NUMBAT)
    │   ├── alleloscope/
    │   ├── numbat/
    │   └── [analysis plots]
    │
    ├── comparison/                         # Comparative analysis plots
    │   ├── variant_qc/                     # Variant QC comparisons
    │   │   ├── chromR_plot_*.pdf
    │   │   ├── chromoqc_*.pdf
    │   │   ├── dp_distribution_*.pdf
    │   │   ├── genotype_quality_*.pdf
    │   │   ├── qual_distribution_*.pdf
    │   │   ├── variant_counts_comparison.pdf
    │   │   ├── variant_type_*.pdf
    │   │   ├── phased/
    │   │   └── png_preview/
    │   │
    │   ├── somatic_characterization/       # Somatic variant characterization
    │   │   └── [somatic analysis plots]
    │   │
    │   └── somatic_comparison_old/         # Legacy somatic comparisons
    │       └── [old comparison plots]
    │
    ├── README.md (plots organization guide)
    └── ORGANIZATION_SUMMARY.md (this file)
```

## Pipeline Organization

### 1. Build Tissue (Barcode Filtering)
```
Input: Fragment files (01_inputs/fragments/)
Process: build_tissue_barcodes_edge_nfrags_plots.R
Output: 
  - Barcodes: Data/01_inputs/barcodes/tissue_barcodes/{object}/
  - Plots: plots/edge_effect/{object}/
  - Summary: edge_effect_nfrags_thresholds.tsv
```

### 2. CNV Calling
```
Input: Fragments + Barcodes + (optional BAM)
Tools: Alleloscope, NUMBAT
Output: plots/cnv_analysis/
```

### 3. Variant Analysis
```
Input: VCFs, BAMs, Fragments
Tools: ArchR, MonoPogen, NUMBAT
Output: plots/comparison/variant_qc/
```

### 4. Somatic SNV Comparison
```
Input: Deepseq & Lowseq variant calls
Process: somatic_snv_comparison.R
Output: plots/comparison/somatic_characterization/
```

## Job Log Categories

| Folder | Purpose | Count |
|--------|---------|-------|
| build_tissue | Edge-effect barcode filtering | 9 |
| diagnostics | Edge checks, variant QC, tests | 14 |
| archived | Debug runs, test scripts | 12 |

## Key Findings (Latest Run)

### Build Tissue Job (ID: 5612496)
- **Status**: ✅ SUCCESS
- **Date**: 2026-05-13 21:03 - 21:34:59
- **Objects processed**: 4/4 (all complete)
- **Output files**: 
  - Barcodes: 17 TSV files + 1 summary
  - Plots: 12 PNG files (4 samples × 3 plot types)

### Edge Effect Statistics
- Total barcodes by object: 11,645 / 4,671
- Average retention: 98.8%
- Edge barcodes removed: 49-178 per sample

## Script Usage Guidelines

### For Alleloscope Runs:
```bash
# Use barcodes from:
Data/01_inputs/barcodes/tissue_barcodes/{object}/{object}.no_edge_effect.barcodes.tsv

# Source fragments:
Data/01_inputs/fragments/{object}/{object}.fragments.sort.filtered.bed.gz
```

### For NUMBAT:
```bash
# Reference BAM:
Data/01_inputs/bam/{object}.bam.lnk

# Phasing output → plots/cnv_analysis/numbat/
```

### For Somatic Comparison:
```bash
# Combine deepseq and lowseq results
# Output → plots/comparison/somatic_characterization/
```

## Maintenance

### Archive Old Runs
When starting new analyses, archive old plots:
```bash
mv analysis/plots/comparison analysis/plots/comparison_$(date +%Y%m%d_%H%M%S)
```

### Clean Job Logs
Periodically archive old job logs to `qsub_logs/archived/`

### Update Documentation
Keep this file in sync with changes to:
- Script locations
- Pipeline modifications
- New analysis types

