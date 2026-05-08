# Alleloscope scATAC Pipeline (Lowseq + Deepseq)

This folder contains the Alleloscope preparation and analysis scripts for the spatial scATAC data.

## Scope

- Scripts: `analysis/src/alleloscope`
- Inputs and outputs: `Data/alleloscope`
- Cluster logs: `analysis/qsub_logs`

## Working scripts (authoritative)

### Data organization wrappers

- `analysis/src/data_org/build_total_tissue_bam_fragments.qsub.sh`
  - Builds per-tissue BAM and fragments files from total tissue barcodes.
  - Uses D1942 fragment sources configured in `build_tissue_files.sh`.
- `analysis/src/data_org/run_save_archr_tissue_no_edge.qsub.sh`
  - Runs `analysis/src/data_org/save_archr_tissue_no_edge.R`.
  - Saves ArchR projects for each dataset+tissue after removing edge barcodes.
- `analysis/src/data_org/run_prepare_lowseq_alleloscope_tissue_from_existing.qsub.sh`
  - Runs `analysis/src/data_org/prepare_lowseq_alleloscope_tissue_from_existing.R`.
  - Reuses existing lowseq/deepseq Alleloscope inputs to build tissue-specific lowseq inputs (no full rerun).

### Input preparation

- `alleloscope_prep_helpers.R`
  - Core helper functions used by prep scripts:
  - barcode normalization and ordering
  - chromosome-aware VCF/BAM matching
  - VarTrix execution and matrix assembly
  - generation of `seg_table.rds`
- `prepare_alleloscope_lowseq.R`
  - Builds lowseq Alleloscope inputs in `Data/alleloscope/lowseq/`:
  - `alt_all.mtx`, `ref_all.mtx`, `var_all.vcf`, `barcodes.tsv`, `chr1000k_fragments.tsv`, `seg_table.rds`
- `prepare_alleloscope_deepseq.R`
  - Same as above for deepseq.
- `run_prepare_alleloscope_lowseq.qsub.sh`
  - SCC submission wrapper for lowseq input prep.
- `run_prepare_alleloscope_deepseq.qsub.sh`
  - SCC submission wrapper for deepseq input prep.

### Full analysis (no matched DNA)

- `run_alleloscope_lowseq.R`
  - End-to-end, reproducible lowseq analysis pipeline:
  - create object + matrix filter
  - parallelized `Est_regions` chromosome-level (with resume)
  - normal cell estimation with `Select_normal(pre_sel=TRUE)`
  - segmentation using normal-cell pseudobulk
  - parallelized `Est_regions` segment-level
  - final `Select_normal(pre_sel=FALSE)`
  - `Genotype_value` + `Genotype`
  - Step-6 style CNV coverage heatmap via `plot_scATAC_cnv`
  - saves checkpoint RDS objects after each major stage
- `run_alleloscope_lowseq.qsub.sh`
  - SCC submission wrapper for `run_alleloscope_lowseq.R`.
- `run_alleloscope_deepseq.R`
  - End-to-end, reproducible deepseq analysis pipeline with the same 10-step flow as lowseq:
  - parallelized `Est_regions`
  - no-DNA normal estimation via `Select_normal`
  - segmentation with normal-cell pseudobulk
  - `Genotype_value` + `Genotype`
  - Step-6 CNV coverage heatmap via `plot_scATAC_cnv`
  - checkpoint RDS outputs at major stages
- `run_alleloscope_deepseq.qsub.sh`
  - SCC submission wrapper for `run_alleloscope_deepseq.R`.

## Legacy / non-authoritative scripts

- `prepare_alleloscope_by_tissue.R`, `build_tissue_barcode_lists.R`
  - Useful utilities for tissue barcode derivation and subset prep.
- `smoke_test_chr1.R`
  - Small test workflow.
- `prepare_alleloscope_lowseq.R` currently contains partial exploratory object-analysis code at the bottom.
  - Use `run_alleloscope_lowseq.R` for full reproducible analysis.

## Runtime dependencies

- `module load R`
- R package library path for Alleloscope package:
  - `/projectnb/paxlab/presh/Rlibs/4.5`
- Alleloscope source/resources directory:
  - `/projectnb/paxlab/presh/software/Alleloscope`

Prep stage also requires (configured in helper scripts):

- VarTrix binary
- per-chromosome BAMs
- per-chromosome VCFs
- reference FASTA
- chromosome sizes file

## Reproducible runbook

Run from project root:

```bash
cd /projectnb/paxlab/presh/projects/spatial_atac
```

### 1) Prepare lowseq inputs (if not already present)

```bash
qsub analysis/src/alleloscope/run_prepare_alleloscope_lowseq.qsub.sh
```

### 2) Run full lowseq Alleloscope analysis

```bash
qsub analysis/src/alleloscope/run_alleloscope_lowseq.qsub.sh
```

Provenance in this workspace:

- Full lowseq analysis submission: job `4700297`

### 3) Prepare deepseq inputs (if not already present)

```bash
qsub analysis/src/alleloscope/run_prepare_alleloscope_deepseq.qsub.sh
```

### 4) Run full deepseq Alleloscope analysis

```bash
qsub analysis/src/alleloscope/run_alleloscope_deepseq.qsub.sh
```

Provenance in this workspace:

- Deepseq prep submission: job `4700348`
- Deepseq full analysis submission (held on prep): job `4700349`

### 5) Monitor

```bash
qstat -u preshita

tail -f analysis/qsub_logs/allelo_lowseq.<JOBID>.log
tail -f analysis/qsub_logs/allelo_deepseq.<JOBID>.log
```

### 6) Tissue data organization and smart tissue correction

```bash
qsub analysis/src/data_org/build_total_tissue_bam_fragments.qsub.sh
qsub analysis/src/data_org/run_save_archr_tissue_no_edge.qsub.sh
qsub analysis/src/data_org/run_prepare_lowseq_alleloscope_tissue_from_existing.qsub.sh
```

## Expected output tree

Lowseq primary directory:

- `Data/alleloscope/lowseq/output/`

Deepseq primary directory:

- `Data/alleloscope/deepseq/output/`

Key analysis outputs for either dataset (inside that dataset's `output/`):

- `output/rds/Obj_after_EM_chr.rds`
- `output/rds/Obj_after_seg.rds`
- `output/rds/Obj_after_EM_seg.rds`
- `output/rds/Obj_after_select_normal.rds`
- `output/rds/Obj_after_gtv.rds`
- `output/rds/Obj_final.rds`
- `output/rds/cov_obj.rds`

Plot outputs:

- `output/plots/step6_CNV_coverage_heatmap.png`
- `output/plots/gtype_scatter_ref_<region>.pdf`
- `output/plots/EMresults/`

Segment-level EM results are written separately to avoid clobbering chromosome-level files:

- `output/seg/rds/EMresults/`

## Notes for reproducibility

- The full lowseq analysis is intended to be restarted safely from existing outputs.
- Chromosome-level EM reuse is enabled via continuation logic in the analysis script.
- Normal cell estimation is data-driven and does not require matched DNA.
- Keep qsub logs with job IDs for provenance and troubleshooting.
