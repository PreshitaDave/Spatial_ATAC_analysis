# data_org workflow

This folder contains scripts that organize tissue-level files and avoid expensive reruns.

## Scripts

- `build_total_tissue_bam_fragments.qsub.sh`
  - Builds per-dataset/per-tissue BAM + fragments using total tissue barcodes.
  - Uses `analysis/src/alleloscope/build_tissue_files.sh` with:
    - `BUILD_BAMS=1`
    - `BUILD_FRAGMENTS=1`
    - no edge filtering in BAM/fragment inputs.

- `save_archr_tissue_no_edge.R`
  - Reads deepseq/lowseq ArchR projects.
  - Uses `*.no_edge_effect.barcodes.tsv` (or derives from all-edge).
  - Saves filtered ArchR projects to:
    - `Data/archr_tissue_no_edge/{dataset}_{tissue}_no_edge/`

- `run_save_archr_tissue_no_edge.qsub.sh`
  - qsub wrapper for `save_archr_tissue_no_edge.R`.

- `prepare_lowseq_alleloscope_tissue_from_existing.R`
  - Smart tissue correction for lowseq Alleloscope without full prep rerun.
  - Reuses existing `Data/alleloscope/{deepseq,lowseq}` matrices.
  - Builds tissue-specific lowseq inputs using no-edge barcodes and deepseq tissue SNP keys.
  - Writes outputs under:
    - `Data/alleloscope/lowseq_tissue_from_existing/{488B,489}/`
  - Writes deepseq tissue SNP VCFs to:
    - `Data/alleloscope/deepseq_tissue_snvs/`

- `run_prepare_lowseq_alleloscope_tissue_from_existing.qsub.sh`
  - qsub wrapper for lowseq tissue correction script.

## Typical run order

1. Build total-barcode tissue BAM/fragments:

```bash
qsub analysis/src/data_org/build_total_tissue_bam_fragments.qsub.sh
```

2. Save ArchR tissue no-edge projects:

```bash
qsub analysis/src/data_org/run_save_archr_tissue_no_edge.qsub.sh
```

3. Build lowseq tissue-specific Alleloscope inputs from existing outputs:

```bash
qsub analysis/src/data_org/run_prepare_lowseq_alleloscope_tissue_from_existing.qsub.sh
```
