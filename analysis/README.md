# Spatial_ATAC_analysis

Spatial ATAC-seq analysis pipeline for deepseq and lowseq (tissues 488B and 489).

## Current Organization

### Scripts (`analysis/src`)

- `analysis/src/pipeline/`: numbered end-to-end pipeline stages.
    - `1_Lib_install.Rmd`
    - `2_scATAC_giotto_obj_creation.Rmd`
    - `3_ArchR_ATAC_analysis.Rmd`
    - `4_CNV_calling_ATAC.Rmd`
    - `5_Xenium_giotto_analysis.Rmd`
    - `6_compare_atac_rna.Rmd`
    - `7_geseca_analysis_atac.Rmd`
    - `8_deepseq_variant_analysis.Rmd`
    - `8_variant_qc_comparison.R`
    - `9_somatic_snv_comparison.R`
    - `10_archr_variant_plotting_deepseq.R`
    - `10b_archr_variant_plotting_lowseq.R`
    - `11_comparing_tissue_variants.R`
    - `12_somatic_snv_characterization.R`
- Domain folders:
    - `analysis/src/alleloscope/`
    - `analysis/src/numbat/`
    - `analysis/src/data_org/`
    - `analysis/src/build_tissue/`
    - `analysis/src/pycistopic/`
    - `analysis/src/monopgen/`
- Backward-compatible symlinks are kept in `analysis/src/` for previously used script paths.

### Qsub Launchers (`analysis/qsub`)

All qsub scripts are centralized under `analysis/qsub/`:

- `analysis/qsub/pipeline/`
    - `run_11_comparing_tissue_variants.qsub.sh`
    - `run_12_somatic_char.qsub.sh`
- `analysis/qsub/alleloscope/`
    - deepseq, lowseq, and lowseq tissue-specific Alleloscope runs
- `analysis/qsub/data_org/`
    - tissue BAM/fragments build and lowseq Alleloscope tissue prep helpers
- `analysis/qsub/build_tissue/`
    - tissue barcode and edge/nFrags jobs
- `analysis/qsub/numbat/`
    - numbat launchers and resume scripts

Backward-compatible symlinks are retained at legacy locations in `analysis/src/**`.

### Data (`Data`)

- `Data/fragments/lowseq/tissue/`: tissue-level lowseq BAM/fragments/whitelists
- `Data/fragments/deepseq/tissue/`: tissue-level deepseq BAM/fragments/whitelists
- `Data/variant_calling/lowseq/tissue` and `Data/variant_calling/deepseq/tissue`: compatibility symlinks to the new `Data/fragments/*/tissue/` paths
- `Data/alleloscope/`: Alleloscope inputs/outputs

## Suggested Script Order

1. `analysis/src/pipeline/1_Lib_install.Rmd`
2. `analysis/src/pipeline/2_scATAC_giotto_obj_creation.Rmd`
3. `analysis/src/pipeline/3_ArchR_ATAC_analysis.Rmd`
4. `analysis/src/pipeline/4_CNV_calling_ATAC.Rmd`
5. `analysis/src/pipeline/5_Xenium_giotto_analysis.Rmd`
6. `analysis/src/pipeline/6_compare_atac_rna.Rmd`
7. `analysis/src/pipeline/7_geseca_analysis_atac.Rmd`
8. `analysis/src/pipeline/8_deepseq_variant_analysis.Rmd`
9. `analysis/src/pipeline/8_variant_qc_comparison.R`
10. `analysis/src/pipeline/9_somatic_snv_comparison.R`
11. `analysis/src/pipeline/10_archr_variant_plotting_deepseq.R`
12. `analysis/src/pipeline/10b_archr_variant_plotting_lowseq.R`
13. `analysis/src/pipeline/11_comparing_tissue_variants.R`
14. `analysis/src/pipeline/12_somatic_snv_characterization.R`

## HPC Guidance

- Submit full jobs through qsub launchers in `analysis/qsub/`.
- Use login nodes only for quick checks (small commands, no full pipelines).
- Inspect `analysis/qsub_logs/*.out` and `analysis/qsub_logs/*.err` routinely.

## Contact

- Project owner: `preshita@bu.edu`
