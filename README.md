# Spatial_ATAC_analysis

Spatial ATAC-seq analysis for deepseq and lowseq cohorts (tissues 488B and 489).

## Repository Layout

- `analysis/`: analysis scripts, qsub launchers, outputs, and logs.
- `Data/`: input/intermediate data, tissue-level fragment files, and variant-calling outputs.
- `documentation/`: project documentation.
- `.github/`: copilot and workflow instructions.

## Key Notes

- Pipeline scripts are organized under `analysis/src/pipeline/` by stage number.
- Qsub launchers are consolidated under `analysis/qsub/` by domain.
- Backward-compatibility symlinks remain in `analysis/src/` so old commands still work.
- Tissue BAM/fragments are organized under `Data/fragments/{lowseq,deepseq}/tissue/`.
	Compatibility symlinks are preserved at `Data/variant_calling/{lowseq,deepseq}/tissue`.

See `analysis/README.md` for full script tree, run order, and qsub mapping.