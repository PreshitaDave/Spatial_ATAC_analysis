# NUMBAT Multiome (ATAC-bin) Pipeline

This folder contains the working, end-to-end NUMBAT ATAC-bin pipeline used in this project.

## Scope

- Scripts: `analysis/src/numbat`
- Inputs and outputs: `Data/numbat`
- Cluster logs: `analysis/qsub_logs`

## Working scripts (authoritative)

- `numbat_common.sh`
	- Shared env/path helpers, logging, dry-run behavior, and required file checks.
- `prepare_numbat_atac_inputs.sh`
	- Creates all prerequisite inputs for NUMBAT:
	- ATAC bin matrix (`*_atac_bin.rds`)
	- Aggregated ATAC reference (`lambdas_ATAC_bincnt.rds`)
	- Allele count table from pileup + phasing (`*_atac_allele_counts.tsv.gz`)
	- Includes all autosomes (`chr1`-`chr22`) by default.
	- Uses tissue-specific barcodes from `Data/alleloscope/barcodes/<dataset>_<TISSUE>.barcodes.tsv` and appends `-1`.
	- Writes the required NUMBAT annotation file with exact columns `cell` and `group`.
- `run_numbat_atac_bin.sh`
	- Runs NUMBAT in ATAC-bin mode using `run_numbat_multiome.R`.
	- Writes runtime params to `run_numbat_params.rds` and post-processes outputs.
- `postprocess_numbat_results.R`
	- Produces vignette-style summary plots from NUMBAT output.
- `run_numbat_lowseq.qsub.sh`
	- Full lowseq cluster entrypoint (prep + run + plots).
- `run_numbat_deepseq.qsub.sh`
	- Full deepseq cluster entrypoint (prep + run + plots).

## Setup/utility script

- `setup_numbat_tools.sh`
	- Optional helper to clone/update NUMBAT repo and install NUMBAT R package.

## Required runtime dependencies

On SCC compute nodes, the qsub scripts load:

- `module load R`
- `module load samtools`

They also export PATH for:

- `cellsnp-lite`: `/projectnb/paxlab/presh/env/calicost_env/bin`
- `eagle`: `/projectnb/paxlab/presh/software/external/Eagle_v2.4.1`

NUMBAT bin scripts are resolved in this order:

1. `NUMBAT_BIN_DIR`
2. `${NUMBAT_REPO}/inst/bin`
3. `system.file('bin', package='numbat')`

Phasing resources used by `pileup_and_phase.R`:

- `PHASE_PANEL`
- `VCF_GENOME1K`
- `GMAP_GZ`

Defaults are auto-resolved in `prepare_numbat_atac_inputs.sh`, with fallback paths already wired for this workspace.

## Reproducible runbook

Run from project root:

```bash
cd /projectnb/paxlab/presh/projects/spatial_atac
```

Optional tool setup:

```bash
bash analysis/src/numbat/setup_numbat_tools.sh
```

Lowseq full run:

```bash
qsub analysis/qsub/pipeline/numbat/run_numbat_lowseq.qsub.sh
```

Deepseq full run:

```bash
qsub analysis/qsub/pipeline/numbat/run_numbat_deepseq.qsub.sh
```

Provenance in this workspace:

- Lowseq full wrapper submission: job `4700033`
- Deepseq full wrapper submission: job `4700034`

Monitor:

```bash
qstat -u preshita
tail -f analysis/qsub_logs/numbat_low.<JOBID>.out
tail -f analysis/qsub_logs/numbat_low.<JOBID>.err
tail -f analysis/qsub_logs/numbat_dep.<JOBID>.out
tail -f analysis/qsub_logs/numbat_dep.<JOBID>.err
```

## Critical implementation notes

- `aggregate_counts()` requires the annotation columns to be named exactly `cell` and `group`.
- The script uses direct `Rscript` calls for key stages to avoid shell-escaping/eval issues.
- Chromosomes default to all autosomes (`chr1`-`chr22`) via `CHROMS`.
- Existing outputs are reused when present to support restart/resume runs.

## Expected outputs

Per dataset (`lowseq` or `deepseq`):

- Input artifacts:
	- `Data/numbat/inputs/<dataset>_atac_bin.rds`
	- `Data/numbat/inputs/alleles/<dataset>_atac_allele_counts.tsv.gz`
- Final result directory:
	- `Data/numbat/results/<dataset>/atac_only`
- Plots directory:
	- `Data/numbat/results/<dataset>/atac_only/plots`

Expected plot files include:

- `phylogeny_heatmap.png`
- `single_cell_phylogeny.png`
- `mutation_history.png`
- `consensus_segments.png`
- `bulk_clone_profiles.png`
- `clone_profiles.png`
- `expression_roll.png`
- `phylogeny_cutree_grid.png`
