---
description: "Use when preparing or running Alleloscope lowseq/deepseq analyses for 488b tissue, comparing result similarity, debugging failed jobs, or building reproducible rerun workflows. Enforces checkpoint resume, optimized Est_regions usage, root-cause diagnostics, and structured todo logging."
name: "Alleloscope 488b Repro Workflow"
applyTo: ["analysis/src/alleloscope/**/*.R", "analysis/src/alleloscope/**/*.sh", "analysis/qsub_logs/**"]
---
# Alleloscope 488b Repro Workflow

## Scope

- Primary goal: compare how Alleloscope results look for 488b tissue in lowseq vs deepseq.
- Required outcome: reproducible preparation and analysis workflow that can restart from last valid checkpoint instead of re-running completed work.

## Dataset And Pairing Rules

- Keep lowseq and deepseq workflows structurally parallel (same ordered steps, same checkpoint policy, same logging detail).
- Use consistent sample naming in outputs: `lowseq` and `deepseq` only.
- When presenting comparisons, always include both:
  - qualitative similarity statement (for example, CNV pattern concordance)
  - quantitative summary (for example, number of cells, number of segments, overlap metrics)

## Reproducibility Requirements

- Every runnable script must:
  - print start/end timestamps
  - print key input paths and output directories
  - print core parameters (`NSLOTS`, threads, `cont`, filters)
- Save checkpoint RDS files after each major step with stable names in `output/rds/`.
- Before expensive steps, check if prior outputs exist and are complete; if yes, skip and log `[resume]`.
- Never delete successful checkpoint files during retries.

## Resume-From-Last-Best-Input Policy

- On rerun after failure, start from the last successfully written checkpoint.
- Prefer this sequence for restart decisions:
  1. existing step-level RDS checkpoint
  2. existing per-chromosome EM RDS files
  3. raw prepared matrices (`alt_all.mtx`, `ref_all.mtx`, `var_all.vcf`)
- If partial outputs are detected for a single step, clean only that step's incomplete artifacts and re-run only that step.

## Est_regions Parallelization Policy

- Prefer the optimized parallel Est_regions wrapper for chromosome/segment runs.
- For `cont=TRUE` runs, skip chromosomes that already have valid `chr*.rds` files.
- Validate outputs from parallel workers before downstream steps:
  - count expected vs completed chromosomes
  - capture and report failed chromosome IDs explicitly
  - abort downstream `Select_normal` if any required chromosome result is missing/corrupt

## Error Root-Cause Analysis Requirements

When any job fails, provide root cause and evidence, not only symptom text.

- Mandatory triage order:
  1. scheduler/accounting status (`qstat`, `qacct`)
  2. job stderr/stdout tail (or `.log`)
  3. input existence/format checks
  4. checkpoint integrity checks
- Classify failure into one of:
  - environment/tooling (module, PATH, binary missing)
  - input mismatch (barcode format, file missing, incompatible dimensions)
  - algorithm/runtime (worker crash, memory, NA/object type error)
- For each failed job, report:
  - failing step
  - exact error line/message
  - probable root cause
  - minimal corrective action

## Todo List And Action Plan Requirements

- Maintain a concise single-level todo list for recovery, with statuses:
  - `not-started`
  - `in-progress`
  - `completed`
- After each diagnostic pass, update the todo list with delta only (what changed).
- Recovery plans must prioritize smallest safe rerun scope first.

## Logging Requirements

- Each qsub wrapper must emit one dedicated run log under `analysis/qsub_logs/`.
- Include these fields at minimum:
  - job id, host, start time, end time, NSLOTS
  - script path and command invoked
  - key step boundaries (`Step 1`, `Step 2`, ...)
  - resume/skip decisions per chromosome or segment
- Use consistent tags in logs:
  - `[start]`, `[step]`, `[resume]`, `[skip]`, `[warn]`, `[error]`, `[done]`

## Guardrails For Comparison Readouts

- Do not claim lowseq/deepseq similarity without citing concrete output artifacts.
- If one dataset failed mid-pipeline, clearly mark comparison as provisional.
- Prefer explicit "what is comparable now" and "what is blocked" sections.

## Deepseq Barcode Format And VarTrix Requirements

The deepseq dataset has a **critical barcode mismatch** between what the somatic pipeline writes and what the BAMs contain:

- **BAM CB tags**: 8 bp barcodes (e.g. `TCTCGGAA`)
- **Somatic cell ID files** (`*.cell_snv.cellID.filter.csv`): 16 bp barcodes (e.g. `TCTCGGAATGAGTAGC`)
- **Consequence**: VarTrix receives 16 bp barcodes via `barcodes.tsv`, finds zero CB-tag matches in the BAM, and produces all-zero matrices. This was the confirmed root cause of job 4729913 failing on chr1 with message `Number of alignments skipped due to not being associated with a cell barcode: 33096718`.

### Barcode Mapping Files

Two files in `Data/variant_calling/` provide the 8 bp ↔ 16 bp mapping, keyed on a shared numeric `id` column:

| File | Barcode length | Key column |
|------|---------------|------------|
| `Data/variant_calling/deepseq_cell_data.csv` | 16 bp (full) | `id` |
| `Data/variant_calling/deepseq_cell_data_8bp.csv` | 8 bp (truncated) | `id` |

The 8 bp value is always the first 8 characters of the corresponding 16 bp barcode; `id` is the unique join key.

### Fix Applied In `alleloscope_prep_helpers.R` And `prepare_alleloscope_deepseq.R`

`prepare_alleloscope_inputs()` now accepts `vartrix_barcode_length` in config. When set to `8L`:

1. Truncates 16 bp barcodes to 8 bp (`substr(barcodes, 1L, 8L)`) for VarTrix input
2. Writes 8 bp barcodes to `barcodes_vartrix8bp.tsv` (passed to VarTrix `-c` flag)
3. `barcodes.tsv` always contains 16 bp barcodes (used by downstream Alleloscope analysis)
4. After all VarTrix chromosome calls, translates matrix column names back to 16 bp via `match()`

`prepare_alleloscope_deepseq.R` sets `vartrix_barcode_length = 8L` in config.

### VarTrix All-Zero Diagnostic Signature

If VarTrix log shows:
```
Number of alignments skipped due to not being associated with a cell barcode: <N>
```
where N equals all evaluated alignments → barcode length mismatch is the cause.
**Do not retry with the same `barcodes.tsv`.** Delete all stale `Data/alleloscope/deepseq/vartrix/chrN/` directories first, then fix and rerun.

### Lowseq vs Deepseq Barcode Summary

| Sample | BAM CB tag length | Somatic cell ID length | VarTrix barcode | Post-VarTrix translation |
|--------|-------------------|----------------------|-----------------|--------------------------|
| lowseq | 16 bp | 16 bp | 16 bp | none needed |
| deepseq | **8 bp** | 16 bp | **8 bp** | translate to 16 bp via `match()` |

Always delete all `Data/alleloscope/deepseq/vartrix/chrN/` directories before rerunning prep after a barcode fix, as all-zero cached files pass size checks but fail the `nnzero == 0` guard only when rows > 1000.

## Required Deliverables For 488b Comparison Tasks

- Reproducible rerun command list (prep + analysis for both lowseq/deepseq).
- Current checkpoint inventory (what exists and can be resumed).
- Failure table for latest jobs with root cause and fix.
- Next-step todo list to complete lowseq vs deepseq comparability.
