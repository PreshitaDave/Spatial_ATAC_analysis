---
description: "Use when running or debugging NUMBAT in ATAC-bin mode for lowseq samples, especially pileup/phasing failures, 0-variant VCF errors, and resumable reruns. Enforces step-by-step checks, barcode-format validation, todo tracking, and smallest-scope reruns."
name: "NUMBAT Lowseq ATAC-Bin Repro Workflow"
applyTo: ["analysis/src/numbat/**/*.sh", "analysis/src/numbat/**/*.R", "analysis/qsub_logs/numbat_*", "Data/numbat/**"]
---
# NUMBAT Lowseq ATAC-Bin Repro Workflow

## Scope

- Primary goal: run lowseq NUMBAT in ATAC-bin mode reproducibly and compare-ready with deepseq outputs.
- Required behavior: diagnose failures with evidence, then rerun from the smallest safe checkpoint.

## Required Step Order

Run and verify in this order. Do not skip checks.

1. Verify scheduler/job status (`qstat`, `qacct`) and record job id.
2. Inspect stderr/stdout tail and identify exact failing step.
3. Validate input existence and non-emptiness:
   - BAM + BAI
   - barcode list
   - SNP VCF + index
   - phasing panel and genetic map
4. Validate barcode compatibility before pileup:
   - sample BAM `CB` tags
   - compare overlap with barcode list in raw and `-1` formats
   - choose format with higher overlap
   - hard-fail if overlap is 0 in both formats
5. Run pileup/phasing and immediately verify pileup output:
   - `cellSNP.base.vcf` exists
   - non-header variant row count > 0
   - AD/DP matrix dimensions have non-zero variant rows
6. Only then run `run_numbat_multiome.R` in ATAC-bin mode.
7. Generate postprocess plots and summarize outputs.

## Resume And Rerun Policy

- Reuse existing expensive artifacts when valid:
  - merged BAM
  - ATAC bin matrix
  - aggregated ATAC reference
- If pileup/phasing failed, clean only stale allele outputs for that dataset and rerun prep.
- Never delete successful NUMBAT result directories unless user asks.

## Failure Classification

Classify each failure as one of:

- environment/tooling (missing binary, module, PATH)
- input mismatch (barcode format mismatch, missing file, incompatible naming)
- runtime/algorithm (R error, phasing crash, memory/time limit)

For each failed run, report:

- failing step
- exact error line
- root cause hypothesis with evidence
- minimal corrective action

## Known High-Risk Failure: 0-Variant Pileup VCF

If you see:

- `Pileup VCF for sample ... has 0 variants`

Then treat barcode mismatch as first suspect:

- confirm BAM `CB` tags format (with or without `-1`)
- ensure pileup barcode file matches BAM format, not fragment format
- verify overlap count before rerun

## Logging And Tags

Use consistent log tags in wrappers and diagnostic notes:

- `[start]`, `[step]`, `[check]`, `[resume]`, `[warn]`, `[error]`, `[done]`

Each run summary must include:

- job id
- dataset
- chosen barcode format for pileup
- pileup variant row count
- whether run resumed or was fresh

## Todo Discipline

Maintain a concise single-level todo list for active recovery with statuses:

- `not-started`
- `in-progress`
- `completed`

After each pass, update only delta changes.

## Deliverables For Lowseq ATAC-bin Run

- exact rerun command(s)
- checkpoint inventory reused vs regenerated
- failure table for latest failed job
- current run status and next blocking step (if any)
