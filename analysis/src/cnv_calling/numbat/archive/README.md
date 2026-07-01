# Archive - Removed Scripts

**Date Cleaned**: June 8, 2026

This folder contains scripts that were removed from production because they are:
- Test/debug scripts (`.test.sh` suffix)
- Deprecated/experimental workflows
- Older versions that have been superseded

## test_scripts/ (9 scripts)

Test and debug scripts created during development:

- `check_params.qsub.sh` - Parameter validation test
- `test_params.sh` - Parameter testing
- `test_deepseq_*.sh` - Deepseq workflow tests
- `test_run_numbat.qsub.sh`, `test_run_numbat_call.R` - NUMBAT test runs
- `generate_*_references.test.sh` - Reference generation tests

**Action Taken**: Removed from production to keep scripts clean
**Reason**: These were intermediate development experiments

## old_versions/ (23 scripts)

Deprecated and alternative workflows:

### Reference Panel Experiments
- `*_refhca*.qsub.sh` (5 variants) - HCA reference panel experiments
- `*_refhca*.R` - HCA-based analysis scripts

### Older Script Versions
- `*_v2.qsub.sh` - Version 2 attempts (superseded by current scripts)
- `run_numbat_tissue_*.qsub.sh` - Early tissue-specific attempts
- `run_numbat_prep_tissue_*.qsub.sh` - Old prep workflows

### Combined Tissue Analyses
- `*_combined.qsub.sh` - Multi-tissue analyses (deprecated)
- `*_multiome_*.qsub.sh` - Early multiome integration attempts

### Alternative Workflows
- `run_numbat_lowseq_analysis_only.qsub.sh` - Analysis-only variant
- `run_numbat_lowseq_proper.qsub.sh` - Previous "proper" version
- `run_numbat_lowseq.qsub.sh` - Generic lowseq script
- `run_numbat_deepseq.qsub.sh` - Generic deepseq script
- `resume_phasing_lowseq.qsub.sh` - Resume script for stalled jobs

**Action Taken**: Moved to archive as reference only
**Reason**: These were experimental or superseded by simpler, tissue-specific scripts

## Should I Keep Experimenting?

If you need to test a new workflow:

1. Create a new script with `.test.sh` or `.debug.sh` suffix
2. Run tests on compute node (never on login node)
3. After testing is successful, move it to `production/` and rename (remove suffix)
4. Move any test outputs and logs to `archive/` with dated timestamps

## Restoration

To restore a script from archive:

```bash
# Find what you need
ls archive/*/

# Copy to appropriate location
cp archive/old_versions/run_numbat_tissue_deepseq.qsub.sh production/deepseq_488B/

# Update and test before using
```

## See Also
- `../ORGANIZATION.md` - Full folder structure
- `../production/` - Current production scripts
- `../lib/` - Reusable utilities
