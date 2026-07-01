Spatial ATAC-seq Analysis Instructions

Project Path: /projectnb/paxlab/presh/projects/spatial_atac
Data: Spatial ATAC-seq (deepseq vs. lowseq) for patients 448B and 489.
Structure: Scripts are located in analysis/src. Key subfolders include Data (organization)

Coding & Documentation Standards

Documentation: Maintain a master README and individual READMEs for each analysis step.
Best Practices: Use Git version control, descriptive file naming, and thorough code commenting.
HPC (qsub) Guidelines
**CRITICAL PATH CHECKS BEFORE JOB SUBMISSION:**
- **ALWAYS verify ALL input paths exist and are readable** - check fragments, barcodes, BAM files, reference files
- **ALWAYS use absolute paths** (never relative paths in production scripts)
- **ALWAYS check output directory exists and is writable** before qsub submission
- **ALWAYS ensure output directory is EMPTY/CLEAN before VarTrix jobs** - VarTrix will fail if output dir contains existing files from previous runs. Clean with: `rm -rf {output_dir}/* && mkdir -p {output_dir}`
- **ALWAYS verify function parameters match invocations** - Check that all R/Python function parameters match how they're called in scripts BEFORE executing/submitting. Mismatched parameters cause silent failures or incorrect outputs.
- **Organize input files** logically in `Data/01_inputs/` with proper subfolders (fragments/, barcodes/, references/)
- **Organize output files** in appropriate `Data/0X_outputs/` folders with descriptive naming
- **Use consistent naming conventions** across all scripts and data folders
- Move anything old and unused to `Data/01_inputs/archive/` with dated subfolder names

**Script Placement & Organization (Mandatory)**:
- **NEVER create or edit scripts in root directory or main `analysis/` folder**
- **All qsub scripts**: `analysis/src/cnv_calling/{tool}/{script_name}.qsub.sh`
- **All helper scripts**: Colocated with qsub script in same subdirectory
- **All R/Python analysis**: `analysis/src/{analysis_type}/{script_name}.R` or `.py`
- **After creating/editing ANY script**: Verify it's in the correct location
- **If script is in wrong location**: Move it immediately, do NOT leave duplicates
- **RULE**: Each production script has ONE canonical location, never multiple copies

**Testing & Debugging Policy**:
- When creating test/debug scripts during development, use `.test.sh` or `.debug.sh` suffix for naming
- After testing is complete, IMMEDIATELY move test scripts and all test-related files to archive:
  - Location: `Data/01_inputs/archive/test_files_archive_YYYYMMDD/`
  - Include all: test scripts, test logs, error files, temporary outputs
- Never leave test files in main analysis/src, root directory, or primary qsub_logs folders
- This keeps the project clean and prevents confusion between production and test code

**CRITICAL - Automatic File Cleanup After Job Completion**:
- **YOU MUST ENFORCE THIS AUTOMATICALLY - Do NOT require user to ask for cleanup**
- **After ANY job completes successfully (check qstat or log files), IMMEDIATELY:**
  1. Move completed job logs: `mv analysis/qsub_logs/{jobname}.{JOBID}.{err,out} analysis/qsub_logs/archive/completed_jobs_YYYYMMDD/`
  2. Move test files: `mv {name}.test.* Data/01_inputs/archive/test_files_archive_YYYYMMDD/`
  3. Move temporary files: `rm` any `.tmp`, `.debug.log`, or temporary outputs not in approved folders
  4. Verify output files are in correct locations per File Organization Standards
  5. **Report** what was cleaned up and moved, don't make user ask
- **Cleanup Pattern** (after confirming job success):
  ```bash
  mkdir -p analysis/qsub_logs/archive/completed_jobs_$(date +%Y%m%d)
  mkdir -p Data/01_inputs/archive/test_files_archive_$(date +%Y%m%d)
  # Move job logs
  mv analysis/qsub_logs/*.$JOBID.err analysis/qsub_logs/archive/completed_jobs_$(date +%Y%m%d)/ 2>/dev/null || true
  mv analysis/qsub_logs/*.$JOBID.out analysis/qsub_logs/archive/completed_jobs_$(date +%Y%m%d)/ 2>/dev/null || true
  # Move test files
  mv ./*.test.* Data/01_inputs/archive/test_files_archive_$(date +%Y%m%d)/ 2>/dev/null || true
  ```
- **What STAYS in main directories**: Only production/active files, never test or old files
- **Rule**: If user doesn't see cleanup happening, you're not doing your job - this is non-negotiable

**Module & R Script Initialization Fixes**:
- **Problem**: `module: command not found` or `Rscript: command not found` in SGE jobs
- **Solution**: ALL qsub scripts MUST include proper module initialization BEFORE using module or Rscript commands
- **Required Pattern**:
  ```bash
  #!/bin/bash
  set -eo pipefail
  
  # CRITICAL: Initialize module system before using 'module' command
  set +u
  for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
    if [[ -f "$profile_file" ]]; then
      . "$profile_file" 2>/dev/null || true
      break
    fi
  done
  set -u
  
  # NOW module command is available
  module load R
  which Rscript  # Verify it worked
  ```
- **Why this is needed**: SGE job scripts run in minimal shell environments where `module` command is not in PATH
- **Always verify**: After `module load R`, run `which Rscript && Rscript --version` to confirm module loaded correctly

**Progress Indicators in Long-Running Jobs (CRITICAL)**:
- **Purpose**: Enable real-time monitoring and early error detection without waiting for job completion
- **Required for**: Jobs > 1 hour, VarTrix runs, matrix operations, any I/O intensive tasks
- **Implementation in Bash**:
  ```bash
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] Starting task description" >&2
  # ... task code ...
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] Task completed: result summary" >&2
  ```
- **Implementation in R** (CRITICAL for long operations):
  ```R
  message(sprintf("[%s] Step N: Description starting", Sys.time()))
  # ... R code ...
  message(sprintf("[%s] Step N: Completed with N items processed", Sys.time()))
  ```
- **Checklist for adding progress indicators**:
  1. Add timestamped message at START of each major step
  2. Add timestamped message at END of each major step with result summary
  3. For loops/iterations > 100: Add progress every 10-20% 
  4. For matrix operations: Log dimensions before/after
  5. For file I/O: Log file size and paths being created
  6. Always use `message()` in R (not `print()`) - messages go to stderr, visible in logs
  7. Always use `>&2` in bash - redirects to stderr for immediate visibility
- **Monitoring with progress indicators**:
  ```bash
  # Real-time log monitoring while job runs
  tail -f analysis/qsub_logs/jobname.$JOB_ID.err
  
  # Check progress at any time
  tail -20 analysis/qsub_logs/jobname.$JOB_ID.err
  ```
- **Benefits**:
  - Detect stuck/hung jobs within minutes, not hours
  - Identify performance bottlenecks (which step is slow?)
  - Validate that jobs are processing expected data volumes
  - Enable early termination if job is going wrong

**Path Validation & Symlink Resolution**:
- **BEFORE every script run**: Implement path checking to handle broken/relocated symlinks
- **Required validation pattern** (Bash):
  ```bash
  # Function to check and resolve paths
  check_and_resolve_path() {
    local requested_path="$1"
    local lookup_file="$2"  # Path to symlink_mapping_AFTER.txt or similar
    
    if [[ -L "$requested_path" && ! -e "$requested_path" ]]; then
      # Symlink is broken - look up correct path
      echo "[WARN] Broken symlink detected: $requested_path" >&2
      
      if [[ -f "$lookup_file" ]]; then
        local corrected=$(grep "^$requested_path" "$lookup_file" | cut -d'|' -f2)
        if [[ -n "$corrected" && -e "$corrected" ]]; then
          echo "[INFO] Found corrected path: $corrected" >&2
          echo "$corrected"
          return 0
        fi
      fi
      
      echo "[ERROR] Could not resolve broken symlink: $requested_path" >&2
      return 1
    elif [[ -e "$requested_path" ]]; then
      echo "$requested_path"
      return 0
    else
      echo "[ERROR] Path does not exist: $requested_path" >&2
      return 1
    fi
  }
  
  # Usage in scripts:
  # EMBEDDINGS_PATH=$(check_and_resolve_path "Data/Embeddings" "Data/symlink_mapping_AFTER.txt") || exit 1
  ```
- **Required validation pattern** (R):
  ```R
  # Function to check and resolve paths
  check_and_resolve_path <- function(requested_path, lookup_file = NULL) {
    # First check if path exists and is valid
    if (file.exists(requested_path)) {
      return(requested_path)
    }
    
    # If it's a symlink, check if target exists
    if (Sys.readlink(requested_path) != "") {
      target <- Sys.readlink(requested_path)
      if (file.exists(target)) {
        return(target)
      }
      # Symlink is broken - try lookup file
      if (!is.null(lookup_file) && file.exists(lookup_file)) {
        mappings <- read.table(lookup_file, sep="|", stringsAsFactors=FALSE)
        corrected <- mappings[mappings[,1] == requested_path, 2]
        if (length(corrected) > 0 && file.exists(corrected)) {
          cat(sprintf("[INFO] Found corrected path: %s\n", corrected))
          return(corrected)
        }
      }
      stop(sprintf("[ERROR] Broken symlink with no mapping: %s", requested_path))
    }
    
    stop(sprintf("[ERROR] Path does not exist: %s", requested_path))
  }
  
  # Usage in scripts:
  # embeddings_path <- check_and_resolve_path("Data/Embeddings", "Data/symlink_mapping_AFTER.txt")
  ```
- **Lookup file format** (`symlink_mapping_AFTER.txt`):
  ```
  Data/Embeddings|Data/04_analysis/03_intermediate/archr/artifacts/Embeddings
  Data/variant_calling|Data/04_analysis/03_intermediate/variant_calling/monopogen
  ```
- **When paths change**: 
  1. Verify all symlinks are still valid: `find Data -type l -exec test ! -e {} \; -print`
  2. If broken symlinks found, check mapping files for correct location
  3. Update symlink targets or use correction function in scripts
  4. Update mapping files after changes

Symlink Management Policy:
- **ALWAYS label symlinks with ".lnk" suffix** so they're distinguishable from actual files (e.g., `reference.bam.lnk`)
- **Internal symlinks** (pointing within /projectnb/paxlab/presh/): Remove and replace with absolute paths in scripts
- **External symlinks** (pointing outside /projectnb/paxlab/presh/): Keep in designated location `Data/01_inputs/archive/symlink_external/`
- **Before organizing files**: Create `symlink_mapping_BEFORE.txt` documenting all symlinks with targets
- **After organizing files**: Create `symlink_mapping_AFTER.txt` and verify all paths are correct
- **If scripts break after path changes**: Reference mapping files to diagnose missing paths and update affected scripts
- **VERIFICATION**: After removing symlinks, ensure all scripts have been updated to use absolute paths directly

**Update ORGANIZATION_SUMMARY.md** with new file locations and symlink policies for future reference.

Efficiency: Ensure $\ge$80% utilization of requested CPU cores.

File Organization Standards

Arrow Files:
- Location: `Data/01_inputs/arrow/`
- Format: `.arrow` (ArchR Arrow format for ATAC-seq data)
- Naming convention: `{sample_name}.arrow` (e.g., `Deepseq_488B.arrow`)
- Created from fragment BED files via createArrowFiles()
- Organized in single folder for easy reference

Fragment Files:
- Location: `Data/01_inputs/fragments/{object}/`
- Format: `.bed.gz` (compressed BED files from barcode-filtered fragments)
- Each object should have filtered fragments in its own subfolder
- Index files (`.tbi`) should accompany compressed fragments

BAM Files:
- Location: `Data/01_inputs/bam/` with `.bam.lnk` symlinks pointing to archived BAM files
- Original BAM files preserved in: `Data/01_inputs/archive/fragments_old_nested_source/{depth}/tissue/`
- Naming convention: `{object}.bam.lnk` (e.g., `deepseq_489.bam.lnk`)
- BAM files contain all barcodes before edge-effect filtering

Barcode Files:
- Location: `Data/01_inputs/barcodes/tissue_barcodes/{object}/`
- Files per object:
  - `{object}.barcodes.tsv`: All barcodes before filtering
  - `{object}.no_edge_effect.barcodes.tsv`: Barcodes kept after edge-effect removal
  - `{object}.edge_effect.barcodes.tsv`: Barcodes removed as edge effects
  - `{object}_nFrags_from_fragments.tsv.gz`: Fragment counts per barcode

Cleanup Rules:
- Remove fragment files from any folder except `01_inputs/fragments/`
- Move stale/archive files to `01_inputs/archive/` with dated subfolder names
- External data references (outside paxlab/presh) should use symlinks with `.lnk` suffix
- Archive external symlinks in `01_inputs/archive/symlink_external/`
- **Test Scripts & Files**: After debugging/testing is complete, ALWAYS move test scripts and log files to archive folders:
  - Test scripts (`.sh`, `.R`, etc.) → `Data/01_inputs/archive/test_files_archive_YYYYMMDD/`
  - Test log files → `analysis/qsub_logs/archive/`
  - Testing outputs/intermediates → appropriate archive folders with dated subfolder names
  - **DO NOT** leave test scripts in main analysis/src, root project directory, or main qsub_logs folders
  - Use descriptive names for archive folders (e.g., `test_files_archive_20260513`, `debugging_logs_20260513`)

ArchR Output Structure:
- Arrow Files: `Data/01_inputs/arrow/{sample_name}.arrow` (centralized storage)
- ArchR Projects: `Data/01_outputs/archR_objects/{object}/{object}_archR_project_final/`
- QC Plots: `analysis/plots/archr_obj/archR_qc_{object}.pdf` (4 publication-ready PDFs)
- Processing Summary: `analysis/plots/archr_obj/archR_processing_summary.tsv`

NUMBAT Reference Files:
- Location: `Data/04_analysis/cnv/numbat/reference/`
- Key files:
  - `lambdas_ATAC_bincnt.rds` - Aggregated ATAC reference for NUMBAT ATAC-bin mode (178K)
  - `var220kb.rds` - Genomic bins (220kb windows) for binning (80K)
  - `phased_panel_bcf_links/` - Eagle phasing panel directory
  - `par_numbatm.rds` - NUMBAT analysis parameters
- Note: Lambda file is generated from normal sample ATAC data using `get_binned_atac.R --generateAggRef`

Reorganized File Structure (as of May 14, 2026):
- `Data/01_inputs/` - Input files (fragments, barcodes, references, BAM)
- `Data/01_outputs/` - Initial outputs (archR_objects with RDS files and projects)
- `Data/02_references/` - Genome references and annotations
- `Data/04_analysis/` - **CONSOLIDATED analysis folder**
  - `04_analysis/03_intermediate/` - Intermediate outputs (moved from root)
    - `archr/` - ArchR artifacts, metadata, projects
    - `variant_calling/` - Variant calling results
  - `04_analysis/05_results/` - Final results (moved from root)
  - `04_analysis/cnv/` - CNV analysis (alleloscope, NUMBAT)
  - `04_analysis/multiomic/` - Multi-omic integration (cellwalkr, pycistopic)
- `04_analysis/` symlinks access intermediate/results folders (automatic via symlinks in Data/)
- `06_logs/` - **REMOVED** (logs now in analysis/qsub_logs/)
- Symlinks in Data/ folder point to relocated folders (e.g., `Data/Embeddings` → `Data/04_analysis/03_intermediate/archr/artifacts/Embeddings`)

Interactive Compute Session Policy:
- For any code testing, diagnostics, or script execution, do not stay on the SCC login node.
- **ALWAYS verify hostname before running code**: Run `hostname` - if it starts with "scc1", you're on LOGIN NODE (danger zone).
- **ALWAYS open a visible terminal** for the user to see real-time command execution and output.
- **ALWAYS execute `qrsh` inside a `tmux` session** to ensure persistence and visibility.
- **CRITICAL: Check for existing tmux sessions before creating new ones**: 
  - BEFORE creating a tmux session with `tmux new -As {name}`, FIRST check if one already exists
  - Command to list sessions: `tmux list-sessions` 
  - If a session exists (e.g., "spatial_atac_work"), attach to it: `tmux a -t spatial_atac_work`
  - ONLY create a new session if one doesn't already exist
  - This prevents creating duplicate sessions and losing active compute node allocations
  - **DO NOT** create a new session if user already has an active one with compute resources allocated

**COMPLETE TESTING & SUBMISSION WORKFLOW:**
1. **STEP 1: Verify Compute Node Allocation**
   - Run: `hostname` → should NOT start with "scc1"
   - If on login node, run: `qrsh -l h_rt=16:00:00 -pe omp 1 -P paxlab -l mem_per_core=8G`
   - Inside tmux: `tmux new -As spatial_atac_work` (or attach to existing: `tmux a -t spatial_atac_work`)

2. **STEP 2: Load Required Modules & Verify**
   - Run the module initialization block (see Module & R Script Initialization Fixes above)
   - Verify with: `module list` and `which Rscript && Rscript --version`
   - Verify input files: `ls -lh Data/01_inputs/fragments/*/` 
   - Verify output directory is writable: `mkdir -p Data/0X_outputs/{name}` and `touch Data/0X_outputs/{name}/.test`

3. **STEP 3: Create & Run Small Test Script**
   - Create test script with `.test.sh` or `.debug.sh` suffix
   - Test script should be SMALL and FAST (not full production run)
   - Include diagnostic output: `echo "[TEST STEP N] Description"` at each checkpoint
   - Test ONLY critical functionality: module loading, file access, basic R parsing, key computations

4. **STEP 4: Verify Test Success**
   - Confirm test ran WITHOUT ERRORS
   - Verify output files were created where expected
   - If test FAILED → debug, fix script, re-test (do NOT proceed to step 5)
   - If test PASSED → proceed to step 5

5. **STEP 5: Archive Test Files Immediately**
   - Move test script: `mv {name}.test.sh Data/01_inputs/archive/test_files_archive_YYYYMMDD/`
   - Move test logs: `mv {name}.test.*.log Data/01_inputs/archive/test_files_archive_YYYYMMDD/` (or to `analysis/qsub_logs/archive/`)
   - Move test outputs: `mv {name}.test.* Data/01_inputs/archive/test_files_archive_YYYYMMDD/`
   - **NEVER leave test files in main analysis/src, root project directory, or main qsub_logs/**

6. **STEP 6: Create Production qsub Script**
   - Copy test script patterns → production qsub wrapper script
   - Replace `.test.sh` → `.qsub.sh` naming
   - Include FULL module initialization block
   - Include ALL critical path checks
   - Location: appropriate folder under `analysis/qsub/` or `analysis/src/*/` 
   - Use absolute paths EVERYWHERE

7. **STEP 7: Final Verification Before qsub Submission**
   - Verify script syntax: `bash -n {script}.qsub.sh` (should have no output = no errors)
   - Verify all input paths: `ls -lh {path}` for every input referenced
   - Verify output directory exists: `mkdir -p {output_dir}` 
   - Verify script has correct shebang: `head -1 {script}.qsub.sh` → should be `#!/bin/bash`
   - **CRITICAL: Check that all function parameters match the command invoked**
     - For R scripts: Verify that `Rscript script.R arg1 arg2` matches `commandArgs(trailingOnly=TRUE)` parsing and function definitions
     - For functions: Ensure function parameters in script match how they're called in the qsub wrapper
     - Example: If calling `run_analysis(input_file, output_dir, cores=8)`, verify the qsub script passes all required args
     - Do this check BEFORE submitting or executing anything - mismatched parameters cause silent failures or wrong outputs

8. **STEP 8: Submit Production Job via qsub**
   - Run: `qsub {script}.qsub.sh`
   - Note the Job ID returned
   - Check status: `qstat -j {JobID}`
   - Monitor log file in real-time: `tail -f analysis/qsub_logs/{folder}/{name}_{JobID}.log`

Module Loading Critical:
- R requires: module load R (do NOT assume R is in PATH)
- Check available modules: module avail
- After loading, verify: which Rscript

Key Commands to Always Execute:
- [Login node check] hostname
- [Module check] module list
- [R check] which Rscript && Rscript --version
- [Input check] ls -lh Data/01_inputs/fragments/*/ (verify files exist and are readable)
- [Job check] qstat -u preshita (check job queue)

Terminal Visibility Rule:
- ALWAYS use send_to_terminal with waitForOutput=true to show commands
- ALWAYS include descriptive echo statements: echo "[STEP N] Description: command"
- NEVER run code silently in background
- NEVER assume code works without showing output

CRITICAL - NEVER TEST ON LOGIN NODE (scc1):
- **ALL testing, debugging, and script execution MUST be on compute node only**
- If you run commands on scc1 (login node), you're doing it WRONG
- ALWAYS test on compute node allocated via qrsh inside tmux session
- If you need to test: use the compute node terminal from tmux session
- If you're on scc1 and want to test: attach to tmux spatial_atac_work and verify hostname != scc1*
- This is non-negotiable for reproducibility and resource management
CNV Analysis Workflows (NUMBAT & Alleloscope)

## NUMBAT Analysis Pipeline

**Overview**: NUMBAT performs CNV detection from scATAC-seq data using a generative model that combines:
1. Copy-number agnostic variant calling (pileup stage)
2. Long-range phasing (genetic map + reference panel)
3. ATAC-seq binning (220kb genomic windows)
4. Aggregated reference comparison (ATAC-based CNV baseline)
5. CNV inference (multi-omic integration if available)

**Complete NUMBAT Workflow Steps**:

### Step 1: Prepare Allele Counts & Variant Calls (pileup_and_phase.R)
```bash
# Input: BAM file (all barcodes), barcode list, reference genome
# Output: alleles.csv (variant counts per barcode), pileup.csv

# Command format:
qsub analysis/src/cnv_calling/numbat/run_numbat_pileup_phase_{TISSUE}.qsub.sh

# Critical checks BEFORE submission:
- Verify BAM file path exists and is readable: ls -lh Data/01_inputs/bam/
- Verify barcode file exists: ls Data/01_inputs/barcodes/tissue_barcodes/{TISSUE}/
- Verify reference genome files: ls Data/02_references/
- Verify output directory is writable: mkdir -p Data/04_analysis/cnv/numbat/results/{TISSUE}
```

### Step 2: Generate ATAC-Binned Matrix & Reference (get_binned_atac_fixed.R)
```bash
# Input: Fragment BED.GZ file, barcode list, genomic bins (220kb)
# Output: 
#   - adata_atac.rds (cell x bin matrix)
#   - lambdas_ATAC_bincnt.rds (aggregated reference for CNV detection)

# Key: ALWAYS use get_binned_atac_fixed.R (NOT original get_binned_atac.R)
# Reason: Handles barcode format mismatch (see CRITICAL ISSUE section below)

# For reference generation with --generateAggRef flag:
qrsh -l h_rt=02:00:00 -pe omp 8 -P paxlab -l mem_per_core=8G bash \
  analysis/src/cnv_calling/numbat/generate_all_tissue_references.test.sh

# Within qrsh/tmux session: Run the reference generation script
# This generates lambdas_*_ATAC_bincnt.rds for all tissues
# Output location: Data/04_analysis/cnv/numbat/reference/
```

### Step 3: Run NUMBAT CNV Analysis (run_numbat_multiome.R)
```bash
# Input: 
#   - alleles.csv (from step 1)
#   - adata_atac.rds (from step 2)
#   - lambdas_ATAC_bincnt.rds (from step 2)
#   - Reference parameters (par_numbatm.rds)

# Output: CNV calls, cell-level copy-number profiles, tumor vs normal classification

qsub analysis/src/cnv_calling/numbat/run_numbat_analysis_{TISSUE}.qsub.sh
```

**Reference File Requirements**:
- Location: `Data/04_analysis/cnv/numbat/reference/`
- Required files:
  - `lambdas_ATAC_bincnt.rds` - Aggregated ATAC reference (178K) - **USE TISSUE-SPECIFIC VERSION**
  - `var220kb.rds` - Genomic bins (80K)
  - `par_numbatm.rds` - NUMBAT parameters (46 bytes)
  - `phased_panel_bcf_links/` - Eagle phasing panel directory
  - `genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz` - 1000G SNPs (140M)
  - `genetic_map_hg38_withX.txt.gz` - Genetic map for phasing (512 bytes)

**Key Parameters for qrsh Reference Generation**:
```bash
# Recommended resource allocation:
qrsh -l h_rt=02:00:00 -pe omp 8 -P paxlab -l mem_per_core=8G

# Within tmux session (for persistence):
tmux new-session -s numbat_refs
# Then run: bash analysis/src/cnv_calling/numbat/generate_all_tissue_references.test.sh
```

---

## Alleloscope Analysis Pipeline

**Overview**: Alleloscope infers single-cell haplotypes from aggregated variant data. Useful for:
- Validating CNV calls from other methods
- Phase-aware CNV analysis
- Tumor clonal structure inference

**Complete Alleloscope Workflow Steps**:

### Step 1: Prepare Variant Input Data
```bash
# Input: BAM file (tumor tissue), barcode list, reference genome
# Output: Variant count matrix suitable for Alleloscope

qsub analysis/qsub/alleloscope/{deepseq|lowseq}/prepare_alleloscope_{TISSUE}.qsub.sh
```

### Step 2: Run Alleloscope Haplotype Inference
```bash
# Input: Variant count matrix from step 1
# Output: Single-cell haplotype assignments, phase probabilities

qsub analysis/qsub/alleloscope/{deepseq|lowseq}/run_alleloscope_{TISSUE}.qsub.sh
```

### Step 3: Integration with NUMBAT (Optional)
- Use Alleloscope haplotype calls to validate/refine NUMBAT CNV calls
- Cross-validate tumor vs normal classification
- Compare clonal CNV patterns

---

## CRITICAL ISSUE: Barcode Format Compatibility

**The Problem**:
- 10X Cell Ranger pipeline adds `-1` suffix to fragment barcodes (e.g., `TGGCTTCAAGCCATGC-1`)
- Spatial barcode files contain raw barcodes without suffix (e.g., `AACGTGATTCCGTCTA`)
- Without suffix handling: barcode matching returns 0 fragments → corrupted output
- Example failure: binned ATAC matrix = 182 bytes (should be ~100M)

**Root Cause Analysis**:
```R
# Fragment barcodes (from fragments.bed.gz):
# TGGCTTCAAGCCATGC-1  chr1  1000  2000  1
# AACGTGATTCCGTCTA-1  chr1  2000  3000  1

# Cell barcode file (spatial_barcodes.tsv):
# TGGCTTCAAGCCATGC
# AACGTGATTCCGTCTA

# Matching without fix: 0 matches (suffix mismatch)
# Matching with fix: 2 matches (suffix stripped)
```

**The Fix**:
Use patched script: `analysis/src/cnv_calling/numbat/get_binned_atac_fixed.R`
```R
# Key changes:
# 1. Strip suffix before matching
barcodes_clean <- sub("-[12]$", "", barcodes_from_fragments)
# Now matches work correctly

# 2. When using --generateAggRef, add group column to annotations
annot <- data.table::fread(barcode_file)
colnames(annot) <- "barcode"
annot$group <- "reference"  # Required by aggregate_counts()
# This fixes: Error object 'group' not found in aggregate_counts()
```

**Prevention Checklist**:
1. **Always verify barcode compatibility** before running binning:
   ```bash
   # Check fragment barcodes:
   zcat Data/01_inputs/fragments/{TISSUE}/fragments.bed.gz | head -5 | cut -f4 | sort -u
   # Output should show: {16bp}-1
   
   # Check cell barcode file:
   head Data/01_inputs/barcodes/tissue_barcodes/{TISSUE}/{TISSUE}.barcodes.tsv | cut -f1 | sort -u
   # Output should show: {16bp} (no suffix)
   ```

2. **Always use get_binned_atac_fixed.R** for NUMBAT binning (never original)

3. **Verify output size** after binning:
   ```bash
   # Check ATAC matrix:
   ls -lh Data/04_analysis/cnv/numbat/results/{TISSUE}/adata_atac.rds
   # Should be ~100-300M, NOT 182 bytes
   ```

4. **Check fragment matching logs**:
   ```bash
   tail Data/04_analysis/cnv/numbat/reference/{TISSUE}_reference_generation.log
   # Should show: "N fragments matched" (not 0)
   ```

---

## CRITICAL ISSUE 2: Barcode File Consistency Between Pipeline Stages

**The Problem**:
- NUMBAT requires cell-by-cell correspondence between ATAC bin matrix and allele counts
- If different barcode files are used in each stage, inputs will have mismatched cell sets
- Causes: "No matching cell names between count_mat and df_allele" → silent failure
- Result: "Filtering out all X cells with 0 coverage" → no CNV calls

**Root Cause Example**:
```
Stage 1 (pileup_and_phase.R):
  Uses: lowseq_489.barcodes.tsv (4,671 cells - FULL set)
  Creates: allele_counts with 4,671 cells

Stage 2 (get_binned_atac_fixed.R):
  Uses: lowseq_489.no_edge_effect.barcodes.tsv (4,211 cells - FILTERED)
  Creates: ATAC matrix with 4,211 cells

Result when NUMBAT runs:
  ✗ Mismatch: 4,211 cells in ATAC, 4,671 in allele counts
  ✗ 460 allele-only cells have 0 ATAC coverage
  ✗ NUMBAT filters them all out
  ✗ Analysis fails with no results
```

**Prevention & Fix**:

1. **Always validate before running NUMBAT analysis**:
   ```bash
   # Check barcode consistency
   Rscript analysis/src/cnv_calling/numbat/validate_numbat_inputs.R {TISSUE}
   
   # Example for lowseq_489:
   Rscript analysis/src/cnv_calling/numbat/validate_numbat_inputs.R lowseq_489
   ```

2. **If validation FAILS**, regenerate ATAC matrix with correct barcode file:
   ```bash
   # Regenerate using SAME barcode file as pileup stage
   qsub analysis/src/cnv_calling/numbat/regenerate_atac_with_correct_barcodes.qsub.sh
   
   # This script:
   # - Loads the pileup barcode file (used for allele counts)
   # - Regenerates ATAC matrix with that file
   # - Ensures 100% barcode match between inputs
   # - Creates backup of original matrix
   ```

3. **Built-in validation in NUMBAT scripts**:
   - All NUMBAT analysis scripts now include automatic validation as Step 0
   - Scripts check consistency BEFORE running analysis
   - If validation fails, analysis won't start (prevents wasted compute time)
   - Location: Lines checking `Rscript analysis/src/cnv_calling/numbat/validate_numbat_inputs.R`

4. **Using the validation as a general check across tissues**:
   ```bash
   # Validate ANY tissue before analysis
   Rscript analysis/src/cnv_calling/numbat/validate_numbat_inputs.R lowseq_489
   Rscript analysis/src/cnv_calling/numbat/validate_numbat_inputs.R lowseq_488B
   Rscript analysis/src/cnv_calling/numbat/validate_numbat_inputs.R deepseq_489
   
   # All uses same validation script - works for any tissue structure
   # Output shows:
   #   - ATAC cells count
   #   - Allele cells count
   #   - Overlap percentage
   #   - Cells-only-in-ATAC count
   #   - Cells-only-in-allele count
   #   - PASS or FAIL verdict
   ```

---

## Important Workflow Points

### For Both NUMBAT and Alleloscope:
- ✓ Always test on **compute node** (qrsh) NOT login node (scc1*)
  - Login node is FOR FILE INSPECTION ONLY
  - All code execution must be on compute allocation
  - Use tmux for persistence: `tmux new -As spatial_atac_work`

- ✓ Always use **absolute paths** (no relative paths in production scripts)

- ✓ Always verify **input file existence** before job submission:
  ```bash
  # Fragment files
  ls -lh Data/01_inputs/fragments/*/
  
  # Barcode files
  ls -lh Data/01_inputs/barcodes/tissue_barcodes/*/
  
  # Reference files (for NUMBAT)
  ls -lh Data/04_analysis/cnv/numbat/reference/
  ```

- ✓ Always create **output directories** and verify writability:
  ```bash
  mkdir -p Data/04_analysis/cnv/numbat/results/{TISSUE}
  touch Data/04_analysis/cnv/numbat/results/{TISSUE}/.test && rm .test
  ```

- ✓ Always include **module initialization** in qsub scripts:
  ```bash
  # Initialize module system
  set +u
  for profile_file in /etc/profile /etc/profile.d/modules.sh /usr/share/modules/init/bash; do
    if [[ -f "$profile_file" ]]; then
      . "$profile_file" 2>/dev/null || true
      break
    fi
  done
  set -u
  
  module load R
  which Rscript  # Verify module loaded
  ```

### For NUMBAT Specifically:
- Generate **separate reference matrices per tissue** using --generateAggRef flag
  - Each tissue needs its own lambdas_ATAC_bincnt.rds
  - Use qrsh compute allocation: `qrsh -l h_rt=02:00:00 -pe omp 8 ...`
  - Reference generation takes ~10-30 min per tissue

- Use **fixed binning script** to avoid barcode mismatch errors
  - Location: `analysis/src/cnv_calling/numbat/get_binned_atac_fixed.R`
  - Includes diagnostic output showing fragment match counts
  - Logs saved to: `Data/04_analysis/cnv/numbat/reference/{TISSUE}_reference_generation.log`

- Validate **allele counts file** after pileup stage:
  ```bash
  # Check file exists with expected size
  ls -lh Data/04_analysis/cnv/numbat/results/{TISSUE}/alleles.csv
  # Should be 10-100M, not <1K
  
  # Check line count (variants per cell):
  wc -l Data/04_analysis/cnv/numbat/results/{TISSUE}/alleles.csv
  # Should be >1M lines for good coverage
  ```

### For Alleloscope Specifically:
- Prepare inputs from **tumor tissue ONLY** (not normal reference)
- Aggregate variant data across **spatial barcodes** before input
- Use **consistent barcode format** across all input files (no suffix mixing)
- Check **variant count distribution** in output:
  ```bash
  # Verify haplotype assignments exist
  ls -lh Data/04_analysis/cnv/alleloscope/results/{TISSUE}/
  # Should contain: haplotype_assignments.csv, phase_probabilities.csv
  ```

---

## Complete Example: Running NUMBAT for lowseq_489

```bash
# 1. On login node (scc1) - file verification ONLY
ls -lh Data/01_inputs/fragments/lowseq_489/
ls -lh Data/01_inputs/barcodes/tissue_barcodes/lowseq_489/
ls -lh Data/04_analysis/cnv/numbat/reference/

# 2. Allocate compute node with tmux
qrsh -l h_rt=04:00:00 -pe omp 8 -P paxlab -l mem_per_core=8G
tmux new-session -s spatial_atac_work  # or attach if exists: tmux a -t spatial_atac_work

# 3. Verify compute allocation (should NOT be scc1*)

hostname

# 4. Generate reference (if not already done)
cd /projectnb/paxlab/presh/projects/spatial_atac
bash analysis/src/cnv_calling/numbat/generate_all_tissue_references.test.sh

# 5. Submit pileup job (from login node after compute is done)
qsub analysis/src/cnv_calling/numbat/run_numbat_pileup_phase_lowseq_489.qsub.sh

# 6. Monitor progress
qstat -j {JOBID}
99 x 24
I'm currently in "Ask" mode (read-only), so I can't directly write files. However, I can provide you with the complete Rmarkdown notebook content that you can use to update the file.

To apply this notebook, you have two options:

Option 1: Use VS Code to paste the content directly
Open the file: analyse_alleloscope_lowseq_488B.Rmd
Select all (Ctrl+A)
Paste the complete content I provided above
Option 2: Create the file via terminal
Save this command to create/update the file:

2. Load Alleloscope Results & Input Data
3. Clean & Filter Heatmap Data (Variance + Blacklist + Coverage)
4. Re-plot Step 6 CNV Heatmap with Cleaned Data
5. Extract & Visualize Genotype Information
6. Clonal Deconvolution & Lineage Analysis
7. Cell Type & Cluster-Specific Analysis
8. Summary & Output Report
Session Info
EOF

Or knit in interactive session (on compute node via qrsh):

Use qrsh to allocate compute resources
Open RStudio or use R --vanilla
rmarkdown::render("path/to/analyse_alleloscope_lowseq_488B.Rmd")
Output will be generated at:

/projectnb/paxlab/presh/projects/spatial_atac/Data/alleloscope/lowseq_tissue_from_existing/488B/output/analysis_downstream/
Would you like me to help you run this on a compute 
tail -f analysis/qsub_logs/numbat/{JOBNAME}_{JOBID}.log

# 7. Once pileup completes, submit analysis job
qsub analysis/src/cnv_calling/numbat/run_numbat_analysis_lowseq_489.qsub.sh
```

---

## Archiving Test Files After Completion

**CRITICAL**: After testing completes successfully:
```bash
# Archive test scripts
mkdir -p Data/01_inputs/archive/test_files_archive_$(date +%Y%m%d)
mv analysis/src/cnv_calling/numbat/*.test.sh Data/01_inputs/archive/test_files_archive_$(date +%Y%m%d)/

# Archive logs
mkdir -p analysis/qsub_logs/archive/test_logs_$(date +%Y%m%d)
mv analysis/qsub_logs/*test*.log analysis/qsub_logs/archive/test_logs_$(date +%Y%m%d)/ 2>/dev/null || true

# Verify only production scripts remain in src
ls analysis/src/cnv_calling/numbat/*.qsub.sh  # Should only see .qsub.sh files
```

---

Job Monitoring: Redirect large outputs to files rather than the console. 
If job runs successfully, you can stop working and I'll ask you to check the output files and logs when I see it finishes. 
Resource Management: Immediately terminate and debug any stalled or failing jobs to conserve computational resources.

