Spatial ATAC-seq Analysis Instructions

Project Path: /projectnb/paxlab/presh/projects/spatial_atac
Data: Spatial ATAC-seq (deepseq vs. lowseq) for patients 448B and 489.
Structure: Scripts are located in analysis/src. Key subfolders include Data (organization)

Coding & Documentation Standards

Documentation: Maintain a master README and individual READMEs for each analysis step.
Best Practices: Use Git version control, descriptive file naming, and thorough code commenting.
HPC (qsub) Guidelines
Always check if the paths point to the correct files and directories before submitting jobs. Use absolute paths. And check if the output dir exists and is correct. 
Organize new scripts properly in the correct folder and also make a new one if necessary to avoid clutter. Organize the data and output also in folders that make sense. Move anything old and unused to an archive folder.
Symlink Management Policy:
- ALWAYS label symlinks with ".lnk" suffix so they're distinguishable from actual files.
- **Internal symlinks** (pointing within /projectnb/paxlab/presh/): Remove and replace with absolute paths in scripts.
- **External symlinks** (pointing outside /projectnb/paxlab/presh/): Keep in designated location `Data/01_inputs/archive/symlink_external/`.
- **Before organizing files**: Create `symlink_mapping_BEFORE.txt` documenting all symlinks.
- **After organizing files**: Create `symlink_mapping_AFTER.txt` and verify all paths are correct.
- **If scripts break after path changes**: Reference mapping files to diagnose missing paths and update affected scripts.
- **Update ORGANIZATION_SUMMARY.md** with new file locations and symlink policies for future reference.

Efficiency: Ensure $\ge$80% utilization of requested CPU cores.

File Organization Standards

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

Interactive Compute Session Policy:
- For any code testing, diagnostics, or script execution, do not stay on the SCC login node.
- **ALWAYS verify hostname before running code**: Run `hostname` - if it starts with "scc1", you're on LOGIN NODE (danger zone).
- **ALWAYS open a visible terminal** for the user to see real-time command execution and output.
- **ALWAYS execute `qrsh` inside a `tmux` session** to ensure persistence and visibility.

WORKFLOW (DO THIS AT THE START OF EVERY SESSION):
1. Open a terminal with async mode so user sees live execution
2. Check hostname: if login node (scc1*), proceed to step 3
3. Start persistent session: tmux new -As spatial_atac_work (or screen -S spatial_atac_work)
4. Request compute node: qrsh -l h_rt=16:00:00 -pe omp 1 -P paxlab -l mem_per_core=8G
5. Verify allocation: hostname should NOT be scc1* (if still scc1, qrsh failed - use qsub instead)
6. Load required modules: module load R (or other dependencies)
7. Execute code on compute node
8. Show user EVERY step in terminal with clear [STEP N] labels

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

Job Monitoring: Redirect large outputs to files rather than the console. 
If job runs successfully, you can stop working and I'll ask you to check the output files and logs when I see it finishes. 
Resource Management: Immediately terminate and debug any stalled or failing jobs to conserve computational resources.





Current focus:
I'm currently trying to run alleloscope on lowseq data, for all tissues separately, and together (total 3 runs), but tissue 489 is giving me grief. Parallely, I run somatic_chr to compare the somatic snv variants from monopogen to compare lowseq and deepseq since technically thet should be the same variants.