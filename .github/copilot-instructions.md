Spatial ATAC-seq Analysis Instructions

Project Path: /projectnb/paxlab/presh/projects/spatial_atac
Data: Spatial ATAC-seq (deepseq vs. lowseq) for patients 448B and 489.
Structure: Scripts are located in analysis/src. Key subfolders include Data (organization), alleloscope, numbat, and archr.

Coding & Documentation Standards

Documentation: Maintain a master README and individual READMEs for each analysis step.
Best Practices: Use Git version control, descriptive file naming, and thorough code commenting.
HPC (qsub) Guidelines

Efficiency: Ensure $\ge$80% utilization of requested CPU cores.
Login Nodes: Use only for brief logic tests (10–20 lines); never execute full pipelines.
Job Monitoring: Redirect large outputs to files rather than the console. Regularly inspect .out and .err logs.
Resource Management: Monitor job progress frequently. Immediately terminate and debug any stalled or failing jobs to conserve computational resources.

I'm currently trying to run alleloscope on lowseq data, for all tissues separately, and together (total 3 runs), but tissue 489 is giving me grief. Parallely, I run somatic_chr to compare the somatic snv variants from monopogen to compare lowseq and deepseq since technically thet should be the same variants. 