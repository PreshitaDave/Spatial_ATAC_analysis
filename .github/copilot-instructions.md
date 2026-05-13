Spatial ATAC-seq Analysis Instructions

Project Path: /projectnb/paxlab/presh/projects/spatial_atac
Data: Spatial ATAC-seq (deepseq vs. lowseq) for patients 448B and 489.
Structure: Scripts are located in analysis/src. Key subfolders include Data (organization)

Coding & Documentation Standards

Documentation: Maintain a master README and individual READMEs for each analysis step.
Best Practices: Use Git version control, descriptive file naming, and thorough code commenting.
HPC (qsub) Guidelines
Always check if the paths point to the correct files and directories before submitting jobs. Use absolute paths. And check if the output dir exists and is correct. 
Organize new scripts properly in the correct folder and also make a new one if necessary to avoid clutter. Organize the data and output also in folders that make sense. 

Efficiency: Ensure $\ge$80% utilization of requested CPU cores.
Login Nodes: Use only for brief logic tests (10–20 lines); never execute full pipelines. I want you to run the full pipelines on the compute nodes, not on the login nodes, or submit a qsub job. If terminal fails on login node, immediately move to qsub job submission. Try to estimate the time it'll take to run it on login node and if it's more than 2 minutes, submit a qsub job instead.
Job Monitoring: Redirect large outputs to files rather than the console. Regularly inspect .out and .err logs.
Resource Management: Immediately terminate and debug any stalled or failing jobs to conserve computational resources.

Current focus:
I'm currently trying to run alleloscope on lowseq data, for all tissues separately, and together (total 3 runs), but tissue 489 is giving me grief. Parallely, I run somatic_chr to compare the somatic snv variants from monopogen to compare lowseq and deepseq since technically thet should be the same variants. 