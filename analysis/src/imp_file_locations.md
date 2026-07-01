Important file locations : 

FOR ORIGINAL DATA RECEIVED - 
1. The original cram -> bam file for deepseq (both tissues) 
/projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D1942/D1942_deepseq/processed_data/D1942_deepseq.bam
This has ba and bb barcode tags of 8bp each 

2. Lowseq bam file (both tissues) I had to align from the fastq files and the final bam file is present at -> 
/projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D1942/D01942_NG05827_lowseq/lowseq.sorted.bam 
In the lowseq bam file you have  CB:Z:CTGAGCCACCTATGTC (all 16bp together)


All these bam files have be recreated separately for each tissue object downstream!

The bam files located here : /projectnb/paxlab/presh/projects/spatial_atac/Data/01_inputs/bam
Have deepseq bams with CB:Z tag with 16bp! 
The fragments files have the barcode ending with '-1' suffix 

