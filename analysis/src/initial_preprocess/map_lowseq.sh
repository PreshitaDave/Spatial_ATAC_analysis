#!/bin/bash -l
#$ -P  paxlab        # Set SCC project to charge
#$ -wd /projectnb/paxlab/Projects/CodeSpaghetti/tools/ATX_epigenomics-main/fastq2frags # Specify the working directory
#$ -pe omp 32        # Request cores
#$ -l h_rt=24:00:00  # Specify hard time limit for the job
#$ -N map       # Name job
#$ -o /projectnb/paxlab/Projects/CodeSpaghetti/results/00_preprocess/map
#$ -e /projectnb/paxlab/Projects/CodeSpaghetti/results/00_preprocess/map_error
#$ -m bea            # Send an email when the job finishes or aborts
#$ -M yetingli@bu.edu
#$ -j y              # Join error and output streams in one file


# Keep track of information related to the current job
echo "# -------------------------------------------------"
echo "Start date: $(date)"
echo "Job name: $JOB_NAME"
echo "Job ID: $JOB_ID  $SGE_TASK_ID"
echo "Running in directory: $PWD"
echo "# -------------------------------------------------"

module load miniconda
# conda activate ATX
chromap -t 32 \
        --preset atac \
        -x Refseq/index \
        -r Refseq/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
        -1 D01596_NG03920_linker2_R1.fastq.gz \
        -2 D01596_NG03920_linker2_R2.fastq.gz \
        -o D01596_NG03920_bc96_aln.bed \
        -b D01596_NG03920_linker2_R2.fastq.gz \
        --barcode-whitelist bc96.txt.gz \
        --read-format bc:22:29,bc:60:67,r1:0:-1,r2:117:-1
awk 'BEGIN{FS=OFS=" "}{$4=$4"-1"}4' D01596_NG03920_bc96_aln.bed > D01596_NG03920_bc96_temp.bed
sed 's/ /\t/g' D01596_NG03920_bc96_temp.bed > D01596_NG03920_bc96_fragments.tsv
bgzip -c D01596_NG03920_bc96_fragments.tsv > D01596_NG03920_bc96_fragments.tsv.gz

echo "Job finished: $(date +%F)"