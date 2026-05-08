#!/bin/bash -l

#$ -P paxlab            # Specify the SCC project name you want to use
#$ -l h_rt=48:00:00     # Specify the hard time limit for the job
#$ -N somatic_deep        # Give job a name
#$ -l mem_per_core=8G  # Request memory per core (adjust as needed)
#$ -pe omp 16         # Request 12 cores for parallel processing
#$ -o /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs   # Specify the path for standard output logs
#$ -e /projectnb/paxlab/presh/projects/spatial_atac/analysis/qsub_logs  # Specify the path for error logs

module load python3
module load R

MONOPOGEN_PATH="/projectnb/paxlab/presh/software/Monopogen"
REGION_LST="${MONOPOGEN_PATH}/resource/GRCh38.region.lst"
REF_FA="/projectnb/paxlab/presh/projects/spatial_atac/Data/hg38_resources/Homo_sapiens_assembly38.fasta"


python  "${MONOPOGEN_PATH}/src/Monopogen.py"  somatic  \
    -a   "${MONOPOGEN_PATH}/apps" -r /projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/chr_region.lst   -t 8 \
    -i  "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq"  -l  /projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq_cell_data.csv   -s featureInfo     \
    -g   "${REF_FA}"



python  "${MONOPOGEN_PATH}/src/Monopogen.py"  somatic  \
    -a  "${MONOPOGEN_PATH}/apps"  -r /projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/chr_region.lst   -t 8 \
    -i  "/projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq"  -l  /projectnb/paxlab/presh/projects/spatial_atac/Data/variant_calling/deepseq_cell_data.csv   -s cellScan     \
    -g   "${REF_FA}"



