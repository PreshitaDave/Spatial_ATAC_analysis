#!/bin/bash
#$ -P paxlab        # Set SCC project to charge
#$ -wd /projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D01942_NG05827/lowseq_alignment # Specify the working directory
#$ -pe omp 8       # Request cores
#$ -l h_rt=12:00:00  # Specify hard time limit for the job
#$ -N map       # Name job
#$ -o /projectnb/paxlab/presh/src/Spatial_ATAC_analysis/qsub_logs
#$ -e /projectnb/paxlab/presh/src/Spatial_ATAC_analysis/qsub_logs
#$ -M preshita@bu.edu


# Keep track of information related to the current job
echo "# -------------------------------------------------"
echo "Start date: $(date)"
echo "Job name: $JOB_NAME"
echo "Job ID: $JOB_ID  $SGE_TASK_ID"
echo "Running in directory: $PWD"
echo "# -------------------------------------------------"


CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"

# module load bbmap/38.16
module load java
# module unload python3
# module load miniconda
conda activate '/projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D01942_NG05827/chromap_env'

cd /projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D01942_NG05827/lowseq_alignment
# Adjust these paths as needed
REF_DIR="/projectnb/paxlab/Projects/CodeSpaghetti/tools/ATX_epigenomics-main/fastq2frags/Refseq"
OUTPUT_DIR="/projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D01942_NG05827/lowseq_alignment/"
BARCODE_FILE="/projectnb/paxlab/Projects/CodeSpaghetti/tools/ATX_epigenomics-main/fastq2frags/bc96.txt.gz"

# The files are here 
SAMPLE_PREFIX="D01942_NG05827_200403"

  L1_R1="${OUTPUT_DIR}${SAMPLE_PREFIX}_linker1_R1.fastq.gz"
  L1_R2="${OUTPUT_DIR}${SAMPLE_PREFIX}_linker1_R2.fastq.gz"
  L2_R1="${OUTPUT_DIR}${SAMPLE_PREFIX}_linker2_R1.fastq.gz"
  L2_R2="${OUTPUT_DIR}${SAMPLE_PREFIX}_linker2_R2.fastq.gz"
  # ALN_BED="${SAMPLE_PREFIX}_aln.bed"
  ALN_SAM="${OUTPUT_DIR}${SAMPLE_PREFIX}_aln.sam"

echo "--- Verifying input files ---"
ls -l "$L2_R1" "$L2_R2" "$BARCODE_FILE" "$REF_DIR/index" "$REF_DIR/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
echo "--- End file verification ---"

# Step 3: Chromap alignment
echo "--- Starting Chromap Alignment ---"
echo "Chromap command:"
echo "chromap -t $NSLOTS \\"
echo "  --preset atac \\"
echo "  -x \"$REF_DIR/index\" \\"
echo "  -r \"$REF_DIR/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa\" \\"
echo "  -1 \"$L2_R1\" \\"
echo "  -2 \"$L2_R2\" \\"
echo "  --SAM -o \"$ALN_SAM\" \\"
echo "  -b \"$L2_R2\" \\"
echo "  --barcode-whitelist \"$BARCODE_FILE\" \\"
echo "  --read-format bc:22:29,bc:60:67,r1:0:-1,r2:117:-1"
echo "------------------------------"

# Use NSLOTS to automatically pick up the number of requested cores from SGE
# (Your -pe omp 4 sets NSLOTS=4)
chromap -t $NSLOTS \
  --preset atac \
  -x "$REF_DIR/index" \
  -r "$REF_DIR/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa" \
  -1 "$L2_R1" \
  -2 "$L2_R2" \
  --SAM -o "$ALN_SAM" \
  -b "$L2_R2" \
  --barcode-whitelist "$BARCODE_FILE" \
  --read-format bc:22:29,bc:60:67,r1:0:-1,r2:117:-1

# Check Chromap exit status
if [ $? -eq 0 ]; then
  echo "Chromap alignment completed successfully."
else
  echo "Chromap alignment FAILED. Check error logs."
  exit 1 # Exit with error status
fi


  # # Step 4: Convert BED to fragments
  # awk 'BEGIN{FS=OFS=" "}{$4=$4"-1"}4' "$ALN_BED" > "$TEMP_BED"
  # sed 's/ /\t/g' "$TEMP_BED" > "$FRAG_TSV"
  # bgzip -c "$FRAG_TSV" > "$FRAG_TSVGZ"
  # rm "$TEMP_BED"
  
echo "Job finished: $(date +%F)"
  
