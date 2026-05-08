#!/bin/bash
#$ -P paxlab        # Set SCC project to charge
#$ -wd /projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D01942_NG05827/lowseq_alignment # Specify the working directory
#$ -pe omp 16        # Request cores
#$ -l h_rt=12:00:00  # Specify hard time limit for the job
#$ -N trim       # Name job
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

module load bbmap/38.16
module load java
module unload python3
module load miniconda
conda activate /projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D01942_NG05827/chromap_env
  

cd /projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D01942_NG05827/lowseq_alignment
# Adjust these paths as needed
REF_DIR="/projectnb/paxlab/Projects/CodeSpaghetti/tools/ATX_epigenomics-main/fastq2frags/Refseq/Refdata_scATAC_MAESTRO_GRCh38_1.1.0"
FASTQ_DIR="/projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D01942_NG05827/preprocessing/fastq"
BARCODE_FILE="/projectnb/paxlab/Projects/CodeSpaghetti/tools/ATX_epigenomics-main/fastq2frags/bc96.txt.gz"
CORES=16

# SAMPLE=$(basename "$R1" | sed 's/_R1.*//')

LINKER1="GTGGCCGATGTTTCGCATCGGCGTACGACT"
LINKER2="ATCCACGTGCTTGAGAGGCCAGAGCATTCG"

for SAMPLE_PREFIX in $(ls "$FASTQ_DIR/"*R1*fastq.gz | xargs -n 1 basename | sed 's/_R1.*//' | sort -u); do
  echo "Processing sample $SAMPLE_PREFIX..."
  
  # Find R1/R2 files using the sample prefix
  R1=$(ls "$FASTQ_DIR/${SAMPLE_PREFIX}"_R1*fastq.gz | head -1)
  R2=$(ls "$FASTQ_DIR/${SAMPLE_PREFIX}"_R2*fastq.gz | head -1)

  # Check if R1/R2 files were found
  if [[ -z "$R1" || -z "$R2" ]]; then
    echo "Warning: Could not find paired R1/R2 files for sample $SAMPLE_PREFIX. Skipping."
    continue # Skip to next sample
  fi

  echo "  R1 file: $R1"
  echo "  R2 file: $R2"
  
  
  # Output filenames
  L1_R1="${SAMPLE_PREFIX}_linker1_R1.fastq.gz"
  L1_R2="${SAMPLE_PREFIX}_linker1_R2.fastq.gz"
  L2_R1="${SAMPLE_PREFIX}_linker2_R1.fastq.gz"
  L2_R2="${SAMPLE_PREFIX}_linker2_R2.fastq.gz"
  # ALN_BED="${SAMPLE}_aln.bed"
  STATS_L1="${SAMPLE_PREFIX}_stats.linker1.txt"
  STATS_L2="${SAMPLE_PREFIX}_stats.linker2.txt"
  # TEMP_BED="${SAMPLE}_temp.bed"
  # FRAG_TSV="${SAMPLE}_fragments.tsv"
  # FRAG_TSVGZ="${SAMPLE}_fragments.tsv.gz"

  # Step 1: Filter linker 1
  /projectnb/paxlab/Projects/CodeSpaghetti/tools/ATX_epigenomics-main/fastq2frags/bbmap/bbduk.sh \
    in1="$R1" \
    in2="$R2" \
    outm1="$L1_R1" \
    outm2="$L1_R2" \
    k=30 \
    mm=f \
    rcomp=f \
    restrictleft=103 \
    skipr1=t \
    hdist=3 \
    stats="$STATS_L1" \
    threads="$CORES" \
    literal="$LINKER1"

  # Step 2: Filter linker 2
  /projectnb/paxlab/Projects/CodeSpaghetti/tools/ATX_epigenomics-main/fastq2frags/bbmap/bbduk.sh \
    in1="$L1_R1" \
    in2="$L1_R2" \
    outm1="$L2_R1" \
    outm2="$L2_R2" \
    k=30 \
    mm=f \
    rcomp=f \
    restrictleft=65 \
    skipr1=t \
    hdist=3 \
    stats="$STATS_L2" \
    threads="$CORES" \
    literal="$LINKER2"

  # # Step 3: Chromap alignment
  # chromap/chromap \
  #   -t "$CORES" \
  #   --preset atac \
  #   -x "$REF_DIR/GRCh38_genome.index" \
  #   -r "$REF_DIR/GRCh38_genome.fa" \
  #   -1 "$L2_R1" \
  #   -2 "$L2_R2" \
  #   -o "$ALN_BED" \
  #   -b "$L2_R2" \
  #   --barcode-whitelist "$BARCODE_FILE" \
  #   --read-format bc:22:29,bc:60:67,r1:0:-1,r2:117:-1
  # 
  # # Step 4: Convert BED to fragments
  # awk 'BEGIN{FS=OFS=" "}{$4=$4"-1"}4' "$ALN_BED" > "$TEMP_BED"
  # sed 's/ /\t/g' "$TEMP_BED" > "$FRAG_TSV"
  # bgzip -c "$FRAG_TSV" > "$FRAG_TSVGZ"
  # rm "$TEMP_BED"
  

  # Step 3: Chromap alignment
  chromap -t 16 \
    --preset atac \
    -x "$REF_DIR/index" \
    -r "$REF_DIR/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa" \
    -1 "$L2_R1" \
    -2 "$L2_R2" \
    --SAM -o "$ALN_SAM" \
    -b "$L2_R2" \
    --barcode-whitelist "$BARCODE_FILE" \
    --read-format bc:22:29,bc:60:67,r1:0:-1,r2:117:-1
  

  echo "Done $SAMPLE_PREFIX"
  
done


