#!/bin/bash
#$ -P paxlab        # Set SCC project to charge
#$ -pe omp 4       # Request cores
#$ -l h_rt=12:00:00  # Specify hard time limit for the job
#$ -N download       # Name job
#$ -o /projectnb/paxlab/presh/src/Spatial_ATAC_analysis/qsub_logs
#$ -e /projectnb/paxlab/presh/src/Spatial_ATAC_analysis/qsub_logs
#$ -M preshita@bu.edu


# --- Configuration ---
MANIFEST_FILE="/projectnb/paxlab/presh/Data/variant_calling/phased_hg38_1000G_ref/phased-manifest_July2021.tsv"
# IMPORTANT: Replace this with the actual base URL where your files are hosted.
# Example: BASE_URL="https://example.com/data/"
BASE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/" # <--- EDIT THIS LINE

# --- Script Start ---
echo "Starting download from manifest: $MANIFEST_FILE"
echo "Base URL: $BASE_URL"

# Check if wget or curl is available
if ! command -v wget &> /dev/null && ! command -v curl &> /dev/null; then
    echo "Error: Neither 'wget' nor 'curl' found. Please install one of them."
    exit 1
fi

# Loop through each line in the manifest file, extracting the filename
# We use 'awk' to get the first column (the filename)
while read -r line; do
    FILENAME=$(echo "$line" | awk '{print $1}')

    if [[ -z "$FILENAME" ]]; then # Skip empty lines
        continue
    fi

    FULL_URL="${BASE_URL}${FILENAME}"

    echo -e "\nDownloading: $FILENAME"
    echo "From: $FULL_URL"

    # --- Download the file ---
    # Use wget if available, otherwise fall back to curl
    if command -v wget &> /dev/null; then
        wget -c "$FULL_URL" # -c for continuing partial downloads
    elif command -v curl &> /dev_null; then
        curl -L -O "$FULL_URL" # -L to follow redirects, -O to use original filename
    fi

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Download failed for $FILENAME."
    else
        echo "Download successful for $FILENAME."
    fi

done < "$MANIFEST_FILE"

echo -e "\nDownload process complete."
