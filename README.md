# Spatial_ATAC_analysis

This repository contains a pipeline for the processing and analysis of single-cell ATAC-seq data, potentially integrated with other modalities like RNA-seq and Xenium In Situ data. It covers steps from raw data preparation to downstream analyses such as CNV calling and variant calling.

## File Structure and Description

The following table outlines the main files in this repository and a brief description of their purpose.

| File                              | Description                                                                 |
| :-------------------------------- | :-------------------------------------------------------------------------- |
| `00_create_index.sh`              | Script to generate necessary genomic indices for alignment.                 |
| `01_trim_fastq.sh`                | Script for quality trimming and adapter removal from raw FASTQ reads.       |
| `02_map_fastq.sh`                 | Script to align processed FASTQ reads to a reference genome.                |
| `03_download_ref.sh`              | Script to download reference genome and annotation files.                   |
| `1_Lib_install.Rmd`               | R Markdown script to install all required R packages for the analysis.      |
| `2_scATAC_giotto_obj_creation.Rmd`| R Markdown script for creating a Giotto object from scATAC-seq data.      |
| `3_ArchR_ATAC_analysis.Rmd`       | R Markdown script for performing core scATAC-seq analysis using ArchR.      |
| `4_CNV_calling_ATAC.Rmd`          | R Markdown script to call Copy Number Variants from ATAC-seq data.          |
| `5_Xenium_giotto_analysis.Rmd`    | R Markdown script for analyzing Xenium In Situ data, likely using Giotto.   |
| `5_Xenium_giotto_analysis.nb.html`| HTML output of the Xenium Giotto analysis R Markdown notebook.              |
| `6_compare_atac_rna.Rmd`          | R Markdown script for integrating and comparing scATAC-seq and scRNA-seq data. |
| `6_compare_atac_rna.nb.html`      | HTML output of the ATAC-RNA comparison R Markdown notebook.                 |
| `7_variant_calling_atac.Rmd`      | R Markdown script for calling genetic variants from ATAC-seq data.          |
| `7_variant_calling_atac.nb.html`  | HTML output of the ATAC variant calling R Markdown notebook.                |
| `cram_bam.sh`                     | Script for converting between CRAM and BAM file formats.                    |
| `map_lowseq.sh`                   | Script for mapping low-coverage sequencing data.                            |

## Getting Started

To run this pipeline, please follow these general steps:

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/yourusername/your-repo-name.git
    cd your-repo-name
    ```
2.  **Install dependencies:** Ensure all necessary software (e.g., samtools, BWA) and R packages are installed. You can use `1_Lib_install.Rmd` to install R packages.
3.  **Prepare reference data:** Run `03_download_ref.sh` to get the required reference files.
4.  **Process raw sequencing data:** Start with `00_create_index.sh`, `01_trim_fastq.sh`, and `02_map_fastq.sh` to prepare your raw sequencing data.

## Usage

-   Shell scripts (`.sh`) can be executed directly from your terminal.
-   R Markdown files (`.Rmd`) can be opened and run interactively in RStudio, or knitted to HTML/PDF documents using `rmarkdown::render()`.

## Dependencies

This project relies on a combination of bioinformatics tools and R packages. Key dependencies include:

-   R (version 4.0 or higher recommended)
-   Samtools
-   BWA
-   Specific R packages (details in `1_Lib_install.Rmd`, including ArchR, Giotto, etc.)

## Contact

For any questions or issues, please contact [Your Name/Email Address/GitHub Profile].
