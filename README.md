# Optimizing Gene Integration: Hot and Cold Loci Prediction in *E. coli*

## Short Description

This project provides a robust Snakemake-based bioinformatics pipeline designed to process RNA-seq data from *Escherichia coli* BL21. Its primary goal is to generate a comprehensive, pre-processed dataset essential for training machine learning models that predict "hot" (high expression) and "cold" (low expression) integration loci for foreign genes. 

## Problem Statement

Efficient expression of integrated foreign genes is a critical bottleneck in biotechnology and synthetic biology. Suboptimal integration sites can lead to low expression levels, impacting the viability and efficiency of genetically modified organisms. Current methods often involve trial-and-error, which is resource-intensive, time-consuming, and costly. There is a pressing need for a predictive approach to identify genomic locations that are likely to support desired expression levels upon gene integration.

## Solution & Approach

We propose a novel method to predict integration points in *Escherichia coli* (specifically the BL21 strain) that are associated with varying levels of gene expression (elevated or reduced). This prediction capability enables optimized genetic modifications, leading to significant savings in time and financial resources for biotechnological applications.

Our solution is built around a comprehensive Snakemake workflow that automates the entire RNA-seq data processing pipeline:

1.  **Data Acquisition:** Raw RNA-seq reads are sourced from NCBI SRA.
2.  **Quality Control:** Rigorous quality assessment and trimming ensure clean and reliable data.
3.  **Expression Quantification:** Reads are mapped and quantified against a reference transcriptome to determine gene expression levels.
4.  **Data Normalization & Aggregation:** Expression data is normalized (using RPKM) and aggregated into a unified matrix.
5.  **Genomic Feature Engineering:** The processed data is further enriched with crucial genomic features, such as gene annotations and distances to known *E. coli* replication origins (ori-sites), which are important indicators for gene expression.

The output of this pipeline is a meticulously prepared dataset, ready for machine learning model development to pinpoint optimal (hot) and suboptimal (cold) integration loci.

## Key Features



## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

Before running the pipeline, ensure you have the following installed:

*   **Miniconda3/Anaconda3**: Snakemake heavily relies on `conda` for environment management.
*   **Snakemake**: Install it within your base `conda` environment or via `pip`.
    ```bash
    conda install -c bioconda -c conda-forge snakemake
    ```
*   **Python 3.10**: The project's Python scripts are developed with this version.

### Installation

1.  **Clone the repository:**
    ```bash
    git clone [YOUR_REPOSITORY_URL]
    cd [YOUR_REPOSITORY_NAME]
    ```
2.  **Create Conda environments:**
    The pipeline uses specific `conda` environments for its tools. These are defined in the `envs/` directory.
    ```bash
    conda env create -f envs/fastp.yaml
    conda env create -f envs/kallisto.yaml
    conda env create -f envs/sra_toolkit.yaml # If you uncomment SRA download rules
    # ... and any other environment files you have in 'envs/'
    ```
    *Note: Ensure your `envs` directory contains the necessary `.yaml` files.*
	*Note: .yaml envs must be placed in /pipeline/envs/*

### Data Preparation

This pipeline expects certain data to be pre-existing in specific locations.

1.  **SRA Accession List:**
    Create a file named `rnaseq_paired.tsv` (or as specified in your `cfg.yaml` under `samples_file`) in your `data/` directory. This file should contain a simple list of NCBI SRA Accession IDs, one per line.
    Example `data/rnaseq_paired.tsv`:
    ```
    SRR1234567
    SRR1234568
    SRR1234569
    ```

2.  **Reference Transcriptome and Kallisto Index:**
    You need to provide a reference transcriptome in FASTA format and a pre-built Kallisto index for it.
    *   Place your transcriptome FASTA file (e.g., `ecoli_bl21_transcriptome.fna`) in the `pipeline/data/` directory.
    *   Create a directory for the Kallisto index (e.g., `data/kallisto_index/`) and place your index file (e.g., `kallisto.idx`) inside it.
    *   *How to obtain:*
        *   **Genome/Annotation:** Download genome (`.fasta`) and annotation (`.gff`/`.gtf`) from NCBI.
        *   **Transcriptome:** Use `gffread` (from StringTie/GFFRead) to extract the transcriptome FASTA from the genome and GFF file (e.g., `gffread -w ecoli_bl21_transcriptome.fna -g ecoli_bl21_genome.fasta ecoli_bl21_annotation.gff`).
        *   **Kallisto Index:** Build the index using `kallisto index -i pipeline/data/kallisto_index/kallisto.idx pipeline/data/ecoli_bl21_transcriptome.fna`.
    

3.  **Configuration File (`cfg.yaml`):**
    Ensure your `cfg.yaml` file (located in the project root) is correctly configured with paths to your `samples_file`, `transcriptome_fasta`, `kallisto_index_dir`, and any desired `fastp` or `kallisto` options. An example configuration is provided below:
    ```yaml
    # cfg.yaml
    samples_file: "data/rnaseq_paired.tsv" 
    transcriptome_fasta: "data/ecoli_bl21_transcriptome.fna"
    kallisto_index_dir: "data/kallisto_index/" 

    fastp_options: "--detect_adapter_for_pe -q 20 -l 30"
    kallisto_options: "--bias" 
    kallisto_bootstrap: 100
    ```
	*Note: config file must be placed in  /pipeline/scripts/*

## Pipeline

SRA-toolkit fetches .fasta files and zips them from .tsv in `/pipeline/data/`.
Files go into `/pipeline/raw/` and then getting trimmed by fastp. 
Trimmed reads are put in `/pipeline/data/`.
Trimmed files are getting quantified by kallisto. Subdirectories named by reads are created in ?, each of those contains `abundance.tsv` file with all required data. 
Then kallisto merges them in one matrix, which is placed at `/pipeline/results/kallisto_tpm.matrix.tsv`
## Usage

Once all prerequisites are met and data is prepared, you can run the Snakemake pipeline from the project's root directory.

```bash
snakemake --rerun-incomplete --printshellcmds --conda-frontend conda --keep-going > snakemake.log 2>&1 &
```

- --rerun-incomplete: Ensures incomplete jobs are rerun.
    
- --printshellcmds: Prints the shell commands executed by Snakemake.
    
- --conda-frontend conda: Ensures Snakemake uses conda to manage environments.
    
- --keep-going: Continues with other jobs if one job fails.
    
- > snakemake.log 2>&1 &: Redirects all output to snakemake.log and runs the process in the background.
	   
-  *additional:*  --cores N: Specifies the number of CPU cores to use. Adjust as needed.

## Output

The primary output of this pipeline is the aggregated gene expression matrix, kallisto_tpm_matrix.tsv, located in the pipeline/results/ directory.

- `pipeline/results/kallisto_tpm_matrix.tsv`: A tab-separated values (TSV) file containing the RPKM-normalized expression levels for all genes across all processed samples. This matrix serves as the main input for downstream machine learning model training.
    
- `pipeline/qc/fastp/`: Directory containing HTML and JSON reports from fastp for each sample, useful for quality assessment.
    
- `pipeline/kallisto_results/`: Subdirectories for each sample, containing abundance.tsv files generated by Kallisto.
    
- `pipeline/logs/`: Directory for Snakemake and tool-specific logs.

## Postprocessing
### Normalisation
Aggregated matrix must be normalised. Normalisation in this project is done via RPKM python script *normalisation.py*

### Annotating
Output matrix contains "tech-y" gene names. In order to get proper, human-understandable gene names, we offer you a two-step solution.

Step 1. Create a comparison map using *names.py* using original .gfa genome annotation. Please note that it works for prokaryotes only currently and consider to tweak this script for your needs.

Step 2. With comparison map, and original kallisto matrix, 
