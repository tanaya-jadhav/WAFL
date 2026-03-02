
# WAFL: Workflow for Annotating & Filtering Long Read Data

WAFL is a comprehensive bioinformatics workflow designed for the annotation and filtering of genetic variants from long-read sequencing data. It integrates Single Nucleotide Variants (SNVs), Indels, Structural Variants (SVs), and repeat expansions to provide a detailed analysis of a sample's genetic makeup. The workflow aims to identify and prioritize potentially pathogenic variants by leveraging various annotation databases and computational tools.

## Workflow Overview

The WAFL workflow performs the following key steps:

1.  **Structural Variant (SV) Processing:**
    *   Annotates SVs using SVAFotate with population frequencies (gnomAD, HPRC, CoLoRSdb) and platform controls.
    *   Further annotates SVs using the Ensembl Variant Effect Predictor (VEP) for predicted effects.
    *   Filters SVs based on allele frequency cutoffs from various databases.

2.  **SNV/Indel Processing:**
    *   Filters SNVs/Indels to relevant genomic regions (e.g., chromosomes 1-22, X, Y) and normalizes them.
    *   Annotates SNVs/Indels using ECHTVAR with population frequencies and filters them based on thresholds.
    *   Filters for genic variants using a RefSeq BED file.
    *   Annotates genic variants with VEP, SpliceAI, REVEL, LoFtool, CADD, and Clinvar information.
    *   Filters variants based on Clinvar significance and predicted impact.

3.  **Variant Inheritance and Ranking:**
    *   Identifies compound heterozygous variants by combining SNV/Indel and SV information, considering inheritance patterns (trio, duo, singleton).
    *   Ranks variants based on multiple evidence sources, including population frequency (gnomAD, HPRC, rCGL), pathogenicity predictions (LoFtool, Revel, CADD, SpliceAI), Clinvar/HGMD annotations, gene constraint (LOEUF, regional constraint), de novo status, and inheritance patterns.
    *   Separates variants into dominant and recessive categories based on their predicted inheritance models and pathogenicity features.

4.  **Methylation Analysis:**
    *   Creates methylation profiles for samples.
    *   Identifies variants located within Differentially Methylated Regions (DMRs) and analyzes their overlap with imprinted genes.

5.  **Repeat Expansion Analysis:**
    *   Identifies and reports on potential repeat expansions using TRGT and a specialized repeat database.

6.  **Final Reporting:**
    *   Compiles all processed variant information (dominant, recessive, SV+SNV comphets, filtered SVs, variants in DMRs, repeat expansions) into a single, comprehensive Excel report.

## Key Files

*   `Snakefile`: The main Snakemake workflow definition file, orchestrating all the analysis steps.
*   `config.yaml`: Central configuration file where all input file paths, sample information, filtering thresholds, and tool parameters are defined. **This file must be customized for each analysis.**
*   `config_template.yaml`: A template for creating the `config.yaml` file.
*   `environment.yml`: Defines the Conda environment required to run the workflow, listing all necessary software and their versions.
*   `snakemake.script.sh`: A shell script to execute the Snakemake workflow with specific configurations (e.g., using Singularity containers, setting job parallelism).
*   `submissionscript.sh`: An example Slurm submission script to run the `snakemake.script.sh` on a cluster.
*   `workflow_scripts/`: A directory containing various Python scripts used for specific annotation, filtering, and analysis tasks within the workflow.

## Required Inputs

To run the WAFL workflow, the following inputs must be specified in the `config.yaml` file:

*   **Sample Information:**
    *   `Sample`: The identifier for the sample being analyzed (e.g., "236734WN").
*   **Input Variant Files:**
    *   `SNV_vcf`: Path to the input SNV/Indel VCF file.
    *   `SV_vcf`: Path to the input Structural Variant (SV) VCF file.
    *   `trgt_vcf`: Path to the TRGT VCF file for repeat expansion analysis.
*   **Pedigree File:**
    *   `ped`: Path to the pedigree file, defining family relationships (proband, parents).
*   **Methylation Analysis Files:**
    *   `cpgfile_prefix`: Prefix for CpG files.
    *   `cpg_island_regions`: Path to the CpG island regions BED file.
    *   `controlmethylationprofile`: Path to the control methylation profile TSV.
    *   `imprintedgenesbed`: Path to the imprinted genes BED file.
*   **Annotation and Reference Files:**
    *   `hg38_ref`: Path to the human reference genome (hg38) FASTA file.
    *   Various annotation archive files (e.g., `gnomad_archive`, `clinvar_archive`, `hprc_archive`, `spliceAI_snvarchive`, `revel_archive`, `lof_archive`, `cadd_snvarchive`, `clinvar_processed`, `repeatsdb`). Refer to `config.yaml` for a complete list.
    *   Reference gene annotation files (e.g., `refseq_genes_bed`, `refseq_cache_dir`).
*   **Filtering Thresholds:**
    *   `gnomad_cutoff`, `hprc_cutoff`, `colorsdb_cutoff`: Allele frequency cutoffs for SV filtering.
    *   `gnomad_popmax_af_threshold`, `hprc_FRQ_threshold`: Population frequency thresholds for ECHTVAR annotation.
*   **Tool Configurations:**
    *   Container paths for tools like VEP, ECHTVAR, SLIVAR.

## Setup and Execution

1.  **Environment Setup:**
    Create and activate the Conda environment using the provided `environment.yml` file:
    ```bash
    conda env create -f environment.yml
    conda activate wafl_env
    ```
    *(Note: Ensure you have Conda installed and configured.)*

2.  **Configuration:**
    *   Copy `config_template.yaml` to `config.yaml`.
    *   Edit `config.yaml` to specify the correct paths for all input files, the sample ID, and any necessary parameters according to your dataset.

3.  **Workflow Execution:**
    The workflow can be executed using Snakemake. It is recommended to use the provided shell scripts for execution, which handle Singularity container usage and job submission.

    *   **Local Execution (using `snakemake.script.sh`):**
        ```bash
        bash snakemake.script.sh
        ```
        *(This command will run Snakemake with default settings as defined in the script and config.yaml)*

    *   **Cluster Execution (e.g., Slurm using `submissionscript.sh`):**
        ```bash
        bash submissionscript.sh
        ```
        *(This script submits the Snakemake job to a Slurm cluster. Adjustments to the `submissionscript.sh` might be needed based on your cluster's configuration.)*

## Output

The primary output of the WAFL workflow is an Excel file named `[SampleID].CompiledReport.xlsx` located in the `output/[SampleID]/` directory. This file consolidates the results from various analysis modules, providing a summarized view of the prioritized variants. Intermediate VCF and TSV files are generated throughout the workflow.
