
# WAFL: Workflow for Annotating & Filtering Long Read Data

WAFL is a Snakemake workflow for annotating and filtering genetic variants from long-read sequencing data. It integrates Single Nucleotide Variants (SNVs), Indels, Structural Variants (SVs), methylation, and repeat expansions for a single proband, and identifies and prioritizes potentially pathogenic variants using population databases, pathogenicity predictors, and inheritance information from the pedigree.

## Workflow Overview

The WAFL workflow performs the following key steps:

1.  **Structural Variant (SV) Processing:**
    *   Annotates SVs with population allele frequencies using SVAFotate (gnomAD, HPRC, CoLoRSdb, and other public sources), then with internal platform controls (CDLS).
    *   Annotates SVs with AnnotSV.
    *   Tags each AnnotSV record with a refined location category and the nearest intron-exon junction using the MANE RefSeq GTF (`tag_annotsv.py`).
    *   Filters SVs on gnomAD, HPRC, CoLoRSdb, and CDLS frequency cutoffs, writing filtered SVs and inversions separately (`filterSVs.py`).
    *   Matches Sawfish BND breakends against HiTSV insertions and separates the results into breakend/insertion matches, all insertions, filtered CNVs, and non-coding CNVs, flagging segmental duplications (`match_bnd_ins_TJ.py`).

2.  **SNV/Indel Processing:**
    *   Restricts SNVs/Indels to chr1-22, X, and Y, drops records with no ALT allele, then decomposes and normalizes them with bcftools.
    *   Annotates variants with echtvar (gnomAD, ClinVar, HPRC, RIMGC, HGMD) and keeps those below the gnomAD popmax, HPRC, and RIMGC frequency thresholds.
    *   Keeps PASS variants overlapping RefSeq genes.
    *   Annotates genic variants with VEP (RefSeq transcripts) and the SpliceAI, REVEL, LoFtool, CADD, and dbNSFP (AlphaMissense, PrimateAI, BayesDel) plugins, plus custom ClinVar and HGMD annotations.
    *   Uses slivar to remove ClinVar benign/likely benign variants and keep variants that are missense or more severe, splice-region, in ClinVar/HGMD, or have a SpliceAI delta score > 0.5.

3.  **Variant Inheritance and Ranking:**
    *   Tags compound heterozygous SNV/Indel pairs using the pedigree (`getcomphets.py`).
    *   Ranks variants (`variantranker_v3.py`) using population frequency (gnomAD, HPRC, RIMGC), pathogenicity predictions (BayesDel, AlphaMissense, PrimateAI, REVEL, CADD, SpliceAI, LoFtool), ClinVar/HGMD annotations, gene constraint (LOEUF, regional missense constraint), disease-gene lists (OMIM morbid genes, GenCC), and de novo status.
    *   Dominant variants are sorted by ClinVar tier (Pathogenic > Likely pathogenic > VUS > other), then by a D_score logistic model fit on ClinVar. Recessive variants are reported separately.
    *   Identifies SNV + SV compound heterozygotes by combining filtered SNVs/Indels with insertions, inversions, and CNVs (`getcomphets_SVs.py`).

4.  **Methylation Analysis:**
    *   Builds a CpG island methylation profile for the sample with methbat.
    *   Intersects the profile with imprinted genes, compares it to a control profile, and reports SNVs/Indels and SVs falling in differentially methylated regions (`get_methylation_output.py`).

5.  **Repeat Expansion Analysis:**
    *   Parses the TRGT VCF against a repeats database and a pathogenic repeats BED file to report potential repeat expansions (`trgt_parser.py`).

6.  **Final Reporting:**
    *   Compiles the results into a single Excel report with one sheet each for: Dominant Variants, Recessive Variants, SNV+SV CompHets, Large Insertions, CNVs, Non-Coding CNVs, Inversions, Variants in DMR, and Repeat Expansion Report.

## Key Files

*   `Snakefile`: The main Snakemake workflow definition file, orchestrating all the analysis steps.
*   `config.yaml`: Central configuration file defining input files, sample information, filtering thresholds, containers, annotation resources, and compute resources. The sample input fields are blank and **must be filled in for each analysis.**
*   `config_template.yaml`: A template for creating the `config.yaml` file, with placeholder paths for the annotation resources.
*   `environment.yaml`: Defines the `wafl_env` Conda environment (Snakemake, bcftools, bedtools, AnnotSV, methbat, SVAFotate, and the Python packages used by `workflow_scripts/`), with exact versions.
*   `snakemake.script.sh`: Runs the workflow with Singularity containers, `/mnt/isilon` bound into the containers, up to 500 parallel jobs, and `--keep-going`.
*   `submissionscript.sh`: An example Slurm submission script that runs `snakemake.script.sh` on a cluster.
*   `workflow_scripts/`: Python scripts used for the annotation, filtering, ranking, and reporting steps within the workflow.

## Required Inputs

To run the WAFL workflow, the following inputs must be specified in the `config.yaml` file:

*   **Sample Information:**
    *   `Sample`: The proband identifier. It must match the proband ID in the pedigree file.
*   **Input Variant Files:**
    *   `SNV_vcf`: Path to the input SNV/Indel VCF file.
    *   `SV_vcf`: Path to the input Structural Variant (SV) VCF file (Sawfish).
    *   `hitsv_vcf`: Path to the HiTSV VCF file, used to match breakends to insertions.
    *   `trgt_vcf`: Path to the TRGT VCF file for repeat expansion analysis.
*   **Pedigree File:**
    *   `ped`: Path to the pedigree file, defining family relationships (proband, parents).
*   **Methylation Input:**
    *   `cpgfile_prefix`: Prefix for the sample's CpG pileup files; the workflow expects `<cpgfile_prefix>.combined.bed.gz`.

The remaining settings point to annotation resources and defaults. In `config_template.yaml`, the resource paths are placeholders that must be set to local copies of each file:

*   **Annotation and Reference Files:** `hg38_ref`, echtvar archives (`gnomad_archive`, `clinvar_archive`, `hprc_archive`, `rimgc_archive`, `hgmd_archive`), VEP resources (`refseq_cache_dir`, `spliceAI_*`, `revel_archive`, `lof_archive`, `cadd_*`, `dbnsfp_file`, `clinvar_processed`, `hgmd_processed`), SV resources (`svafotate_bed`, `svafotatecontrols_bed`, `annotSV_annotations`, `refseq_gtf`, `segduplications`), ranking resources (`loeuf_file`, `morbidgenes_file`, `gencc_file`, `regionalconstraint_file`), methylation resources (`cpg_island_regions`, `controlmethylationprofile`, `imprintedgenesbed`), and repeat resources (`repeatsdb`, `pathogenicrepeats_bed`). Refer to `config.yaml` for the complete list.
*   **Filtering Thresholds:**
    *   `gnomad_popmax_af_threshold`, `hprc_FRQ_threshold`, `rimgc_FRQ_threshold`: Population frequency thresholds for SNVs/Indels.
    *   `gnomad_cutoff`, `hprc_cutoff`, `colorsdb_cutoff`, `cdls_cutoff`: Allele frequency cutoffs for SV filtering.
*   **Tool Containers:** Singularity/Docker images for echtvar, VEP, and slivar.
*   **Resources:** `memory`, `cpu`, `svafotate_mem`, `vep_cpu`, `vep_mem`.

## Setup and Execution

1.  **Environment Setup:**
    Create and activate the Conda environment using the provided `environment.yaml` file:
    ```bash
    conda env create -f environment.yaml
    conda activate wafl_env
    ```
    *(Note: Ensure you have Conda installed and configured. Singularity must also be available to run the containerized steps.)*

2.  **Configuration:**
    *   Copy `config_template.yaml` to `config.yaml`.
    *   Edit `config.yaml` to fill in the sample ID, the paths to the sample's input files listed under **Required Inputs**, and the paths to the annotation resources.

3.  **Workflow Execution:**
    Run the workflow from the directory containing the `Snakefile`:
    *   **Dry run (recommended first):**
        ```bash
        snakemake -n --configfile config.yaml
        ```
    *   **Full run:**
        ```bash
        bash snakemake.script.sh
        ```
    *   **Cluster run (Slurm):**
        ```bash
        bash submissionscript.sh
        ```
        *(Adjust the partition, time, and memory in `submissionscript.sh` for your cluster.)*
    *   **Workflow diagram:**
        ```bash
        snakemake --dag | dot -Tsvg > dag.svg
        ```

## Output

All outputs are written to `output/[SampleID]/`. The primary output is `[SampleID].CompiledReport.xlsx`, which consolidates the results from all analysis modules. Intermediate files are kept alongside it, including:

*   `[SampleID].DominantVariants.tsv`, `[SampleID].RecessiveVariants.tsv`, `[SampleID].Features.tsv`: Ranked SNVs/Indels and the features used to rank them.
*   `[SampleID].Filtered.CompHetsTagged.vcf(.gz)`: Filtered SNVs/Indels with compound heterozygote tags.
*   `[SampleID].FilteredSVs.tsv`, `[SampleID].FilteredInversions.tsv`, `[SampleID].filtered_cnv.tsv`, `[SampleID].noncoding_cnv.tsv`, `[SampleID].all_insertions.tsv`, `[SampleID].bnd_ins_matches.tsv`: Filtered SV results.
*   `[SampleID].snv_sv_comphets.tsv`: SNV + SV compound heterozygotes.
*   `[SampleID].methbatprofile.tsv`, `[SampleID].VarsInDMRs.tsv`: Methylation profile and variants in DMRs.
*   `[SampleID].repeat_expansion_report.tsv`: Repeat expansion report.
