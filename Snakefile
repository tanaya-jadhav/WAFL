# WAFL: Workflow for Annotating & Filtering Long Read Data - Tanaya Jadhav
# Run in conda environment wafl_env
# To make new conda environment use environment.yaml
# To generate a DAG run: snakemake --dag | dot -Tsvg > dag.svg

configfile: "config.yaml"

# Define the main rule
rule all:
    input:  expand(["output/{sample}/{sample}.CompiledReport.xlsx"], sample=config["Sample"])

# Filter extra contigs to only keep chr1-22,X,Y. Decompose & Normalize filtered vcf
# In Case of Error: [W::vcf_parse_info] INFO 'DB' is not defined in the header, assuming Type=String
# [E::vcf_format] Invalid BCF, the INFO tag id=54 is too large at 
# run the following prior to re-running this workflow: bcftools annotate -h templine.txt <inputfile> > <outputvcffile>
# then bgzip vcf file and use as new input for workflow

# SV filtering module
# annotate with gnomad, 1000 genomes, and other public data allele frequencies
rule svafotate_annotate:
    input:
        vcf=config["SV_vcf"]
    output:
        annotatedvcf="output/{sample}/{sample}.svafotate.vcf"
    params:
        svafotate_bed=config["svafotate_bed"]
    resources:
            mem=config["svafotate_mem"]
    threads:
            config["cpu"]
    shell:
        """
        svafotate annotate \
        -v {input.vcf} \
        -o {output.annotatedvcf} \
        -b {params.svafotate_bed} \
        -f 0.5 -a best
        """

rule svafotate_platformcontrolannotation:
    input:
        vcf="output/{sample}/{sample}.svafotate.vcf"
    output:
        out_vcf="output/{sample}/{sample}.svafotate.withplatformcontrols.vcf"
    params:
        svafotatecontrols_bed=config["svafotatecontrols_bed"]
    shell:
        """
        svafotate annotate \
        -v {input.vcf} \
        -o {output.out_vcf} \
        -b {params.svafotatecontrols_bed} \
        -f 0.5 -a best
        """

# annotate with vep
rule vep_annotate_SVs:
    input:
        VEP_input_VCF="output/{sample}/{sample}.svafotate.withplatformcontrols.vcf"
    output:
        VEP_output_VCF="output/{sample}/{sample}.svafotate.vep.vcf"
    params:
        refseq_cache_dir=config["refseq_cache_dir"],
        revel_archive=config["revel_archive"],
        lof_archive=config["lof_archive"],
        cadd_svarchive=config["cadd_svarchive"],
        hg38_ref=config["hg38_ref"]
    log:
        "output/{sample}/{sample}.vep.log"
    container:
        config["VEP_container"]
    resources:
        mem_mb=config["vep_mem"]
    threads:
        config["vep_cpu"]
    shell:
        """
        vep \
        --verbose \
        --no_stats \
        --offline \
        --dir {params.refseq_cache_dir} \
        --cache_version 112 \
        --species homo_sapiens \
        --assembly GRCh38 \
        --force_overwrite \
        --fork 16 \
        --buffer_size 20000 \
        --vcf \
        --refseq \
        --fasta {params.hg38_ref} \
        --exclude_predicted \
        --distance 0 \
        --minimal \
        --allele_number \
        --hgvs \
        --pick \
        -o {output.VEP_output_VCF} \
        -i {input.VEP_input_VCF} \
        --plugin REVEL,{params.revel_archive} \
        --plugin LoFtool,{params.lof_archive} \
        --plugin CADD,sv={params.cadd_svarchive}
        """

rule filter_SVs:
    input:
        input_VCF="output/{sample}/{sample}.svafotate.vep.vcf"
    output:
        output_VCF="output/{sample}/{sample}.FilteredSVs.vcf",
        output_TSV="output/{sample}/{sample}.FilteredSVs.vcf.tsv"
    params:
        gnomad_cutoff=config["gnomad_cutoff"],
        hprc_cutoff=config["hprc_cutoff"],
        colorsdb_cutoff=config["colorsdb_cutoff"]
    shell:
        """
        ./workflow_scripts/filterSVs.py \
        -in_vcf {input.input_VCF} \
        -gnomad_cutoff {params.gnomad_cutoff} \
        -hprc_cutoff {params.hprc_cutoff} \
        -colorsdb_cutoff {params.colorsdb_cutoff} \
        -o_vcf {output.output_VCF} \
        -o_tsv {output.output_TSV}
        """

# SNV Indel Filtering Module 
rule bcftools_filterAndNormalize:
    input:
        vcf=config["SNV_vcf"],
        # vcf="inputs/{sample}.hard-filtered.vcf.gz",
        ref="/mnt/isilon/rajagopalan_lab/reference/Hg38.fasta"
    output:
        filteredvcf="output/{sample}/{sample}.extracontigsfiltered.vcf",
        normalizedvcf="output/{sample}/{sample}.normalized.vcf"
    shell:
        """
        bcftools view {input.vcf} --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY -O v -o {output.filteredvcf} --threads 10
        bcftools norm -m - {output.filteredvcf} -w 10000 -f {input.ref} -O v -o {output.normalizedvcf} --threads 10
        """

# Reformat to replace empty fields for echtvar
rule reformatvcf:
    input:
        vcf="output/{sample}/{sample}.normalized.vcf"
    output:
        vcf="output/{sample}/{sample}.normalized.reformatted.vcf"
    shell:
        """
        sed -e "s/\.;//g" -e "s/Flag/String/g" {input.vcf} > {output.vcf}
        """

# Annotate variants with gnomad, clinvar, hprc databases and filter based on population frequency
rule echtvar_annotate:
    input:
        ECHTVAR_input_VCF="output/{sample}/{sample}.normalized.reformatted.vcf"
    output:
        ECHTVAR_output_VCF="output/{sample}/{sample}.echtvar_annotated.vcf.gz",
        ECHTVAR_output_index="output/{sample}/{sample}.echtvar_annotated.vcf.gz.csi"
    params:
        gnomad_archive=config["gnomad_archive"],
        clinvar_archive=config["clinvar_archive"],
        hprc_archive=config["hprc_archive"],
        gnomad_popmax_af_threshold=config["gnomad_popmax_af_threshold"],
        hprc_FRQ_threshold=config["hprc_FRQ_threshold"]
    log:
        "output/{sample}/{sample}.echtvar.log"
    container: 
        config["ECHTVAR_container"]
    resources:
        mem_mb=config["memory"]
    threads:
        config["cpu"]
    shell:
        """
        echtvar anno \
        -e {params.gnomad_archive} \
        -e {params.clinvar_archive} \
        -e {params.hprc_archive} \
        -i "(gnomad_popmax_af < {params.gnomad_popmax_af_threshold} && hprc_FRQ < {params.hprc_FRQ_threshold})" \
        {input.ECHTVAR_input_VCF} \
        {output.ECHTVAR_output_VCF}
        bcftools index {output.ECHTVAR_output_VCF}
        """

# Filter against BED file from refseq to keep only genic variants
rule bcftools_genicfilter:
    input:
        GENIC_input_vcf="output/{sample}/{sample}.echtvar_annotated.vcf.gz",
        GENIC_input_index="output/{sample}/{sample}.echtvar_annotated.vcf.gz.csi"
    output:
        GENIC_output_vcf="output/{sample}/{sample}.GenicVariantsOnly.vcf"
    params:
        refseq_bed=config["refseq_genes_bed"]
    resources:
        mem_mb=config["memory"]
    threads:
        config["cpu"]
    shell:
        """
        bcftools view -R {params.refseq_bed} -f PASS {input.GENIC_input_vcf} -O v -o {output.GENIC_output_vcf}
        """

# Annotate genic variants with VEP 
rule vep_annotate:
    input:
        VEP_input_VCF="output/{sample}/{sample}.GenicVariantsOnly.vcf"
    output:
        VEP_output_VCF="output/{sample}/{sample}.vep_annotated.vcf"
    params:
        refseq_cache_dir=config["refseq_cache_dir"],
        spliceAI_snvarchive=config["spliceAI_snvarchive"],
        spliceAI_indelarchive=config["spliceAI_indelarchive"],
        revel_archive=config["revel_archive"],
        lof_archive=config["lof_archive"],
        cadd_snvarchive=config["cadd_snvarchive"],
        cadd_indelarchive=config["cadd_indelarchive"],
        clinvar_processed=config["clinvar_processed"],
        hg38_ref=config["hg38_ref"]
    log:
        "output/{sample}/{sample}.vep.log"
    container:
        config["VEP_container"]
    resources:
        mem_mb=config["vep_mem"]
    threads:
        config["vep_cpu"]
    shell:
        """
        vep \
        --verbose \
        --no_stats \
        --offline \
        --dir {params.refseq_cache_dir} \
        --cache_version 112 \
        --species homo_sapiens \
        --assembly GRCh38 \
        --force_overwrite \
        --fork 16 \
        --buffer_size 20000 \
        --vcf \
        --refseq \
        --fasta {params.hg38_ref} \
        --exclude_predicted \
        --distance 0 \
        --minimal \
        --allele_number \
        --hgvs \
        -o {output.VEP_output_VCF} \
        -i {input.VEP_input_VCF} \
        --plugin SpliceAI,snv={params.spliceAI_snvarchive},indel={params.spliceAI_indelarchive} \
        --plugin REVEL,{params.revel_archive} \
        --plugin LoFtool,{params.lof_archive} \
        --plugin CADD,snv={params.cadd_snvarchive},indels={params.cadd_indelarchive} \
        --custom file={params.clinvar_processed},short_name=CLINVAR,format=vcf,type=overlap,fields=Clinvar_Pathogenic%Clinvar_Benign%Clinvar_Uncertain.significance
        """

# filter benign/likely benign variants
rule slivar_filter:
    input:
        SLIVAR_input_VCF="output/{sample}/{sample}.vep_annotated.vcf"
    output:
        SLIVAR_output_VCF="output/{sample}/{sample}.SlivarHighestImpactClinVarHGMD.vcf"
    params:
        csq_script=config["csq_script"]
    container:
        config['SLIVAR_container']
    resources:
        mem_mb=config["memory"]
    threads:
        config["cpu"]
    shell:
        """
        slivar expr --js {params.csq_script} \
        --vcf {input.SLIVAR_input_VCF} \
        --out-vcf {output.SLIVAR_output_VCF} \
        --info "!(INFO.CLNSIG == 'Benign' || INFO.CLNSIG == 'Benign/Likely_benign' || INFO.CLNSIG == 'Likely_benign') && (INFO.highest_impact_order <= ImpactOrder.missense || !(INFO.CLNSIG == 'MISSING' || INFO.PHEN == 'MISSING') || CSQs(INFO.CSQ, VCF.CSQ, ['SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL']).some(function(csq) {{ return csq.SPLICEAI_PRED_DS_AG > 0.5 || csq.SPLICEAI_PRED_DS_AL > 0.5 || csq.SPLICEAI_PRED_DS_DG > 0.5 || csq.SPLICEAI_PRED_DS_DL > 0.5 }}))"
        """


# Add a tag in INFO field for compound heterozygous variants
rule comphets:
    input:
        COMPHET_input_VCF="output/{sample}/{sample}.SlivarHighestImpactClinVarHGMD.vcf",
        COMPHET_input_PED=config["ped"]
    output:
        COMPHET_output_VCF="output/{sample}/{sample}.Filtered.CompHetsTagged.vcf",
        COMPHET_output_gzippedVCF="output/{sample}/{sample}.Filtered.CompHetsTagged.vcf.gz"
    resources:
        mem_mb=config["memory"]
    threads:
        config["cpu"]
    shell:
        """
        ./workflow_scripts/getcomphets.py \
        -i {input.COMPHET_input_VCF} \
        -o {output.COMPHET_output_VCF} \
        -ped {input.COMPHET_input_PED} \
        -proband {wildcards.sample}
        grep "^#" {output.COMPHET_output_VCF} > output/{wildcards.sample}/temp.vcf
        grep -v "^#" {output.COMPHET_output_VCF} | sort -k1,1V -k2,2g >> output/{wildcards.sample}/temp.vcf
        bcftools view -Oz output/{wildcards.sample}/temp.vcf -o {output.COMPHET_output_gzippedVCF}
        bcftools index {output.COMPHET_output_gzippedVCF}
        rm output/{wildcards.sample}/temp.vcf
        """

# rank dominant and recessive variants
rule rank_variants:
    input:
        Rank_input_VCF="output/{sample}/{sample}.Filtered.CompHetsTagged.vcf",
        Rank_input_PED=config["ped"]
    output:
        Rank_output_Features="output/{sample}/{sample}.Features.tsv",
        Rank_output_Dominant="output/{sample}/{sample}.DominantVariants.tsv",
        Rank_output_Recessive="output/{sample}/{sample}.RecessiveVariants.tsv"
    params:
        loeuf_file=config["loeuf_file"],
        morbidgenes_file=config["morbidgenes_file"],
        gencc_file=config["gencc_file"],
        regionalconstraint_file=config["regionalconstraint_file"]
    resources:
        mem_mb=config["memory"]
    threads:
        config["cpu"]
    shell:
        """
        ./workflow_scripts/variantranker.py \
        -InVCF {input.Rank_input_VCF} \
        -ped {input.Rank_input_PED} \
        -proband {wildcards.sample} \
        -loeuf {params.loeuf_file} \
        -morbidgenes {params.morbidgenes_file} \
        -gencc {params.gencc_file} \
        -regconst {params.regionalconstraint_file} \
        -OutDir output/{wildcards.sample}/
        """

# create sample methylation profile
rule create_methprofile:
    input:
        f"{config['cpgfile_prefix']}.combined.bed.gz"
    output:
        profile_TSV="output/{sample}/{sample}.methbatprofile.tsv"
    params:
        cpgfile_prefix=config["cpgfile_prefix"],
        cpg_island_regions=config["cpg_island_regions"]
    resources:
        mem_mb=config["memory"]
    threads:
        config["cpu"]
    shell:
        """
        methbat profile \
        --input-prefix {params.cpgfile_prefix} \
        --input-regions {params.cpg_island_regions} \
        --output-region-profile {output.profile_TSV}
        """

rule find_DMR_Variants:
    input:
        profile_TSV="output/{sample}/{sample}.methbatprofile.tsv",
        snv_VCF="output/{sample}/{sample}.Filtered.CompHetsTagged.vcf",
        sv_VCF=get_input_file
    output:
        out_imprintedprofile="output/{sample}/{sample}.imprintedgenesmethylation.tsv",
        out_TSV="output/{sample}/{sample}.VarsInDMRs.tsv"
    params:
        controlmethylationprofile=config["controlmethylationprofile"],
        imprintedgenesbed=config["imprintedgenesbed"]
    shell:
        """
        bedtools intersect -a {input.profile_TSV} -b {params.imprintedgenesbed} -wo > {output.out_imprintedprofile}
        ./workflow_scripts/get_methylation_output.py \
        -case {output.out_imprintedprofile} \
        -control {params.controlmethylationprofile} \
        -sv {input.sv_VCF} \
        -snv {input.snv_VCF} \
        -o {output.out_TSV}
        rm temp_differentiallymethylatedregions.bed
        """

# Combine SNVs and SVs
rule SV_SNV_comphets:
    input:
        SNV_file="output/{sample}/{sample}.Filtered.CompHetsTagged.vcf",
        SV_file="output/{sample}/{sample}.FilteredSVs.vcf",
        input_PED=config["ped"]
    output:
        OUT_tsv="output/{sample}/{sample}.snv_sv_comphets.tsv"
    params:
        loeuf_file=config["loeuf_file"],
        morbidgenes_file=config["morbidgenes_file"],
        gencc_file=config["gencc_file"]
    resources:
        mem_mb=config["memory"]
    threads:
        config["cpu"]
    shell:
        """
        ./workflow_scripts/getcomphets_SVs.py \
        -sv {input.SV_file} \
        -snv {input.SNV_file} \
        -o {output.OUT_tsv} \
        -loeuf {params.loeuf_file} \
        -morbidgenes {params.morbidgenes_file} \
        -gencc {params.gencc_file} \
        -proband {wildcards.sample} \
        -ped {input.input_PED}
        """

# run repeat expansion workflow
rule find_repeatexpansions:
    input:
        trgt_VCF=config["trgt_vcf"]
    output:
        repeat_TSV="output/{sample}/{sample}.repeat_expansion_report.tsv"
    params:
        repeatsdb=config["repeatsdb"],
        pathogenicrepeats_bed=config["pathogenicrepeats_bed"]
    shell:
        """
        bcftools view {input.trgt_VCF} -O v -o temptrgtfile.vcf
        ./workflow_scripts/trgt_parser.py \
        -trgt_vcf temptrgtfile.vcf \
        -repeatsdb {params.repeatsdb} \
        -bed {params.pathogenicrepeats_bed} \
        -o {output.repeat_TSV}
        rm temptrgtfile.vcf
        """

# Output Module

# write output excel file
# if you would like to delete intermediate files, you can add this code to the following rule
rule write_excel:
    input:
        D_TSV="output/{sample}/{sample}.DominantVariants.tsv",
        R_TSV="output/{sample}/{sample}.RecessiveVariants.tsv",
        snv_sv="output/{sample}/{sample}.snv_sv_comphets.tsv",
        SV_TSV="output/{sample}/{sample}.FilteredSVs.vcf.tsv",
        dmr="output/{sample}/{sample}.VarsInDMRs.tsv",
        tr="output/{sample}/{sample}.repeat_expansion_report.tsv"
    output:
        out_xl="output/{sample}/{sample}.CompiledReport.xlsx"
    shell:
        """
        ./workflow_scripts/write_excel.py \
        -D_tsv {input.D_TSV} \
        -R_tsv {input.R_TSV} \
        -snv_sv {input.snv_sv} \
        -SV {input.SV_TSV} \
        -DMR {input.dmr} \
        -TR {input.tr} \
        -o {output.out_xl}
        """
