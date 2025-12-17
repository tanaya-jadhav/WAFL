#!/usr/bin/env python

import pandas as pd
import numpy as np
from peds import open_ped
import sys
import argparse


def read_pedfile(in_ped, proband_ID):
    families = open_ped(in_ped)
    family = families[0]
    proband = family.get_proband()
    # print(proband)
    if proband.id == proband_ID:
        father = family.get_father(proband)
        mother = family.get_mother(proband)
        print('Input family from ped file\nProband: ', proband.id, 'Father: ', father, 'Mother: ', mother)
        return proband.id, father, mother
    else:
        print("Proband in ped file does not match input proband\n" + proband.id + '≠' + proband_ID + '\nExiting...')
        sys.exit(1)


def maxspliceaiscore(colAG, colAL, colDG, colDL):
    if len({colAG, colAL, colDG, colDL}) == 1:
        return np.NAN
    elif max(colAG, colAL, colDG, colDL) == colAG:
        return 'Delta score for acceptor gain'
    elif max(colAG, colAL, colDG, colDL) == colAL:
        return 'Delta score for acceptor loss'
    elif max(colAG, colAL, colDG, colDL) == colDG:
        return 'Delta score for donor gain'
    elif max(colAG, colAL, colDG, colDL) == colDL:
        return 'Delta score for donor loss'
    else:
        return np.NAN


def getregionalconsscore(regionalconstraint_df, varchrom, varpos):
    varchrom = varchrom.split('chr')[1]
    try:
        vardf = regionalconstraint_df[(regionalconstraint_df['hg38_chr'] == varchrom) &
                                      (regionalconstraint_df['hg38_startpos'] <= varpos) &
                                      (regionalconstraint_df['hg38_endpos'] >= varpos)]
        score = list(vardf['oe'])[0]
        # print(varchrom, varpos)
        # print(vardf)
        return score
    except:
        return 0


def createfeaturedf(proband_ID, in_df, sampletype, outdir, loeuf_df, morbidgenes_df, gencc_df, regionalconstraint_df):
    # pull out features for scoring
    in_df = in_df.rename(columns={"#CHROM": "CHROM"})
    in_df['CSQ'] = in_df['INFO'].str.split('CSQ=').str[1].str.split(',').str[0].str.split(';').str[0]
    in_df['Gene'] = in_df['CSQ'].str.split('|').str[3]
    in_df['Exon'] = in_df['CSQ'].str.split('|').str[8]
    in_df['Consequence'] = in_df['CSQ'].str.split('|').str[1]
    in_df['Impact'] = in_df['CSQ'].str.split('|').str[2]
    in_df['Highest_Impact_Order'] = in_df['INFO'].str.split('highest_impact_order=').str[1].str.split(';').str[0]
    in_df['Highest_Impact_Order'] = pd.to_numeric(in_df['Highest_Impact_Order'], errors='coerce')
    in_df['HGVSc'] = in_df['CSQ'].str.split('|').str[10]
    in_df['HGVSp'] = in_df['CSQ'].str.split('|').str[11]
    in_df['gnomad_popmax_af'] = in_df['INFO'].str.split('gnomad_popmax_af=').str[1].str.split(';').str[0]
    in_df['hprc_FRQ'] = in_df['INFO'].str.split('hprc_FRQ=').str[1].str.split(';').str[0]
    in_df['rimgc_FRQ'] = in_df['INFO'].str.split('rimgc_FRQ=').str[1].str.split(';').str[0]
    in_df['Hgmd_class'] = in_df['INFO'].str.split('HGVS_CLASS=').str[1].str.split(';').str[0]
    in_df['Hgmd_rankscore'] = in_df['INFO'].str.split('HGMD_RANKSCORE=').str[1].str.split(';').str[0]
    in_df['CLNSIG'] = in_df['INFO'].str.split('CLNSIG=').str[1].str.split(';').str[0]
    in_df['CLNDN'] = in_df['INFO'].str.split('CLNDN=').str[1].str.split(';').str[0]
    in_df['CompHet'] = np.where(in_df['INFO'].str.contains('Compound_Heterozygous'), 1, 0)
    in_df['Homozygous'] = np.where((in_df['Proband_GT'] == '1|1') | (in_df['Proband_GT'] == 1 / 1), 1, 0)
    in_df['SpliceAI_pred_DS_AG'] = in_df['CSQ'].str.split('|').str[37]
    in_df['SpliceAI_pred_DS_AG'] = pd.to_numeric(in_df['SpliceAI_pred_DS_AG'], errors='coerce')
    in_df['SpliceAI_pred_DS_AL'] = in_df['CSQ'].str.split('|').str[38]
    in_df['SpliceAI_pred_DS_AL'] = pd.to_numeric(in_df['SpliceAI_pred_DS_AL'], errors='coerce')
    in_df['SpliceAI_pred_DS_DG'] = in_df['CSQ'].str.split('|').str[39]
    in_df['SpliceAI_pred_DS_DG'] = pd.to_numeric(in_df['SpliceAI_pred_DS_DG'], errors='coerce')
    in_df['SpliceAI_pred_DS_DL'] = in_df['CSQ'].str.split('|').str[40]
    in_df['SpliceAI_pred_DS_DL'] = pd.to_numeric(in_df['SpliceAI_pred_DS_DL'], errors='coerce')
    in_df['Revel'] = in_df['CSQ'].str.split('|').str[42]
    in_df['LoFtool'] = in_df['CSQ'].str.split('|').str[43]
    in_df['CADD_PHRED'] = in_df['CSQ'].str.split('|').str[44]
    in_df['CADD_PHRED'] = pd.to_numeric(in_df['CADD_PHRED'], errors='coerce')
    in_df['CADD_RAW'] = in_df['CSQ'].str.split('|').str[45]
    in_df['Clinvar_path'] = in_df['CSQ'].str.split('|').str[47]
    in_df['Clinvar_path'] = pd.to_numeric(in_df['Clinvar_path'], errors='coerce')
    in_df['Clinvar_benign'] = in_df['CSQ'].str.split('|').str[48]
    in_df['Clinvar_benign'] = pd.to_numeric(in_df['Clinvar_benign'], errors='coerce')
    in_df['Clinvar_vus'] = in_df['CSQ'].str.split('|').str[49]
    in_df['Clinvar_vus'] = pd.to_numeric(in_df['Clinvar_vus'], errors='coerce')
    in_df['SpliceAI_maxscore'] = in_df[
        ['SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL']].max(axis=1)
    # v = np.vectorize(maxspliceaiscore)
    # in_df['SpliceAI_maxtype'] = v(in_df.SpliceAI_pred_DS_AG, in_df.SpliceAI_pred_DS_AL, in_df.SpliceAI_pred_DS_DG,
    #                               in_df.SpliceAI_pred_DS_DL)
    in_df['SpliceAI_maxtype'] = in_df.apply(lambda x: maxspliceaiscore(x.SpliceAI_pred_DS_AG, x.SpliceAI_pred_DS_AL,
                                                                       x.SpliceAI_pred_DS_DG, x.SpliceAI_pred_DS_DL), axis=1)
    # Is variant de novo?
    if sampletype == 'trio':
        in_df['Denovo'] = np.where((((in_df['Dad_GT'] == '0|0') | (in_df['Dad_GT'] == '0/0')) &
                                    ((in_df['Mom_GT'] == '0|0') | (in_df['Mom_GT'] == '0/0'))), 1, 0)
    elif (sampletype == 'duo_dad') | (sampletype == 'duo_mom'):
        in_df['Denovo'] = np.where((in_df['Parent_GT'] == '0|0') | (in_df['Parent_GT'] == '0/0'), 1, 0)
    else:
        in_df['Denovo'] = 0

    # add loeuf score
    in_df = in_df.merge(loeuf_df, left_on='Gene',  right_on='gene', how='left')
    in_df = in_df.rename(columns={"lof.oe": "Loeuf_score"})

    # add morbid gene annotation
    in_df = in_df.merge(morbidgenes_df, left_on='Gene', right_on='gene', how='left')

    # add gencc_morbid gene annotation
    in_df = in_df.merge(gencc_df, left_on='Gene', right_on='Gene', how='left')

    # add regional missense constraint score
    in_df['regional_missense_constraint'] = in_df.apply(lambda  x: getregionalconsscore(regionalconstraint_df, x.CHROM, x.POS), axis=1)

    if sampletype == 'trio':
        feature_df = in_df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'Proband_GT', 'Gene', 'Consequence', 'Impact',
                            'Highest_Impact_Order', 'HGVSc',
                            'HGVSp', 'gnomad_popmax_af', 'Hgmd_class', 'Hgmd_rankscore', 'CLNSIG', 'CLNDN', 'CompHet',
                            'Homozygous', 'Denovo', 'Dad_GT', 'Mom_GT', 'Revel', 'LoFtool', 'Loeuf_score', 'regional_missense_constraint',
                            'SpliceAI_maxscore', 'SpliceAI_maxtype', 'CADD_PHRED', 'CADD_RAW', 'Morbid_Gene', 'GenCC_Submission',
                            'Clinvar_path', 'Clinvar_benign', 'Clinvar_vus', 'gnomad_popmax_af', 'hprc_FRQ', 'rimgc_FRQ']]
    elif (sampletype == 'duo_dad') | (sampletype == 'duo_mom'):
        feature_df = in_df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'Proband_GT', 'Gene', 'Consequence', 'Impact',
                            'Highest_Impact_Order', 'HGVSc',
                            'HGVSp', 'gnomad_popmax_af', 'Hgmd_class', 'Hgmd_rankscore', 'CLNSIG', 'CLNDN', 'CompHet',
                            'Homozygous', 'Denovo', 'Parent_GT', 'Revel', 'LoFtool', 'Loeuf_score', 'regional_missense_constraint', 'SpliceAI_maxscore',
                            'SpliceAI_maxtype', 'CADD_PHRED', 'CADD_RAW', 'Morbid_Gene', 'GenCC_Submission', 'Clinvar_path',
                            'Clinvar_benign', 'Clinvar_vus', 'gnomad_popmax_af', 'hprc_FRQ', 'rimgc_FRQ']]
    else:
        feature_df = in_df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'Proband_GT', 'Gene', 'Consequence', 'Impact',
                            'Highest_Impact_Order', 'HGVSc',
                            'HGVSp', 'gnomad_popmax_af', 'Hgmd_class', 'Hgmd_rankscore', 'CLNSIG', 'CLNDN', 'CompHet',
                            'Homozygous', 'Denovo', 'Revel', 'LoFtool', 'Loeuf_score', 'regional_missense_constraint', 'SpliceAI_maxscore',
                            'SpliceAI_maxtype', 'CADD_PHRED', 'CADD_RAW', 'Morbid_Gene', 'GenCC_Submission', 'Clinvar_path',
                            'Clinvar_benign', 'Clinvar_vus', 'gnomad_popmax_af', 'hprc_FRQ', 'rimgc_FRQ']]
    feature_df['LoFtool'] = feature_df.LoFtool.replace('', 0, regex=True).astype(float)
    feature_df['SpliceAI_maxscore'] = feature_df.SpliceAI_maxscore.replace('', 0, regex=True).astype(float)
    feature_df['Revel'] = feature_df.Revel.replace('', 0, regex=True).astype(float)
    feature_df[['SpliceAI_maxscore', 'Morbid_Gene', 'Clinvar_path', 'Clinvar_benign', 'Clinvar_vus']] = feature_df[['SpliceAI_maxscore', 'Morbid_Gene', 'Clinvar_path', 'Clinvar_benign', 'Clinvar_vus']].fillna(value=0).astype(float)

    clinvar_conds = [feature_df['CLNSIG'].isin(['Pathogenic', 'Pathogenic/Likely_pathogenic', 'Pathogenic/Likely_pathogenic|risk_factor']),
             feature_df['CLNSIG'].isin(['Likely_pathogenic']),
             feature_df['CLNSIG'].isin(['Uncertain_significance', 'Conflicting_classifications_of_pathogenicity'])]
    clinvar_values = [3, 2, 1]
    feature_df['Clinvar_sig'] = np.select(clinvar_conds, clinvar_values, default=0)

    hgmd_conds = [feature_df['Hgmd_class'].isin(['DM']),
                     feature_df['Hgmd_class'].isin(['DM?'])]
    hgmd_values = [2, 1]
    feature_df['hgmd_class'] = np.select(hgmd_conds, hgmd_values, default=0)

    impact_conds = [feature_df['Impact'].isin(['HIGH']),
                    feature_df['Impact'].isin(['MODERATE']),
                    feature_df['Impact'].isin(['LOW', 'MODIFIER'])]
    impact_values = [3, 2, 1]
    feature_df['Impact_significance'] = np.select(impact_conds, impact_values, default=0)

    feature_df.to_csv(outdir + proband_ID + '.Features.tsv', sep='\t', index=False)
    return feature_df


def main(in_vcf, pedfile, proband_ID, outdir, loeuf_file, morbidgenes_file, gencc_file, regionalconstraint_file):
    # get proband, mom, dad from ped file
    proband, father, mother = read_pedfile(pedfile, proband_ID)

    # parse vcf file
    with open(in_vcf, "r") as i:
        vcf_header = i.readlines()
    vcf_header = [line for line in vcf_header if line.startswith('#')]
    headerline = [line for line in vcf_header if line.startswith('#CHROM')][0]
    headerline = headerline.strip('\n').split('\t')

    # determine if parent samples are present
    if (father is not None) and (father.id in headerline) and (mother is not None) and (mother.id in headerline):
        sampletype = 'trio'
    elif (father is not None) and (father.id in headerline) and (mother.id not in headerline):
        sampletype = 'duo_dad'
    elif (mother is not None) and (mother.id in headerline) and (father.id not in headerline):
        sampletype = 'duo_mom'
    else:
        sampletype = 'singleton'
    print('Sample is a ' + sampletype)
    pd.set_option('display.max_columns', 500)
    in_df = pd.read_csv(in_vcf, sep='\t', comment='#', header=None, names=headerline)

    # pull genotypes out as columns for analysis
    in_df['Proband_GT'] = in_df[proband_ID].str.split(':').str[0]
    if sampletype == 'trio':
        in_df['Dad_GT'] = in_df[father.id].str.split(':').str[0]
        in_df['Mom_GT'] = in_df[mother.id].str.split(':').str[0]
    elif sampletype == 'duo_dad':
        in_df['Parent_GT'] = in_df[father.id].str.split(':').str[0]
    elif sampletype == 'duo_mom':
        in_df['Parent_GT'] = in_df[mother.id].str.split(':').str[0]

    loeuf_df = pd.read_csv(loeuf_file, sep='\t', header=0)

    morbidgenes_df = pd.read_csv(morbidgenes_file,sep='\t', header=None, names=['gene', 'Morbid_Gene'])

    gencc_df = pd.read_csv(gencc_file, sep='\t', header=0)

    regionalconstraint_df = pd.read_csv(regionalconstraint_file, sep='\t', header=0)
    # print(regionalconstraint_df.head(5))

    # get features in dataframe
    feature_df = createfeaturedf(proband_ID, in_df, sampletype, outdir, loeuf_df, morbidgenes_df, gencc_df, regionalconstraint_df)
    print(feature_df.shape, feature_df.head(5))

    # separate by dominant and recessive
    recessive_df = feature_df[(feature_df['CompHet'] == 1) | (feature_df['Homozygous'] == 1)]
    dominant_df = feature_df[(feature_df['CompHet'] == 0) & (feature_df['Homozygous'] == 0)]

    # remove variants that are homozygous in parents
    if sampletype == 'trio':
        dominant_df = dominant_df[(dominant_df['Dad_GT'] != '1/1') | (dominant_df['Dad_GT'] != '1|1') |
                                    (dominant_df['Mom_GT'] != '1/1') | (dominant_df['Mom_GT'] != '1|1')]
        recessive_df = recessive_df[(recessive_df['Dad_GT'] != '1/1') | (recessive_df['Dad_GT'] != '1|1') |
                                    (recessive_df['Mom_GT'] != '1/1') | (recessive_df['Mom_GT'] != '1|1')]
    elif (sampletype == 'duo_dad') | (sampletype == 'duo_mom'):
        dominant_df = dominant_df[(dominant_df['Parent_GT'] != '1/1') | (dominant_df['Parent_GT'] != '1|1')]
        recessive_df = recessive_df[(recessive_df['Parent_GT'] != '1/1') | (recessive_df['Parent_GT'] != '1|1')]

    # Scoring dominant variants
    # Variant gets a higher score if it is denovo, is in a LoF gene (continuous linear scale), is predicted to affect splicing
    # Scoring is boosted if it is annotated P/LP/VUS in CLinvar and DM in hgmd and in the morbid genes list
    dominant_df['D_score'] = 2 * feature_df['Denovo'] + (1 - feature_df['Loeuf_score']) * feature_df['SpliceAI_maxscore'] + \
                             10 * feature_df['Clinvar_sig'] * feature_df['hgmd_class'] * (feature_df['Morbid_Gene']*2) + \
                             feature_df['Impact_significance'] * (feature_df['Revel'] + feature_df['CADD_PHRED'])

    # Tag variants for prioritization
    # all missense and deletion variants with regional constraint score <0.6
    missensevars = ['missense_variant', 'frameshift_variant', 'inframe_deletion']
    dominant_df.loc[(dominant_df['Consequence'].isin(missensevars)) & (dominant_df['regional_missense_constraint'] < 0.6)
                    & (dominant_df['regional_missense_constraint'] > 0), 'Variant_in_constrained_region'] = 1

    # write dominant df to file
    dominant_df = dominant_df.sort_values(by='D_score', ascending=False)
    dominant_df.to_csv(outdir + proband_ID + '.DominantVariants.tsv', sep='\t', index=False)

    # scoring recessive variants
    genelist = set(recessive_df['Gene'])
    recessive_df['retain'] = ''
    for gene in genelist:
        gene_df = recessive_df[recessive_df['Gene'] == gene]
        # if at least 1 variant is path or likely path, retain this variant and all trans variants in this gene (if phased)
        gene_df_pathvars = gene_df[gene_df['Clinvar_sig'] >= 2]
        if len(gene_df_pathvars) >= 2:
            for i, r in gene_df_pathvars.iterrows():
                phase = r['Proband_GT']
                if (phase != '1/0') & (phase != '0/1'):  # if phased, only retain variants in trans with path variant
                    recessive_df.iloc[(recessive_df['Gene'] == gene) & (recessive_df['Proband_GT'] != phase), recessive_df.columns.get_loc('retain')] = 1
                    recessive_df.iloc[(recessive_df['#CHROM'] == r['#CHROM']) & (recessive_df['POS'] == r['POS']), recessive_df.columns.get_loc('retain')] = 1
                else:  # if unphased, retain all variants in gene
                    recessive_df.iloc[(recessive_df['Gene'] == gene), recessive_df.columns.get_loc('retain')] = 1
        # 2 loss of function variants, defined by loeuf score, in trans (if phased)
        gene_df_lofvars = gene_df[gene_df['Loeuf_score'] <= 0.4]
        if len(gene_df_lofvars) >= 2:
            for i, r in gene_df_lofvars.iterrows():
                phase = r['Proband_GT']
                if (phase != '1/0') & (phase != '0/1'):  # if phased, only retain variants in trans with path variant
                    recessive_df.iloc[(recessive_df['Gene'] == gene) & (
                                recessive_df['Proband_GT'] != phase), recessive_df.columns.get_loc('retain')] = 1
                    recessive_df.iloc[(recessive_df['CHROM'] == r['CHROM']) & (
                                recessive_df['POS'] == r['POS']), recessive_df.columns.get_loc('retain')] = 1
                else:  # if unphased, retain all variants in gene
                    recessive_df.iloc[(recessive_df['Gene'] == gene), recessive_df.columns.get_loc('retain')] = 1
        # if not lof, retain if vus or higher, or in HGMD
        elif len(gene_df_lofvars) < 1:
            gene_df_vus_hgmd = gene_df[(gene_df['Clinvar_sig'] >= 1) | (gene_df['hgmd_class'] >= 1) | (gene_df['Morbid_Gene'] == 1)]
            if len(gene_df_vus_hgmd) >= 1:
                for i, r in gene_df_vus_hgmd.iterrows():
                    phase = r['Proband_GT']
                    if (phase != '1/0') & (
                            phase != '0/1'):  # if phased, only retain variants in trans with path variant
                        recessive_df.iloc[(recessive_df['Gene'] == gene) & (
                                recessive_df['Proband_GT'] != phase), recessive_df.columns.get_loc('retain')] = 1
                        recessive_df.iloc[(recessive_df['CHROM'] == r['CHROM']) & (
                                recessive_df['POS'] == r['POS']), recessive_df.columns.get_loc('retain')] = 1
                    else:  # if unphased, retain all variants in gene
                        recessive_df.iloc[(recessive_df['Gene'] == gene), recessive_df.columns.get_loc('retain')] = 1
        # if 1 lof variant, check revel score or splice AI or constrained region for other variant. cadd and spliceAI for intronic variant
        gene_df_impactlofvars = gene_df[(gene_df['Highest_Impact_Order'] <= 7)]
        if len(gene_df_impactlofvars) == 1:
            for i, r in gene_df_impactlofvars.iterrows():
                phase = r['Proband_GT']
                if (phase != '1/0') & (
                        phase != '0/1'):  # if phased, only retain variants in trans with path variant
                    recessive_df.iloc[(recessive_df['Gene'] == gene) & (recessive_df['Proband_GT'] != phase) &
                                      ((recessive_df['Revel'] >= 0.5) | (recessive_df['SpliceAI_maxscore'] >= 0.5) |
                                       (recessive_df['CADD_PHRED'] >= 10)), recessive_df.columns.get_loc('retain')] = 1
                else:  # if unphased, retain all variants in gene
                    recessive_df.iloc[(recessive_df['Gene'] == gene) &
                                      ((recessive_df['Revel'] >= 0.5) | (recessive_df['SpliceAI_maxscore'] >= 0.5) |
                                       (recessive_df['CADD_PHRED'] >= 10)), recessive_df.columns.get_loc('retain')] = 1
    # utr-annotator for promoter variant

    # filter out genes in curated noisy genes list:
    # MUC (mucin) genes
    muc_pattern = r"MUC\d+"
    recessive_df = recessive_df[~recessive_df['Gene'].str.contains(muc_pattern, regex=True, na=False)]

    retain_df = recessive_df[recessive_df['retain'] == 1]
    retain_df = retain_df.sort_values(by=['GenCC_Submission', 'Gene', 'POS'])
    retain_df.to_csv(outdir + proband_ID + '.RecessiveVariants.tsv', sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rank variants for final output',
                                     formatter_class=argparse.RawTextHelpFormatter)
    group1 = parser.add_argument_group('Mandatory inputs')
    group1.add_argument('-InVCF', type=str, dest='in_vcf', required=True,
                        help='Input vcf file in unzipped format')
    group1.add_argument('-ped', type=str, dest='ped', required=True,
                        help='Input pedigree file')
    group1.add_argument('-proband', type=str, dest='proband_ID', required=True,
                        help='ID of proband')
    group1.add_argument('-OutDir', type=str, dest='outdir', required=True,
                        help='output dir for output files')
    group1.add_argument('-loeuf', type=str, dest='loeuf', required=True,
                        help='TSV file with precomputed loeuf scores')
    group1.add_argument('-morbidgenes', type=str, dest='morbidgenes', required=True,
                        help='Morbid genes list')
    group1.add_argument('-gencc', type=str, dest='gencc', required=True,
                        help='TSV file of GenCC submissions')
    group1.add_argument('-regconst', type=str, dest='regconst', required=True,
                        help='TSV file with regional missense constraint scores from gnomad')

    args = parser.parse_args()
    in_vcf = args.in_vcf
    ped = args.ped
    proband_ID = args.proband_ID
    outdir = args.outdir
    loeuf = args.loeuf
    morbidgenes = args.morbidgenes
    gencc = args.gencc
    regconst = args.regconst

    main(in_vcf, ped, proband_ID, outdir, loeuf, morbidgenes, gencc, regconst)
