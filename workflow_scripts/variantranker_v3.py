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
        return np.nan
    elif max(colAG, colAL, colDG, colDL) == colAG:
        return 'Delta score for acceptor gain'
    elif max(colAG, colAL, colDG, colDL) == colAL:
        return 'Delta score for acceptor loss'
    elif max(colAG, colAL, colDG, colDL) == colDG:
        return 'Delta score for donor gain'
    elif max(colAG, colAL, colDG, colDL) == colDL:
        return 'Delta score for donor loss'
    else:
        return np.nan


def getregionalconsscore(regionalconstraint_df, varchrom, varpos):
    varchrom = varchrom.split('chr')[1]
    varpos = int(varpos)  # POS is read as str; comparing str to the float position columns raises TypeError
    try:
        vardf = regionalconstraint_df[(regionalconstraint_df['hg38_chr'] == varchrom) &
                                      (regionalconstraint_df['hg38_startpos'] <= varpos) &
                                      (regionalconstraint_df['hg38_endpos'] >= varpos)]
        score = list(vardf['oe'])[0]
        # print(varchrom, varpos)
        # print(vardf)
        return score
    except IndexError:  # variant not in any regional constraint interval
        return 0


# dbNSFP score/prediction pairs added to CSQ by the VEP dbNSFP plugin
DBNSFP_FIELDS = [('AlphaMissense_score', 'AlphaMissense_pred'), ('PrimateAI_score', 'PrimateAI_pred'),
                 ('BayesDel_addAF_score', 'BayesDel_addAF_pred'), ('BayesDel_noAF_score', 'BayesDel_noAF_pred')]
DBNSFP_COLS = [col for pair in DBNSFP_FIELDS for col in pair]


def maxdbnsfpscore(scores, preds):
    # dbNSFP reports one value per Ensembl transcript joined by '&' ('.' = no score for that transcript);
    # keep the highest score and the prediction from the same transcript
    scores = str(scores).split('&')
    preds = str(preds).split('&')
    best_score, best_pred = np.nan, np.nan
    for i, s in enumerate(scores):
        try:
            s = float(s)
        except ValueError:
            continue
        if np.isnan(best_score) or s > best_score:
            best_score = s
            best_pred = preds[i] if i < len(preds) and preds[i] not in ('', '.') else np.nan
    return best_score, best_pred


# D_score weights, fitted in clinvar_testing/clinvar_20251221_snvindels (see computedscore)
DSCORE_WEIGHTS = {
    'BayesDel_addAF': 7.13,      # BayesDel_addAF_score; 0 when not scored
    'AlphaMissense': 3.80,       # AlphaMissense_score; 0.34 (benign/ambiguous cutoff) for missense without a score
    'Revel': 1.95,
    'CADD_PHRED': 0.05,
    'SpliceAI_max': 5.79,        # SpliceAI_maxscore
    'Impact_significance': 4.14,  # HIGH=3, MODERATE=2, LOW/MODIFIER=1
    'gnomad_rarity': 2.41,       # -log10(gnomad_popmax_af); absent from gnomAD counted as AF 1e-6
    'no_missense_score': 4.59,   # 1 when there is no BayesDel score (non-missense variant)
}
DSCORE_INTERCEPT = -26.1
# not estimable from ClinVar data (no genotypes); ln(18.7), the odds for ACMG 'strong' evidence (Tavtigian et al. 2018)
DSCORE_DENOVO_WEIGHT = 2.93


def computedscore(df):
    # gnomad_popmax_af appears twice in the feature table; use the first copy
    af = df['gnomad_popmax_af']
    if isinstance(af, pd.DataFrame):
        af = af.iloc[:, 0]
    af = pd.to_numeric(af, errors='coerce')
    bayesdel = pd.to_numeric(df['BayesDel_addAF_score'], errors='coerce')
    alphamissense = pd.to_numeric(df['AlphaMissense_score'], errors='coerce')
    no_missense_score = bayesdel.isna()
    terms = {
        'BayesDel_addAF': bayesdel.fillna(0),
        'AlphaMissense': alphamissense.where(alphamissense.notna() | no_missense_score, 0.34).fillna(0),
        'Revel': pd.to_numeric(df['Revel'], errors='coerce').fillna(0),
        'CADD_PHRED': pd.to_numeric(df['CADD_PHRED'], errors='coerce').fillna(0),
        'SpliceAI_max': pd.to_numeric(df['SpliceAI_maxscore'], errors='coerce').fillna(0),
        'Impact_significance': pd.to_numeric(df['Impact_significance'], errors='coerce').fillna(1),
        'gnomad_rarity': -np.log10(af.where(af > 0, 1e-6).clip(lower=1e-6)),
        'no_missense_score': no_missense_score.astype(int),
    }
    score = DSCORE_INTERCEPT + DSCORE_DENOVO_WEIGHT * df['Denovo']
    for name, weight in DSCORE_WEIGHTS.items():
        score = score + weight * terms[name]
    return score


# ClinVar tiers that are ranked above every other variant: Pathogenic > Likely pathogenic > VUS.
# A variant takes the best tier among the '/'-, '|'- or ','-separated terms of its CLNSIG
# (e.g. Pathogenic/Likely_pathogenic -> 3, Uncertain_significance/Uncertain_risk_allele -> 1);
# Conflicting_classifications_of_pathogenicity, benign, risk-allele-only and missing CLNSIG stay in tier 0.
CLINVAR_TIERS = {'Pathogenic': 3, 'Likely_pathogenic': 2, 'Uncertain_significance': 1}
# added per tier; the model D_score is a log-odds well inside (-50, 50), so tiers never overlap
CLINVAR_TIER_OFFSET = 100


def clinvartier(clnsig):
    terms = str(clnsig).replace('|', '/').replace(',', '/').split('/')
    return max([CLINVAR_TIERS.get(t, 0) for t in terms])


def getcsqindices(vcf_header):
    # map each VEP CSQ field name to its position, using the ##INFO=<ID=CSQ ... Format: a|b|c"> header line,
    # so the field order can change between VEP runs without silently shifting features
    csq_lines = [line for line in vcf_header if line.startswith('##INFO=<ID=CSQ')]
    if not csq_lines:
        print('No ##INFO=<ID=CSQ header line found in input VCF\nExiting...')
        sys.exit(1)
    csq_format = csq_lines[0].split('Format: ')[1].split('"')[0].strip()
    csq_idx = {name: i for i, name in enumerate(csq_format.split('|'))}
    required = ['Consequence', 'IMPACT', 'SYMBOL', 'EXON', 'HGVSc', 'HGVSp', 'SpliceAI_pred_DS_AG',
                'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL', 'REVEL', 'LoFtool',
                'CADD_PHRED', 'CADD_RAW', 'CLINVAR_Clinvar_Pathogenic', 'CLINVAR_Clinvar_Benign',
                'CLINVAR_Clinvar_Uncertain.significance']
    missing = [name for name in required if name not in csq_idx]
    if missing:
        print('CSQ header is missing required fields: ' + ', '.join(missing) + '\nExiting...')
        sys.exit(1)
    # dbNSFP fields are optional so VCFs annotated before the dbNSFP plugin was added can still be ranked
    missing_dbnsfp = [name for name in DBNSFP_COLS if name not in csq_idx]
    if missing_dbnsfp:
        print('CSQ header is missing dbNSFP fields (will be left empty): ' + ', '.join(missing_dbnsfp))
    return csq_idx


def createfeaturedf(proband_ID, in_df, sampletype, outdir, loeuf_df, morbidgenes_df, gencc_df, regionalconstraint_df, csq_idx):
    # pull out features for scoring
    in_df = in_df.rename(columns={"#CHROM": "CHROM"})
    in_df['CSQ'] = in_df['INFO'].str.split('CSQ=').str[1].str.split(',').str[0].str.split(';').str[0]
    in_df['Gene'] = in_df['CSQ'].str.split('|').str[csq_idx['SYMBOL']]
    in_df['Exon'] = in_df['CSQ'].str.split('|').str[csq_idx['EXON']]
    in_df['Consequence'] = in_df['CSQ'].str.split('|').str[csq_idx['Consequence']]
    in_df['Impact'] = in_df['CSQ'].str.split('|').str[csq_idx['IMPACT']]
    in_df['Highest_Impact_Order'] = in_df['INFO'].str.split('highest_impact_order=').str[1].str.split(';').str[0]
    in_df['Highest_Impact_Order'] = pd.to_numeric(in_df['Highest_Impact_Order'], errors='coerce')
    in_df['HGVSc'] = in_df['CSQ'].str.split('|').str[csq_idx['HGVSc']]
    in_df['HGVSp'] = in_df['CSQ'].str.split('|').str[csq_idx['HGVSp']]
    in_df['gnomad_popmax_af'] = in_df['INFO'].str.split('gnomad_popmax_af=').str[1].str.split(';').str[0]
    in_df['hprc_FRQ'] = in_df['INFO'].str.split('hprc_FRQ=').str[1].str.split(';').str[0]
    in_df['rimgc_FRQ'] = in_df['INFO'].str.split('rimgc_FRQ=').str[1].str.split(';').str[0]
    in_df['Hgmd_class'] = in_df['INFO'].str.split('HGVS_CLASS=').str[1].str.split(';').str[0]
    in_df['Hgmd_rankscore'] = in_df['INFO'].str.split('HGMD_RANKSCORE=').str[1].str.split(';').str[0]
    in_df['CLNSIG'] = in_df['INFO'].str.split('CLNSIG=').str[1].str.split(';').str[0]
    in_df['CLNDN'] = in_df['INFO'].str.split('CLNDN=').str[1].str.split(';').str[0]
    in_df['CompHet'] = np.where(in_df['INFO'].str.contains('Compound_Heterozygous'), 1, 0)
    in_df['Homozygous'] = np.where((in_df['Proband_GT'] == '1|1') | (in_df['Proband_GT'] == 1 / 1), 1, 0)
    in_df['SpliceAI_pred_DS_AG'] = in_df['CSQ'].str.split('|').str[csq_idx['SpliceAI_pred_DS_AG']]
    in_df['SpliceAI_pred_DS_AG'] = pd.to_numeric(in_df['SpliceAI_pred_DS_AG'], errors='coerce')
    in_df['SpliceAI_pred_DS_AL'] = in_df['CSQ'].str.split('|').str[csq_idx['SpliceAI_pred_DS_AL']]
    in_df['SpliceAI_pred_DS_AL'] = pd.to_numeric(in_df['SpliceAI_pred_DS_AL'], errors='coerce')
    in_df['SpliceAI_pred_DS_DG'] = in_df['CSQ'].str.split('|').str[csq_idx['SpliceAI_pred_DS_DG']]
    in_df['SpliceAI_pred_DS_DG'] = pd.to_numeric(in_df['SpliceAI_pred_DS_DG'], errors='coerce')
    in_df['SpliceAI_pred_DS_DL'] = in_df['CSQ'].str.split('|').str[csq_idx['SpliceAI_pred_DS_DL']]
    in_df['SpliceAI_pred_DS_DL'] = pd.to_numeric(in_df['SpliceAI_pred_DS_DL'], errors='coerce')
    in_df['Revel'] = in_df['CSQ'].str.split('|').str[csq_idx['REVEL']]
    in_df['LoFtool'] = in_df['CSQ'].str.split('|').str[csq_idx['LoFtool']]
    in_df['CADD_PHRED'] = in_df['CSQ'].str.split('|').str[csq_idx['CADD_PHRED']]
    in_df['CADD_PHRED'] = pd.to_numeric(in_df['CADD_PHRED'], errors='coerce').fillna(0)
    in_df['CADD_RAW'] = in_df['CSQ'].str.split('|').str[csq_idx['CADD_RAW']]
    in_df['Clinvar_path'] = in_df['CSQ'].str.split('|').str[csq_idx['CLINVAR_Clinvar_Pathogenic']]
    in_df['Clinvar_path'] = pd.to_numeric(in_df['Clinvar_path'], errors='coerce')
    in_df['Clinvar_benign'] = in_df['CSQ'].str.split('|').str[csq_idx['CLINVAR_Clinvar_Benign']]
    in_df['Clinvar_benign'] = pd.to_numeric(in_df['Clinvar_benign'], errors='coerce')
    in_df['Clinvar_vus'] = in_df['CSQ'].str.split('|').str[csq_idx['CLINVAR_Clinvar_Uncertain.significance']]
    in_df['Clinvar_vus'] = pd.to_numeric(in_df['Clinvar_vus'], errors='coerce')
    in_df['SpliceAI_maxscore'] = in_df[
        ['SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL']].max(axis=1)
    # v = np.vectorize(maxspliceaiscore)
    # in_df['SpliceAI_maxtype'] = v(in_df.SpliceAI_pred_DS_AG, in_df.SpliceAI_pred_DS_AL, in_df.SpliceAI_pred_DS_DG,
    #                               in_df.SpliceAI_pred_DS_DL)
    in_df['SpliceAI_maxtype'] = in_df.apply(lambda x: maxspliceaiscore(x.SpliceAI_pred_DS_AG, x.SpliceAI_pred_DS_AL,
                                                                       x.SpliceAI_pred_DS_DG, x.SpliceAI_pred_DS_DL), axis=1)
    # dbNSFP scores (AlphaMissense, PrimateAI, BayesDel)
    for score_col, pred_col in DBNSFP_FIELDS:
        if score_col in csq_idx and pred_col in csq_idx:
            best = [maxdbnsfpscore(s, p) for s, p in zip(in_df['CSQ'].str.split('|').str[csq_idx[score_col]],
                                                         in_df['CSQ'].str.split('|').str[csq_idx[pred_col]])]
            in_df[score_col] = [b[0] for b in best]
            in_df[pred_col] = [b[1] for b in best]
        else:
            in_df[score_col] = np.nan
            in_df[pred_col] = np.nan
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
                            'Clinvar_path', 'Clinvar_benign', 'Clinvar_vus', 'gnomad_popmax_af', 'hprc_FRQ', 'rimgc_FRQ'] + DBNSFP_COLS]
    elif (sampletype == 'duo_dad') | (sampletype == 'duo_mom'):
        feature_df = in_df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'Proband_GT', 'Gene', 'Consequence', 'Impact',
                            'Highest_Impact_Order', 'HGVSc',
                            'HGVSp', 'gnomad_popmax_af', 'Hgmd_class', 'Hgmd_rankscore', 'CLNSIG', 'CLNDN', 'CompHet',
                            'Homozygous', 'Denovo', 'Parent_GT', 'Revel', 'LoFtool', 'Loeuf_score', 'regional_missense_constraint', 'SpliceAI_maxscore',
                            'SpliceAI_maxtype', 'CADD_PHRED', 'CADD_RAW', 'Morbid_Gene', 'GenCC_Submission', 'Clinvar_path',
                            'Clinvar_benign', 'Clinvar_vus', 'gnomad_popmax_af', 'hprc_FRQ', 'rimgc_FRQ'] + DBNSFP_COLS]
    else:
        feature_df = in_df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'Proband_GT', 'Gene', 'Consequence', 'Impact',
                            'Highest_Impact_Order', 'HGVSc',
                            'HGVSp', 'gnomad_popmax_af', 'Hgmd_class', 'Hgmd_rankscore', 'CLNSIG', 'CLNDN', 'CompHet',
                            'Homozygous', 'Denovo', 'Revel', 'LoFtool', 'Loeuf_score', 'regional_missense_constraint', 'SpliceAI_maxscore',
                            'SpliceAI_maxtype', 'CADD_PHRED', 'CADD_RAW', 'Morbid_Gene', 'GenCC_Submission', 'Clinvar_path',
                            'Clinvar_benign', 'Clinvar_vus', 'gnomad_popmax_af', 'hprc_FRQ', 'rimgc_FRQ'] + DBNSFP_COLS]
    feature_df['LoFtool'] = feature_df.LoFtool.replace('', 0, regex=True).astype(float)
    feature_df['SpliceAI_maxscore'] = feature_df.SpliceAI_maxscore.replace('', 0, regex=True).astype(float)
    feature_df['Revel'] = feature_df.Revel.replace('', 0, regex=True).astype(float)
    feature_df[['SpliceAI_maxscore', 'Morbid_Gene', 'Clinvar_path', 'Clinvar_benign', 'Clinvar_vus']] = feature_df[['SpliceAI_maxscore', 'Morbid_Gene', 'Clinvar_path', 'Clinvar_benign', 'Clinvar_vus']].fillna(value=0).astype(float)
    feature_df['Morbid_Gene'] = feature_df['Morbid_Gene'] + 1

    clinvar_conds = [feature_df['CLNSIG'].isin(['Pathogenic', 'Pathogenic/Likely_pathogenic', 'Pathogenic/Likely_pathogenic|risk_factor']),
             feature_df['CLNSIG'].isin(['Likely_pathogenic']),
             feature_df['CLNSIG'].isin(['Uncertain_significance', 'Conflicting_classifications_of_pathogenicity'])]
    clinvar_values = [3, 2, 1]
    feature_df['Clinvar_sig'] = np.select(clinvar_conds, clinvar_values, default=0)

    hgmd_conds = [feature_df['Hgmd_class'].isin(['DM']),
                     feature_df['Hgmd_class'].isin(['DM?'])]
    hgmd_values = [3, 2]
    feature_df['hgmd_class'] = np.select(hgmd_conds, hgmd_values, default=1)

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
    csq_idx = getcsqindices(vcf_header)

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
    # in_df = pd.read_csv(in_vcf, sep='\t', comment='#', header=None, names=headerline)
    # skip ## meta lines instead of using comment="#", which truncates lines with '#' inside INFO (e.g. ClinVar CLNVI)
    n_meta = sum(1 for line in vcf_header if line.startswith('##'))
    in_df = pd.read_csv(in_vcf, sep='\t', skiprows=n_meta, header=0, dtype=str)

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
    feature_df = createfeaturedf(proband_ID, in_df, sampletype, outdir, loeuf_df, morbidgenes_df, gencc_df, regionalconstraint_df, csq_idx)
    # print(feature_df.shape, feature_df.head(5))

    # separate by dominant and recessive
    recessive_df = feature_df[(feature_df['CompHet'] == 1) | (feature_df['Homozygous'] == 1)]
    dominant_df = feature_df[(feature_df['Homozygous'] == 0)]

    # remove variants that are homozygous in parents ~ behavior will be changed in future versions. needs an updated logic
    if sampletype == 'trio':
        dominant_df = dominant_df[(dominant_df['Dad_GT'] != '1/1') | (dominant_df['Dad_GT'] != '1|1') |
                                    (dominant_df['Mom_GT'] != '1/1') | (dominant_df['Mom_GT'] != '1|1')]
        recessive_df = recessive_df[(recessive_df['Dad_GT'] != '1/1') | (recessive_df['Dad_GT'] != '1|1') |
                                    (recessive_df['Mom_GT'] != '1/1') | (recessive_df['Mom_GT'] != '1|1')]
        # testing exception where 1 unaffected parent is homozygous, other is heterozygous, proband is heterozygous and affected
        # recessive_df.loc[((recessive_df['Homozygous'] == 1) & ((recessive_df['Dad_GT'] != '1/1') |
        #                             (recessive_df['Dad_GT'] != '1|1') | (recessive_df['Mom_GT'] != '1/1') |
        #                             (recessive_df['Mom_GT'] != '1|1'))), 'homozygous_exception'] = 1
        # recessive_df = recessive_df[recessive_df['homozygous_exception'] != 1]
    elif (sampletype == 'duo_dad') | (sampletype == 'duo_mom'):
        dominant_df = dominant_df[(dominant_df['Parent_GT'] != '1/1') | (dominant_df['Parent_GT'] != '1|1')]
        recessive_df = recessive_df[(recessive_df['Parent_GT'] != '1/1') | (recessive_df['Parent_GT'] != '1|1')]

    # score splice variants same as missense variants for impact in dominant vars
    dominant_df.loc[(dominant_df['Highest_Impact_Order'] == 20) | (dominant_df['Highest_Impact_Order'] == 40)
        | (dominant_df['Highest_Impact_Order'] == 42), 'Highest_Impact_Order'] = 14
    # Scoring dominant variants
    # D_score is the log-odds that a variant is P/LP rather than B/LB (logistic regression on 26,119 ClinVar P/LP and
    # B/LB SNV/indels in RTAG AD genes, clinvar_20251221; 5-fold CV AUC 0.998), plus a fixed de novo term.
    # D_score_model is that variant-level score (unchanged from variantranker_v2). D_score adds CLINVAR_TIER_OFFSET per
    # ClinVar tier so ClinVar Pathogenic, then Likely pathogenic, then VUS variants rank above all other variants;
    # within each tier (including tier 0) the order is the D_score_model order.
    dominant_df['D_score_model'] = computedscore(dominant_df)
    dominant_df['ClinVar_tier'] = dominant_df['CLNSIG'].map(clinvartier)
    dominant_df['D_score'] = dominant_df['D_score_model'] + CLINVAR_TIER_OFFSET * dominant_df['ClinVar_tier']

    # Tag variants for prioritization
    # all missense and deletion variants with regional constraint score <0.6
    missensevars = ['missense_variant', 'frameshift_variant', 'inframe_deletion']
    dominant_df.loc[(dominant_df['Consequence'].isin(missensevars)) & (dominant_df['regional_missense_constraint'] < 0.6)
                    & (dominant_df['regional_missense_constraint'] > 0), 'Variant_in_constrained_region'] = 1
    recessive_df.loc[
        (recessive_df['Consequence'].isin(missensevars)) & (recessive_df['regional_missense_constraint'] < 0.6)
        & (recessive_df['regional_missense_constraint'] > 0), 'Variant_in_constrained_region'] = 1

    # write dominant df to file
    dominant_df = dominant_df.sort_values(by=['ClinVar_tier', 'D_score_model'], ascending=False)
    dominant_df.to_csv(outdir + proband_ID + '.DominantVariants.tsv', sep='\t', index=False)

    # scoring recessive variants
    genelist = set(recessive_df['Gene'])
    recessive_df['retain'] = ''
    # print(recessive_df.columns)
    for gene in genelist:
        gene_df = recessive_df[recessive_df['Gene'] == gene]
        # if at least 1 variant is path or likely path, retain this variant and all trans variants in this gene (if phased)
        gene_df_pathvars = gene_df[gene_df['Clinvar_sig'] >= 2]
        if len(gene_df_pathvars) >= 2:
            for i, r in gene_df_pathvars.iterrows():
                phase = r['Proband_GT']
                if (phase != '1/0') & (phase != '0/1'):  # if phased, only retain variants in trans with path variant
                    # print(recessive_df.head(5))
                    recessive_df.iloc[(recessive_df['Gene'] == gene) & (recessive_df['Proband_GT'] != phase), recessive_df.columns.get_loc('retain')] = 1
                    recessive_df.iloc[(recessive_df['CHROM'] == r['CHROM']) & (recessive_df['POS'] == r['POS']), recessive_df.columns.get_loc('retain')] = 1
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
            gene_df_vus_hgmd = gene_df[(gene_df['Clinvar_sig'] >= 1) | (gene_df['hgmd_class'] > 1) | (gene_df['Morbid_Gene'] == 2)]
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
        gene_df_impactlofvars = gene_df[(gene_df['Highest_Impact_Order'] <= 7) | (gene_df['Highest_Impact_Order'] == 20) |
                                (gene_df['Highest_Impact_Order'] == 40) | (gene_df['Highest_Impact_Order'] == 42)]
        if len(gene_df_impactlofvars) == 1:
            for i, r in gene_df_impactlofvars.iterrows():
                phase = r['Proband_GT']
                if (phase != '1/0') & (
                        phase != '0/1'):  # if phased, only retain variants in trans with path variant
                    recessive_df.iloc[(recessive_df['Gene'] == gene) & (recessive_df['Proband_GT'] != phase) &
                                      ((recessive_df['Revel'] >= 0.5) | (recessive_df['SpliceAI_maxscore'] >= 0.5) |
                                       (recessive_df['CADD_PHRED'] >= 10) | (recessive_df['Variant_in_constrained_region'] == 1)), recessive_df.columns.get_loc('retain')] = 1
                else:  # if unphased, retain all variants in gene
                    recessive_df.iloc[(recessive_df['Gene'] == gene) &
                                      ((recessive_df['Revel'] >= 0.5) | (recessive_df['SpliceAI_maxscore'] >= 0.5) |
                                       (recessive_df['CADD_PHRED'] >= 10) | (recessive_df['Variant_in_constrained_region'] == 1)), recessive_df.columns.get_loc('retain')] = 1
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
