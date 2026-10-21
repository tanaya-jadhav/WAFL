#!/usr/bin/env python

import sys
import pandas as pd
from peds import open_ped
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


def check_for_parents(sv_df, sv_headerline, father, mother):
    if (father is not None) and (father.id in sv_headerline) and (mother is not None) and (mother.id in sv_headerline):
        sampletype = 'trio'
    elif (father is not None) and (father.id in sv_headerline) and (mother.id not in sv_headerline):
        sampletype = 'duo_dad'
    elif (mother is not None) and (mother.id in sv_headerline) and (father.id not in sv_headerline):
        sampletype = 'duo_mom'
    else:
        sampletype = 'singleton'
    pd.set_option('display.max_columns', 500)
    # pull genotypes out as columns for analysis
    if sampletype == 'trio':
        sv_df['Dad_GT'] = sv_df[father.id].str.split(':').str[0]
        sv_df['Mom_GT'] = sv_df[mother.id].str.split(':').str[0]
    elif sampletype == 'duo_dad':
        sv_df['Parent_GT'] = sv_df[father.id].str.split(':').str[0]
    elif sampletype == 'duo_mom':
        sv_df['Parent_GT'] = sv_df[mother.id].str.split(':').str[0]
    # remove variants that are homozygous in parents
    if sampletype == 'trio':
        sv_df = sv_df[(sv_df['Dad_GT'] != '1/1') | (sv_df['Dad_GT'] != '1|1') |
                      (sv_df['Mom_GT'] != '1/1') | (sv_df['Mom_GT'] != '1|1')]
    elif (sampletype == 'duo_dad') | (sampletype == 'duo_mom'):
        sv_df = sv_df[(sv_df['Parent_GT'] != '1/1') | (sv_df['Parent_GT'] != '1|1')]
    return sv_df


def main():
    parser = argparse.ArgumentParser(description='Assign extended comphets', formatter_class=argparse.RawTextHelpFormatter)
    group1 = parser.add_argument_group('Mandatory inputs')
    group1.add_argument('-snv', type=str, dest='snv_file', required=True,
                        help='Input SNV vcf file in unzipped format')
    group1.add_argument('-ins', type=str, dest='ins_file', required=True,
                        help='Input insertions tsv file in unzipped format')
    group1.add_argument('-inv', type=str, dest='inv_file', required=True,
                        help='Input inversions tsv file in unzipped format')
    group1.add_argument('-cnv', type=str, dest='cnv_file', required=True,
                        help='Input CNV tsv file in unzipped format')
    group1.add_argument('-o', type=str, dest='out_file', required=True,
                        help='Output vcf file in unzipped format')
    group1.add_argument('-proband', type=str, dest='proband_ID', required=True,
                        help='ID of proband')
    group1.add_argument('-ped', type=str, dest='ped', required=True,
                        help='Input pedigree file')


    args = parser.parse_args()

    snv_file = args.snv_file
    ins_file = args.ins_file
    inv_file = args.inv_file
    cnv_file = args.cnv_file
    out_file = args.out_file
    proband_ID = args.proband_ID
    ped = args.ped
    # loeuf = args.loeuf
    # morbidgenes = args.morbidgenes
    # gencc = args.gencc
    #
    # loeuf_df = pd.read_csv(loeuf, sep='\t', header=0)
    # morbidgenes_df = pd.read_csv(morbidgenes, sep='\t', header=None, names=['gene', 'Morbid_Gene'])
    # gencc_df = pd.read_csv(gencc, sep='\t', header=0)
    pd.set_option('display.max_columns', None)
    # parse ped file
    # proband, father, mother = read_pedfile(ped, proband_ID)

    # # add morbid gene annotation
    # sv_df = sv_df.merge(morbidgenes_df, left_on='Gene', right_on='gene', how='left')
    #
    # # add gencc_morbid gene annotation
    # sv_df = sv_df.merge(gencc_df, left_on='Gene', right_on='Gene', how='left')
    # pd.set_option('display.max_columns', None)
    #
    # # determine if parent samples are present in sv vcf
    # sv_df = check_for_parents(sv_df, sv_headerline, father, mother)

    # parse snv vcf file
    with open(snv_file, "r") as i:
        snv_header = i.readlines()
    snv_header = [line for line in snv_header if line.startswith('#')]
    snv_headerline = [line for line in snv_header if line.startswith('#CHROM')][0]
    snv_headerline = snv_headerline.strip('\n').split('\t')

    snv_df = pd.read_csv(snv_file, sep='\t', comment="#", header=None, names=snv_headerline)
    # snv_df = snv_df[snv_df['FILTER'] == 'PASS'].reset_index()  # for variants with PASS only
    snv_df = snv_df.reset_index()
    snv_df['Proband_GT'] = snv_df[proband_ID].str.split(':').str[0]
    snv_df['CSQ'] = snv_df['INFO'].str.split(';CSQ=').str[1].str.split(',').str[0].str.split(';').str[0]
    snv_df['Gene'] = snv_df['CSQ'].str.split('|').str[3]
    snv_df['Consequence'] = snv_df['CSQ'].str.split('|').str[1]
    snv_df['Impact'] = snv_df['CSQ'].str.split('|').str[2]
    snv_df['HGVSc'] = snv_df['CSQ'].str.split('|').str[10]
    snv_df['HGVSp'] = snv_df['CSQ'].str.split('|').str[11]
    snv_df['gnomad_popmax_af'] = snv_df['INFO'].str.split('gnomad_popmax_af=').str[1].str.split(';').str[0]
    snv_df['hprc_FRQ'] = snv_df['INFO'].str.split('hprc_FRQ=').str[1].str.split(';').str[0]
    snv_df['rimgc_FRQ'] = snv_df['INFO'].str.split('rimgc_FRQ=').str[1].str.split(';').str[0]
    snv_df['Hgmd_class'] = snv_df['INFO'].str.split('HGVS_CLASS=').str[1].str.split(';').str[0]
    snv_df['Hgmd_rankscore'] = snv_df['INFO'].str.split('HGMD_RANKSCORE=').str[1].str.split(';').str[0]
    snv_df['CLNSIG'] = snv_df['INFO'].str.split('CLNSIG=').str[1].str.split(';').str[0]
    snv_df['CLNDN'] = snv_df['INFO'].str.split('CLNDN=').str[1].str.split(';').str[0]
    snv_df = snv_df[(snv_df['Proband_GT'] != '0/0') & (snv_df['Proband_GT'] != '0|0')]
    # print(snv_df.head(5))

    # filter out genes in curated noisy genes list:
    # MUC (mucin) genes
    muc_pattern = r"MUC\d+"
    snv_df = snv_df[~snv_df['Gene'].str.contains(muc_pattern, regex=True, na=False)]

    snv_genelist = set(snv_df['Gene'])

    # parse and filter insertions file
    ins_df = pd.read_csv(ins_file, sep='\t', header=0)
    ins_df['Proband_GT'] = ins_df[proband_ID].str.split(':').str[0]
    ins_genelist = set(ins_df['Gene_name'])
    common_genes = ins_genelist & snv_genelist
    inslist = []
    if len(common_genes) > 1:
        for cg in common_genes:
            cg_ins = ins_df[ins_df['Gene_name'] == cg]
            cg_snv = snv_df[snv_df['Gene'] == cg]
            for ins_i, ins_r in cg_ins.iterrows():
                outdict = {'chrom': ins_r['SV_chrom'], 'start': ins_r['SV_start'], 'Gene': ins_r['Gene_name'],
                           'TYPE': 'INS', 'REF': ins_r['REF'], 'ALT': ins_r['ALT'],
                           'Zygosity': ins_r['Proband_GT'],
                           'LEN': ins_r['LEN'], 'GenCC': ins_r['GenCC_disease'],
                           'OMIM_phenotype': ins_r['OMIM_phenotype'],
                           'annotsv_location_tag': ins_r['annotsv_location_tag'],
                           'location_boundary_detail': ins_r['location_boundary_detail'],
                           'cds_utr_context': ins_r['cds_utr_context'], 'nearest_junction': ins_r['nearest_junction'],
                           'nearest_junction_dist': ins_r['nearest_junction_dist'], 'Gnomad_AF': ins_r['Gnomad_AF'],
                           'Best_gnomAD_OFP': ins_r['Best_gnomAD_OFP'], 'HPRC_AF': ins_r['HPRC_AF'],
                           'Best_HPRC_OFP': ins_r['Best_HPRC_OFP'], 'CoLoRSdb_AF': ins_r['CoLoRSdb_AF'],
                           'Best_CoLoRSdb_OFP': ins_r['Best_CoLoRSdb_OFP'], 'CdLS_AF': ins_r['CdLS_AF'],
                           'Best_CdLS_OFP': ins_r['Best_CdLS_OFP'], 'ID': ins_r['ID']}
                inslist.append(outdict)
            for snv_i, snv_r in cg_snv.iterrows():
                outdict = {'chrom': snv_r['#CHROM'], 'start': snv_r['POS'], 'REF': snv_r['REF'],
                           'ALT': snv_r['ALT'], 'Gene': snv_r['Gene'], 'TYPE': 'SNV_Indel',
                           'Zygosity': snv_r['Proband_GT'], 'Gnomad_AF': snv_r['gnomad_popmax_af'],
                           'Consequence': snv_r['Consequence'], 'Impact': snv_r['Impact'], 'HGVSc': snv_r['HGVSc'],
                           'HGVSp': snv_r['HGVSp'], snv_r['CLNSIG']: snv_r['CLNSIG'], snv_r['CLNDN']: snv_r['CLNDN'],
                           'INFO': snv_r['INFO']}
                inslist.append(outdict)
    ins_outdf = pd.DataFrame.from_dict(inslist)
    # print(ins_df.head(5))

    # parse and filter inversions file
    inv_df = pd.read_csv(inv_file, sep='\t', header=0)
    # inv_df = inv_df[['SV_chrom', 'SV_start', 'SV_end', 'SV_length', proband_ID, 'Gene_name', 'Frameshift']]
    inv_df = inv_df[inv_df['annotsv_location_tag'] != 'whole_gene']
    inv_df['Proband_GT'] = inv_df[proband_ID].str.split(':').str[0]
    inv_genelist = set(inv_df['Gene_name'])
    common_genes = inv_genelist & snv_genelist
    invlist = []
    if len(common_genes) > 1:
        for cg in common_genes:
            cg_inv = inv_df[inv_df['Gene_name'] == cg]
            cg_snv = snv_df[snv_df['Gene'] == cg]
            for inv_i, inv_r in cg_inv.iterrows():
                outdict = {'chrom': inv_r['SV_chrom'], 'start': inv_r['SV_start'], 'end': inv_r['SV_end'],
                           'REF': inv_r['REF'], 'ALT': inv_r['ALT'], 'Gene': inv_r['Gene_name'],
                           'TYPE': inv_r['SVTYPE'],
                           'Zygosity': inv_r['Proband_GT'],
                           'LEN': inv_r['LEN'], 'GenCC': inv_r['GenCC_disease'],
                           'OMIM_phenotype': inv_r['OMIM_phenotype'],
                           'annotsv_location_tag': ins_r['annotsv_location_tag'],
                           'location_boundary_detail': ins_r['location_boundary_detail'],
                           'cds_utr_context': ins_r['cds_utr_context'], 'nearest_junction': ins_r['nearest_junction'],
                           'nearest_junction_dist': ins_r['nearest_junction_dist'], 'Gnomad_AF': inv_r['Gnomad_AF'],
                           'Best_gnomAD_OFP': inv_r['Best_gnomAD_OFP'], 'HPRC_AF': inv_r['HPRC_AF'],
                           'Best_HPRC_OFP': inv_r['Best_HPRC_OFP'], 'CoLoRSdb_AF': inv_r['CoLoRSdb_AF'],
                           'Best_CoLoRSdb_OFP': inv_r['Best_CoLoRSdb_OFP'], 'CdLS_AF': inv_r['CdLS_AF'],
                           'Best_CdLS_OFP': inv_r['Best_CdLS_OFP'], 'INFO': inv_r['INFO'], 'ID': inv_r['INV_ID']}
                invlist.append(outdict)
            for snv_i, snv_r in cg_snv.iterrows():
                outdict = {'chrom': snv_r['#CHROM'], 'start': snv_r['POS'], 'REF': snv_r['REF'],
                           'ALT': snv_r['ALT'], 'Gene': snv_r['Gene'], 'TYPE': 'SNV_Indel',
                           'Zygosity': snv_r['Proband_GT'], 'Gnomad_AF': snv_r['gnomad_popmax_af'],
                           'Consequence': snv_r['Consequence'], 'Impact': snv_r['Impact'], 'HGVSc': snv_r['HGVSc'],
                           'HGVSp': snv_r['HGVSp'], snv_r['CLNSIG']: snv_r['CLNSIG'], snv_r['CLNDN']: snv_r['CLNDN'],
                           'INFO': snv_r['INFO']}
                invlist.append(outdict)
    inv_outdf = pd.DataFrame.from_dict(invlist)

    # parse and filter cnv file
    cnv_df = pd.read_csv(cnv_file, sep='\t', header=0)
    # cnv_df = cnv_df[
    #     ['chrom', 'start', 'end', 'LEN', 'SVTYPE', proband_ID, 'Gene_name', 'Frameshift', 'annotsv_location_tag']]
    cnv_df['Proband_GT'] = cnv_df[proband_ID].str.split(':').str[0]
    cnv_df = cnv_df[~((cnv_df['SVTYPE'] == 'DUP') & (cnv_df['annotsv_location_tag'] == 'whole_gene'))]
    cnv_genelist = set(cnv_df['Gene_name'])
    common_genes = cnv_genelist & snv_genelist
    # print(common_genes)
    outlist = []
    if len(common_genes) > 1:
        for cg in common_genes:
            cg_cnv = cnv_df[cnv_df['Gene_name'] == cg]
            cg_snv = snv_df[snv_df['Gene'] == cg]
            for cnv_i, cnv_r in cg_cnv.iterrows():
                outdict = {'chrom': cnv_r['chrom'], 'start': cnv_r['start'], 'end': cnv_r['end'],
                           'REF': cnv_r['REF'], 'ALT': cnv_r['ALT'], 'Gene': cnv_r['Gene_name'], 'TYPE': cnv_r['SVTYPE'],
                           'Zygosity': cnv_r['Proband_GT'],
                            'LEN': cnv_r['LEN'], 'GenCC': cnv_r['GenCC_disease'],
                           'OMIM_phenotype': cnv_r['OMIM_phenotype'],
                           'annotsv_location_tag': ins_r['annotsv_location_tag'],
                           'location_boundary_detail': ins_r['location_boundary_detail'],
                           'cds_utr_context': ins_r['cds_utr_context'], 'nearest_junction': ins_r['nearest_junction'],
                           'nearest_junction_dist': ins_r['nearest_junction_dist'], 'Gnomad_AF': cnv_r['Gnomad_AF'],
                           'Best_gnomAD_OFP': cnv_r['Best_gnomAD_OFP'],	'HPRC_AF': cnv_r['HPRC_AF'],
                           'Best_HPRC_OFP': cnv_r['Best_HPRC_OFP'], 'CoLoRSdb_AF': cnv_r['CoLoRSdb_AF'],
                           'Best_CoLoRSdb_OFP': cnv_r['Best_CoLoRSdb_OFP'], 'CdLS_AF': cnv_r['CdLS_AF'],
                           'Best_CdLS_OFP': cnv_r['Best_CdLS_OFP'], 'INFO': cnv_r['INFO'], 'ID': cnv_r['ID']}
                outlist.append(outdict)
            for snv_i, snv_r in cg_snv.iterrows():
                outdict = {'chrom': snv_r['#CHROM'], 'start': snv_r['POS'], 'REF': snv_r['REF'],
                           'ALT': snv_r['ALT'], 'Gene': snv_r['Gene'], 'TYPE': 'SNV_Indel',
                           'Zygosity': snv_r['Proband_GT'], 'Gnomad_AF': snv_r['gnomad_popmax_af'],
                           'Consequence': snv_r['Consequence'], 'Impact': snv_r['Impact'], 'HGVSc': snv_r['HGVSc'],
                           'HGVSp': snv_r['HGVSp'], snv_r['CLNSIG']: snv_r['CLNSIG'], snv_r['CLNDN']: snv_r['CLNDN'],
                           'INFO': snv_r['INFO']}
                outlist.append(outdict)
    cnv_outdf = pd.DataFrame.from_dict(outlist)

    dfs_to_concat = [x for x in [ins_outdf, cnv_outdf, inv_outdf] if len(x) > 0]
    combined_df = pd.concat(dfs_to_concat, ignore_index=True)
    combined_df = combined_df[['chrom', 'start', 'end', 'REF', 'ALT', 'Gene', 'TYPE', 'Zygosity', 'LEN', 'GenCC',
                           'OMIM_phenotype', 'annotsv_location_tag', 'location_boundary_detail', 'cds_utr_context',
                               'nearest_junction', 'nearest_junction_dist', 'Gnomad_AF', 'Best_gnomAD_OFP', 'HPRC_AF',
                           'Best_HPRC_OFP', 'CoLoRSdb_AF', 'Best_CoLoRSdb_OFP', 'CdLS_AF', 'Best_CdLS_OFP', 'INFO',
                               'ID']]
    combined_df.to_csv(out_file, sep='\t', index=False)


if __name__ == '__main__':
    main()
