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


def get_comphetvars(sv_df, snv_df, out_file):
    sv_genelist = set(sv_df['Gene'])
    # print(sv_genelist)
    snv_genelist = set(snv_df['Gene'])
    # print(snv_genelist)
    common_genes = sv_genelist & snv_genelist
    # print(common_genes)
    # print(len(common_genes))
    sv_df = sv_df[['#CHROM', 'POS', 'REF', 'ALT', 'Proband_GT', 'Gene', 'Consequence', 'Impact', 'Consequence',
                   'Impact', 'HGVSc', 'HGVSp', 'GenCC_Submission', 'Morbid_Gene']]
    snv_df = snv_df[['#CHROM', 'POS', 'REF', 'ALT', 'Proband_GT', 'Gene', 'Consequence', 'Impact', 'Consequence',
                     'Impact', 'HGVSc', 'HGVSp', 'GenCC_Submission', 'Morbid_Gene']]
    sv_df = sv_df[sv_df['Gene'].isin(common_genes)]
    snv_df = snv_df[snv_df['Gene'].isin(common_genes)]
    combineddf = pd.concat([sv_df, snv_df], axis=0, ignore_index=True)
    combineddf = combineddf.sort_values(by=['Gene', 'POS'])
    combineddf.to_csv(out_file, sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser(description='Assign extended comphets', formatter_class=argparse.RawTextHelpFormatter)
    group1 = parser.add_argument_group('Mandatory inputs')
    group1.add_argument('-sv', type=str, dest='sv_file', required=True,
                        help='Input SV vcf file in unzipped format')
    group1.add_argument('-snv', type=str, dest='snv_file', required=True,
                        help='Input SNV vcf file in unzipped format')
    group1.add_argument('-o', type=str, dest='out_file', required=True,
                        help='Output vcf file in unzipped format')
    group1.add_argument('-proband', type=str, dest='proband_ID', required=True,
                        help='ID of proband')
    group1.add_argument('-ped', type=str, dest='ped', required=True,
                        help='Input pedigree file')

    args = parser.parse_args()

    sv_file = args.sv_file
    snv_file = args.snv_file
    out_file = args.out_file
    proband_ID = args.proband_ID
    ped = args.ped

    loeuf_file = '/mnt/isilon/rajagopalan_lab/projects/snv_indel_workflow/workingdir/referencefiles/gnomad.v4.1.constraint_metrics.filtered.tsv'
    loeuf_df = pd.read_csv(loeuf_file, sep='\t', header=0)
    morbidgenes_file = '/mnt/isilon/rajagopalan_lab/projects/snv_indel_workflow/workingdir/referencefiles/morbidmap.genes.txt'
    morbidgenes_df = pd.read_csv(morbidgenes_file, sep='\t', header=None, names=['gene', 'Morbid_Gene'])
    gencc_file = '/mnt/isilon/rajagopalan_lab/projects/snv_indel_workflow/workingdir/referencefiles/gencc_submissions_processed.tsv'
    gencc_df = pd.read_csv(gencc_file, sep='\t', header=0)

    # parse ped file
    proband, father, mother = read_pedfile(ped, proband_ID)

    # parse sv vcf file
    with open(sv_file, "r") as i:
        sv_header = i.readlines()
    sv_header = [line for line in sv_header if line.startswith('#')]
    sv_headerline = [line for line in sv_header if line.startswith('#CHROM')][0]
    sv_headerline = sv_headerline.strip('\n').split('\t')

    sv_df = pd.read_csv(sv_file, sep='\t', comment="#", header=None, names=sv_headerline)
    sv_df = sv_df[sv_df['FILTER'] == 'PASS'].reset_index()  # for variants with PASS only
    sv_df['Proband_GT'] = sv_df[proband_ID].str.split(':').str[0]
    sv_df['CSQ'] = sv_df['INFO'].str.split(';CSQ=').str[1].str.split(',').str[0].str.split(';').str[0]
    sv_df['Gene'] = sv_df['CSQ'].str.split('|').str[3]
    sv_df['Consequence'] = sv_df['CSQ'].str.split('|').str[1]
    sv_df['Impact'] = sv_df['CSQ'].str.split('|').str[2]
    sv_df['HGVSc'] = sv_df['CSQ'].str.split('|').str[10]
    sv_df['HGVSp'] = sv_df['CSQ'].str.split('|').str[11]
    sv_df = sv_df[(sv_df['Proband_GT'] != '0/0') & (sv_df['Proband_GT'] != '0|0')]

    # add morbid gene annotation
    sv_df = sv_df.merge(morbidgenes_df, left_on='Gene', right_on='gene', how='left')

    # add gencc_morbid gene annotation
    sv_df = sv_df.merge(gencc_df, left_on='Gene', right_on='Gene', how='left')
    pd.set_option('display.max_columns', None)

    # determine if parent samples are present in sv vcf
    sv_df = check_for_parents(sv_df, sv_headerline, father, mother)

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
    # print(snv_df[.head(5)])


    # add morbid gene annotation
    snv_df = snv_df.merge(morbidgenes_df, left_on='Gene', right_on='gene', how='left')

    # add gencc_morbid gene annotation
    snv_df = snv_df.merge(gencc_df, left_on='Gene', right_on='Gene', how='left')
    pd.set_option('display.max_columns', None)

    # determine if parent samples are present in sv vcf
    snv_df = check_for_parents(snv_df, snv_headerline, father, mother)

    # tag compound heterozygous variants
    get_comphetvars(sv_df, snv_df, out_file)


if __name__ == '__main__':
    main()
