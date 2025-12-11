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


def get_comphetvars_singleton(vcf_df, headerline):
    # group by gene and look for comphet in singleton
    outdf = pd.DataFrame
    genelist = set(vcf_df['Gene'])
    counter = 0
    for gene in genelist:
        # gene = 'LIG3'
        gene_df = vcf_df[(vcf_df['Gene'] == gene)]
        if len(gene_df) > 1:
            genotypes_list = list(gene_df['Proband_GT'])
            # print(gene_df)
            # for phased variants
            if ('1|0' in genotypes_list) & ('0|1' in genotypes_list):
                # print(gene_df)
                counter += 1
                vardf = gene_df[(gene_df['Proband_GT'] == '0|1') | (gene_df['Proband_GT'] == '1|0')].reset_index()
                if counter == 1:
                    outdf = vardf
                else:
                    outdf = pd.concat([outdf, vardf], ignore_index=True)

            # for unphased variants
            unphasedhet = [x for x in genotypes_list if x == '0/1']
            if len(unphasedhet) >= 2:
                counter += 1
                vardf = gene_df[gene_df['Proband_GT'] == '0/1'].reset_index()
                if counter == 1:
                    outdf = vardf
                else:
                    outdf = pd.concat([outdf, vardf], ignore_index=True)

    # write output
    outdf.drop_duplicates(keep='first')
    outdf['INFO'] = outdf['INFO'].astype(str) + ';Compound_Heterozygous'
    outdf['keycol'] = outdf['#CHROM'].astype(str) + outdf['POS'].astype(str) + outdf['ALT'].astype(str)
    vcf_df['keycol'] = vcf_df['#CHROM'].astype(str) + vcf_df['POS'].astype(str) + vcf_df['ALT'].astype(str)
    new_vcf_df = vcf_df[~vcf_df['keycol'].isin(list(outdf['keycol']))]
    outdf = pd.concat([outdf, new_vcf_df], ignore_index=True)
    outdf = outdf[headerline]
    outdf['INFO'] = outdf['INFO'].map(lambda x: x.replace('"', ''))
    return outdf


def get_comphetvars_trio(vcf_df, headerline):
    # group by gene and look for comphet in trio
    outdf = pd.DataFrame
    genelist = set(vcf_df['Gene'])
    counter = 0
    # print(vcf_df)
    # vcf_df['CompHetStatus'] = ''
    # genelist = ['VPS45']
    for gene in genelist:
        # print(gene)
        gene_df = vcf_df[(vcf_df['Gene'] == gene)]
        # print(gene_df)
        if len(gene_df) > 1:
            proband_genotypes = list(gene_df['Proband_GT'])
            # for phased variants
            if ('1|0' in proband_genotypes) & ('0|1' in proband_genotypes):
                counter += 1
                vardf = gene_df[(gene_df['Proband_GT'] == '0|1') |
                                (gene_df['Proband_GT'] == '1|0')].reset_index()
                if counter == 1:
                    outdf = vardf
                else:
                    outdf = pd.concat([outdf, vardf], ignore_index=True)

            # for unphased variants
            else:
                gene_df = gene_df[(gene_df['Proband_GT'] == '0/1') |
                                  (gene_df['Proband_GT'] == '1/0')]
                mom_genotypes = list(gene_df['Mom_GT'])
                dad_genotypes = list(gene_df['Dad_GT'])
                dadinherited = gene_df[
                    ((gene_df['Dad_GT'] == '0/1') | (gene_df['Dad_GT'] == '1/0'))
                    & (gene_df['Mom_GT'] == '0/0')]
                mominherited = gene_df[
                    ((gene_df['Mom_GT'] == '0/1') | (gene_df['Mom_GT'] == '1/0'))
                    & (gene_df['Dad_GT'] == '0/0')]

                # at least one het variant inherited from each parent
                if ('0/1' in mom_genotypes) and ('0/1' in dad_genotypes):
                    # print('unphased 2 hets', gene_df['INFO'])
                    counter += 1
                    vardf = gene_df[(gene_df['Mom_GT'] == '0/1') | (gene_df['Dad_GT'] == '0/1')].reset_index()
                    if counter == 1:
                        outdf = vardf
                    else:
                        outdf = pd.concat([outdf, vardf], ignore_index=True)
                    # print(outdf)

                # 2 het variants present in proband - 1 inherited, 1 or more denovo
                if len(gene_df) > (len(dadinherited)+len(mominherited)) >= 1:
                    counter += 1
                    if counter == 1:
                        outdf = gene_df.reset_index()
                    else:
                        outdf = pd.concat([outdf, gene_df.reset_index()], ignore_index=True)
                    # print(outdf)

                # 2 het variants - both denovo
                if len(gene_df[(gene_df['Dad_GT'] == '0/0') & (gene_df['Mom_GT'] == '0/0')]) >= 2:
                    vardf = gene_df[(gene_df['Dad_GT'] == '0/0') & (gene_df['Mom_GT'] == '0/0')]
                    counter += 1
                    if counter == 1:
                        outdf = vardf
                    else:
                        outdf = pd.concat([outdf, vardf], ignore_index=True)

    # write output
    outdf.drop_duplicates(keep='first')
    outdf['INFO'] = outdf['INFO'].astype(str) + ';Compound_Heterozygous'
    outdf['keycol'] = outdf['#CHROM'].astype(str) + outdf['POS'].astype(str) + outdf['ALT'].astype(str)
    vcf_df['keycol'] = vcf_df['#CHROM'].astype(str) + vcf_df['POS'].astype(str) + vcf_df['ALT'].astype(str)
    new_vcf_df = vcf_df[~vcf_df['keycol'].isin(list(outdf['keycol']))]
    outdf = pd.concat([outdf, new_vcf_df], ignore_index=True)
    outdf = outdf[headerline]
    outdf['INFO'] = outdf['INFO'].map(lambda x: x.replace('"', ''))
    return outdf


def get_comphetvars_duo(vcf_df, headerline):
    outdf = pd.DataFrame
    genelist = set(vcf_df['Gene'])
    counter = 0
    for gene in genelist:
        gene_df = vcf_df[(vcf_df['Gene'] == gene)]
        # print(gene_df)
        if len(gene_df) > 1:
            proband_genotypes = list(gene_df['Proband_GT'])

            # for phased variants
            if ('1|0' in proband_genotypes) & ('0|1' in proband_genotypes):
                counter += 1
                vardf = gene_df[(gene_df['Proband_GT'] == '0|1') |
                                (gene_df['Proband_GT'] == '1|0')].reset_index()
                if counter == 1:
                    outdf = vardf
                else:
                    outdf = pd.concat([outdf, vardf], ignore_index=True)

            # for unphased variants
            else:
                gene_df = gene_df[(gene_df['Proband_GT'] == '0/1') |
                                  (gene_df['Proband_GT'] == '1/0')]
                inheritedvars = gene_df[(gene_df['Parent_GT'] == '0/1') | (gene_df['Parent_GT'] == '1/0')]

                # 2 het variants parent has 1 or none - 1/2 de novo or from missing parent
                if len(gene_df) > len(inheritedvars):
                    counter += 1
                    if counter == 1:
                        outdf = gene_df.reset_index()
                    else:
                        outdf = pd.concat([outdf, gene_df.reset_index()], ignore_index=True)

    # write output
    outdf.drop_duplicates(keep='first')
    outdf['INFO'] = outdf['INFO'].astype(str) + ';Compound_Heterozygous'
    outdf['keycol'] = outdf['#CHROM'].astype(str) + outdf['POS'].astype(str) + outdf['ALT'].astype(str)
    vcf_df['keycol'] = vcf_df['#CHROM'].astype(str) + vcf_df['POS'].astype(str) + vcf_df['ALT'].astype(str)
    new_vcf_df = vcf_df[~vcf_df['keycol'].isin(list(outdf['keycol']))]
    outdf = pd.concat([outdf, new_vcf_df], ignore_index=True)
    outdf = outdf[headerline]
    outdf['INFO'] = outdf['INFO'].map(lambda x: x.replace('"', ''))
    return outdf


def main():
    parser = argparse.ArgumentParser(description='Assign extended comphets', formatter_class=argparse.RawTextHelpFormatter)
    group1 = parser.add_argument_group('Mandatory inputs')
    group1.add_argument('-i', type=str, dest='in_file', required=True,
                        help='Input vcf file in unzipped format')
    group1.add_argument('-o', type=str, dest='out_file', required=True,
                        help='Output vcf file in unzipped format')
    group1.add_argument('-ped', type=str, dest='ped_file', required=True,
                        help='Pedigree file')
    group1.add_argument('-proband', type=str, dest='proband_ID', required=True,
                        help='ID of proband')

    args = parser.parse_args()

    in_file = args.in_file
    out_file = args.out_file
    in_ped = args.ped_file
    proband_ID = args.proband_ID

    # get proband, mom, dad from ped file
    proband, father, mother = read_pedfile(in_ped, proband_ID)

    # parse vcf file
    with open(in_file, "r") as i:
        vcf_header = i.readlines()
    vcf_header = [line for line in vcf_header if line.startswith('#')]
    out_header = vcf_header[0:41]
    # print(out_header)
    out_headerline = [line for line in vcf_header if line.startswith('#CHROM')][0]
    headerline = out_headerline.strip('\n').split('\t')

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
    vcf_df = pd.read_csv(in_file, sep='\t', comment="#", header=None, names=headerline)
    vcf_df = vcf_df[vcf_df['FILTER'] == 'PASS'].reset_index()  # for variants with PASS only
    vcf_df['CSQ'] = vcf_df['INFO'].str.split('CSQ=').str[1]
    vcf_df['Gene'] = vcf_df['CSQ'].str.split('|').str[3]
    genelist = set(vcf_df['Gene'])
    print(str(len(genelist)) + ' genes found')

    # pull genotypes out as columns for analysis
    vcf_df['Proband_GT'] = vcf_df[proband_ID].str.split(':').str[0]
    if sampletype == 'trio':
        vcf_df['Dad_GT'] = vcf_df[father.id].str.split(':').str[0]
        vcf_df['Mom_GT'] = vcf_df[mother.id].str.split(':').str[0]
    elif sampletype == 'duo_dad':
        vcf_df['Parent_GT'] = vcf_df[father.id].str.split(':').str[0]
    elif sampletype == 'duo_mom':
        vcf_df['Parent_GT'] = vcf_df[mother.id].str.split(':').str[0]

    # tag compound heterozygous variants
    if sampletype == 'singleton':
        comphetdf = get_comphetvars_singleton(vcf_df, headerline)  # done
    elif sampletype == 'trio':
        comphetdf = get_comphetvars_trio(vcf_df, headerline)  # done
    elif sampletype == 'duo_mom' or sampletype == 'duo_dad':
        comphetdf = get_comphetvars_duo(vcf_df, headerline)  # done
    else:
        comphetdf = pd.DataFrame

    # write output vcf file
    if len(comphetdf) >= 1:
        with open(out_file, 'w') as of:
            of.writelines(out_header)
            of.writelines(out_headerline)
        comphetdf.to_csv(out_file, mode='a', sep='\t', index=False, header=False)
    print('Results written to ' + out_file)


if __name__ == '__main__':
    main()
