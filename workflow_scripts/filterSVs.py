#!/usr/bin/env python

import sys
import pandas as pd
from peds import open_ped
import argparse
import numpy as np


def maxAF(colval):
    if colval != np.NAN and isinstance(colval, str) and colval != '':
        afvals = colval.split('&')
        afvals = [float(x) for x in afvals]
        if len(afvals) > 1:
            return max(afvals)
        elif len(afvals) == 1:
            return afvals[0]
    elif isinstance(colval, float):
        return colval
    else:
        return 0


def main():
    parser = argparse.ArgumentParser(description='Find variants overlapping differentially methylated regions', formatter_class=argparse.RawTextHelpFormatter)
    group1 = parser.add_argument_group('Mandatory inputs')
    group1.add_argument('-in_vcf', type=str, dest='in_vcf', required=True,
                        help='Input VCF file')
    group1.add_argument('-gnomad_cutoff', type=str, dest='gnomad_cutoff', required=True,
                        help='Maximum popmax_AF in gnomad')
    group1.add_argument('-hprc_cutoff', type=str, dest='hprc_cutoff', required=True,
                        help='Maximum popmax_AF in hprc controls')
    group1.add_argument('-colorsdb_cutoff', type=str, dest='colorsdb_cutoff', required=True,
                        help='Maximum popmax_AF in CoLoRSdb')
    group1.add_argument('-o', type=str, dest='out_vcf', required=True,
                        help='Output VCF file name')

    args = parser.parse_args()

    in_vcf = args.in_vcf
    gnomad_cutoff = args.gnomad_cutoff
    hprc_cutoff = args.hprc_cutoff
    colorsdb_cutoff = args.colorsdb_cutoff
    out_vcf = args.out_vcf

    # read in sv files
    with open(in_vcf, "r") as i:
        vcf_header = i.readlines()
    vcf_header = [line for line in vcf_header if line.startswith('#')]
    out_headerline = [line for line in vcf_header if line.startswith('#CHROM')][0]
    headerline = out_headerline.strip('\n').split('\t')
    vcf_df = pd.read_csv(in_vcf, sep='\t', comment="#", header=None, names=headerline)
    vcf_df['CSQ'] = vcf_df['INFO'].str.split('CSQ=').str[1].str.split(',').str[0]
    vcf_df['Gnomad_AF'] = vcf_df['INFO'].str.split('Best_gnomAD_PopMax_AF=').str[1].str.split(';').str[0]
    vcf_df['HPRC_AF'] = vcf_df['INFO'].str.split('Best_HPRC_AF=').str[1].str.split(';').str[0]
    vcf_df['CoLoRSdb_AF'] = vcf_df['INFO'].str.split('Best_CoLoRSdb_AF=').str[1].str.split(';').str[0]
    vcf_df.fillna({'Gnomad_AF': 0, 'HPRC_AF': 0, 'CoLoRSdb_AF': 0}, inplace=True)
    pd.set_option('display.max_columns', None)

    # filter variants
    # print(vcf_df.shape)
    out_df = vcf_df[(vcf_df['Gnomad_AF'].astype(float) <= float(gnomad_cutoff)) & (vcf_df['HPRC_AF'].astype(float) <= float(hprc_cutoff)) &
                    (vcf_df['CoLoRSdb_AF'].astype(float) <= float(colorsdb_cutoff))]

    # write tsv output
    out_df.to_csv(out_vcf+".tsv", sep='\t', index=False)
    out_df = out_df[headerline]
    # print(out_df.shape)

    #write vcf output
    if len(out_df) >= 1:
        with open(out_vcf, 'w') as of:
            of.writelines(vcf_header)
        out_df.to_csv(out_vcf, mode='a', sep='\t', index=False, header=False)



if __name__ == '__main__':
    main()
