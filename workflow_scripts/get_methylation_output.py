#!/usr/bin/env python

import sys
import pandas as pd
from peds import open_ped
import argparse
import numpy as np
from pybedtools import BedTool


def max_methstate(uncat, meth, unmeth, allelespec):
    if len({uncat, meth, unmeth, allelespec}) == 1:
        return np.NAN
    elif max(uncat, meth, unmeth, allelespec) == uncat:
        return 'Uncategorized'
    elif max(uncat, meth, unmeth, allelespec) == meth:
        return 'Methylated'
    elif max(uncat, meth, unmeth, allelespec) == unmeth:
        return 'Unmethylated'
    elif max(uncat, meth, unmeth, allelespec) == allelespec:
        return 'AlleleSpecificMethylation'
    else:
        return np.NAN


def main():
    parser = argparse.ArgumentParser(description='Find variants overlapping differentially methylated regions', formatter_class=argparse.RawTextHelpFormatter)
    group1 = parser.add_argument_group('Mandatory inputs')
    group1.add_argument('-case', type=str, dest='case', required=True,
                        help='Sample Methylation Profile')
    group1.add_argument('-control', type=str, dest='control', required=True,
                        help='Control Methylation Profile')
    group1.add_argument('-sv', type=str, dest='sv_vcf', required=False, default=" ",
                        help='vcf file containing proband SVs')
    group1.add_argument('-snv', type=str, dest='snv_vcf', required=True,
                        help='vcf file containing proband SNVs')
    group1.add_argument('-o', type=str, dest='out_file', required=True,
                        help='Output file name')

    args = parser.parse_args()

    case = args.case
    control = args.control
    out_file = args.out_file
    sv_vcf = args.sv_vcf
    snv_vcf = args.snv_vcf

    # read and preprocess control profile
    headercols = ['chrom', 'start', 'end', 'data_category', 'num_phased', 'num_unphased', 'NoData', 'Uncategorized',
                  'Methylated', 'Unmethylated', 'AlleleSpecificMethylation', 'avg_abs_meth_deltas', 'stdev_abs_meth_deltas',
                  'avg_combined_methyls', 'stdev_combined_methyls', 'chr', 'st', 'en', 'Gene', 'ExpressedAllele', 'Overlap']
    controlprofile_df = pd.read_csv(control, sep='\t', header=None, names=headercols)
    controlprofile_df['Control_Status'] = controlprofile_df.apply(lambda x: max_methstate(x.Uncategorized, x.Methylated,
                                                                       x.Unmethylated, x.AlleleSpecificMethylation),
                                            axis=1)

    # read case profile
    headercols = ['chrom', 'start', 'end', 'case_summary_label']
    caseprofile_df = pd.read_csv(case, sep= '\t', usecols=[0,1,2,3], names=headercols)

    # find differentially methylated regions
    merged_df = pd.merge(left=caseprofile_df, right=controlprofile_df, how='left', left_on=['chrom', 'start', 'end'], right_on=['chrom', 'start', 'end'])
    merged_df = merged_df[['chrom', 'start', 'end', 'Control_Status', 'case_summary_label', 'Gene']]
    merged_df = merged_df[merged_df['Control_Status'] != merged_df['case_summary_label']]
    print(merged_df.shape)
    missingdatalabels = ['Uncategorized', 'NoData']
    merged_df = merged_df[~merged_df['Control_Status'].isin(missingdatalabels) & ~merged_df['case_summary_label'].isin(missingdatalabels)]
    print(merged_df.shape)
    pd.set_option('display.max_columns', None)
    merged_df.to_csv('temp_differentiallymethylatedregions.bed', sep='\t', header=False, index=False)

    # intersect with case SNVs and SVs
    dmrs = BedTool('temp_differentiallymethylatedregions.bed')
    outlist = []

    if sv_vcf != " ":   # if SV file is provided
        svs = BedTool(sv_vcf)
        nearby_SVs = svs.closest(dmrs, d=True, stream=True)
        for var in nearby_SVs:
            if 0 <= int(var[-1]) < 1000:
                # print(var.fields)
                outlist.append(var)

    snvs = BedTool(snv_vcf)
    nearby_snvs = snvs.closest(dmrs, d=True, stream=True)
    for var in nearby_snvs:
        if 0 <= int(var[-1]) < 1000:
            # print(var.fields)
            outlist.append(var)

    # write to outfile
    with open(out_file, 'w') as o:
        o.write("Chrom\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\tDMR_chrom\tDMR_start\tDMR_end\t"
                "DMR_Control_Status\tDMR_Case_Summary_Label\tDMR_Gene\tDistance_To_Gene\n")
        for var in outlist:
            o.write(str(var))


if __name__ == '__main__':
    main()
