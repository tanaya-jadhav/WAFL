#!/usr/bin/env python

import pandas as pd
import argparse


def main():
    parser = argparse.ArgumentParser(
        description='Compiles output reports into an easily readable excel file',
        formatter_class=argparse.RawTextHelpFormatter)
    group1 = parser.add_argument_group('Mandatory inputs')
    group1.add_argument('-D_tsv', type=str, dest='D_tsv', required=True,
                        help='TSV file containing Dominant Variants')
    group1.add_argument('-R_tsv', type=str, dest='R_tsv', required=True,
                        help='TSV file containing Recessive Variants')
    group1.add_argument('-snv_sv', type=str, dest='snv_sv', required=True,
                        help='TSV file containing SV + SNV comphet variants')
    group1.add_argument('-SV', type=str, dest='SV', required=True,
                        help='TSV file containing filtered structural variants')
    group1.add_argument('-DMR', type=str, dest='dmr', required=True,
                        help='TSV file containing variants in DMRs')
    group1.add_argument('-TR', type=str, dest='tr', required=True,
                        help='TSV file containing the repeat expansion report')
    group1.add_argument('-o', type=str, dest='out', required=True,
                        help='name of output .xlsx file')

    args = parser.parse_args()
    D_tsv = args.D_tsv
    R_tsv = args.R_tsv
    snv_sv = args.snv_sv
    SV = args.SV
    dmr = args.dmr
    tr = args.tr
    outfile = args.out

    files = [D_tsv, R_tsv, snv_sv, SV, dmr, tr]
    sheets = ['Dominant Variants', 'Recessive Variants', 'SNV+SV CompHets', 'Filtered SVs', 'Variants in DMR',
              'Repeat Expansion Report']
    with pd.ExcelWriter(outfile, engine='openpyxl') as writer:
        for i, file in enumerate(files):
            # print(i)
            df = pd.read_csv(file, sep='\t', header=0)
            df.to_excel(writer, sheet_name=sheets[i], index=False)


if __name__ == '__main__':
    main()