#!/usr/bin/env python

import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description='Parses trgt output to identify potentially pathogenic repeat expansions', formatter_class=argparse.RawTextHelpFormatter)
    group1 = parser.add_argument_group('Mandatory inputs')
    group1.add_argument('-trgt_vcf', type=str, dest='trgt_vcf', required=True,
                        help='trgt.sorted.vcf VCF file')
    group1.add_argument('-repeatsdb', type=str, dest='repeatsdb', required=True,
                        help='database tsv file with pathogenic repeats')
    group1.add_argument('-bed', type=str, dest='bed', required=True,
                        help='Bed file with pathogenic repeat loci')
    group1.add_argument('-o', type=str, dest='out_tsv', required=True,
                        help='Output tsv file name')
    args = parser.parse_args()
    trgt_vcf = args.trgt_vcf
    repeatsdb = args.repeatsdb
    bed = args.bed
    out_tsv = args.out_tsv

    # load repeat database tsv
    repeatdf = pd.read_csv(repeatsdb, sep='\t', header=0)
    repeatdf = repeatdf.fillna('')

    # load sample vcf file
    with open(trgt_vcf, "r") as i:
        vcf_header = i.readlines()
    vcf_header = [line for line in vcf_header if line.startswith('#')]
    headerline = [line for line in vcf_header if line.startswith('#CHROM')][0]
    headerline = headerline.strip('\n').split('\t')
    infocol = headerline[-1]
    vcf_df = pd.read_csv(trgt_vcf, sep='\t', comment="#", header=None, names=headerline)
    vcf_df[['GT','AL','ALLR','SD','MC','MS','AP','AM']] = vcf_df[infocol].str.split(':', n=7, expand=True)
    vcf_df['TRID'] = vcf_df['INFO'].str.split('TRID=').str[1].str.split(';').str[0]
    vcf_df['END'] = vcf_df['INFO'].str.split('END=').str[1].str.split(';').str[0]
    vcf_df['MOTIFS'] = vcf_df['INFO'].str.split('MOTIFS=').str[1].str.split(';').str[0]
    vcf_df['STRUC'] = vcf_df['INFO'].str.split('STRUC=').str[1].str.split(';').str[0]
    vcf_df = vcf_df[vcf_df['AL'] != '.']

    with open(out_tsv, 'w') as o:
        o.write('Gene\tHaplotype_copies\tPathogenicity_Threshold\tPotentially_Pathogenic_Motif\tMotifs\tStructure\t'
                'Motif_Copies\tMotif_Support\tPotentially_Pathogenic\tDiseaseName\tDiseaseAlias\tDiseaseInheritance\t'
                'DiseaseEvidence\n')
        for i, r in vcf_df.iterrows():
            gene = r['TRID']
            subdf = repeatdf[repeatdf['g_symbol'] == gene]
            if len(subdf) >= 1:
                disease = list(subdf['d_name'])[0]
                disease_alias = list(subdf['d_name_alias'])[0]
                disease_inheritance = list(subdf['d_inheritance'])[0]
                disease_evidence = list(subdf['d_evidence'])[0]
                repeatnums = r['MC'].split(',')

                # check if pathogenic motif is present
                motifs = r['MOTIFS'].split(',')
                pathogenic_motif = list(subdf['pathogenic_motifs'])[0]
                path_repeats = []
                if pathogenic_motif in motifs:
                    pathindices = [i for i, s in enumerate(motifs) if s == pathogenic_motif]
                    for pathindex in pathindices:
                        path_repeats.append(int(repeatnums[0].split('_')[pathindex]))
                        try:
                            path_repeats.append(int(repeatnums[1].split('_')[pathindex]))
                        except:
                            pass
                path_threshold = list(subdf['pathogenic_low'])[0]

                # check if pathogenic motif repeats are above the pathogenic repeat threshold
                if any(x > path_threshold for x in path_repeats):
                    o.write(gene + '\t' + str(path_repeats) + '\t' + str(path_threshold)
                            + '\t' + pathogenic_motif + '\t' + r['MOTIFS'] + '\t' + r['STRUC'] + '\t' + r['MC'] + '\t'
                            + r['SD'] + '\t' + 'yes' + '\t' + disease + '\t' + disease_alias +
                            '\t' + disease_inheritance + '\t' + disease_evidence + '\n')
                else:
                    o.write(gene + '\t' + str(path_repeats) + '\t' + str(path_threshold)
                            + '\t' + pathogenic_motif + '\t' + str(r['MOTIFS']) + '\t' + str(r['STRUC']) + '\t' +
                            str(r['MC']) + '\t' + str(r['SD']) + '\t' + 'no' + '\t' + disease + '\t' + disease_alias +
                            '\t' + disease_inheritance + '\t' + disease_evidence + '\n')


if __name__ == '__main__':
    main()




if __name__ == '__main__':
    main()
