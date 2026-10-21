#!/usr/bin/env python

import sys
import pandas as pd
from peds import open_ped
import argparse
import numpy as np
import bioframe as bf
import re
import warnings



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

def filter_common_variants(df, af_thresholds=None, ofp_threshold=0.5):
    af_cols  = ["Gnomad_AF",       "HPRC_AF",       "CoLoRSdb_AF",       "CdLS_AF"]
    ofp_cols = ["Best_gnomAD_OFP", "Best_HPRC_OFP", "Best_CoLoRSdb_OFP", "Best_CdLS_OFP"]

    if af_thresholds is None:
        af_thresholds = {col: 0.005 for col in af_cols}

    fail_mask = pd.Series(False, index=df.index)
    for af_col, ofp_col in zip(af_cols, ofp_cols):
        threshold = af_thresholds.get(af_col, 0.005)
        high_af  = df[af_col].astype(float) > float(threshold)
        high_ofp = df[ofp_col].astype(float) > ofp_threshold
        fail_mask |= (high_af & high_ofp)

    return df[~fail_mask]

_INTRONIC_RE = re.compile(r"^intron\d+-intron\d+$")

def _annotsv_tag(ann_mode, loc):
    """Return location tag for one AnnotSV row."""
    """
        Adds one location-tag column to an AnnotSV TSV:

          annotsv_location_tag  — coarse tag derived from AnnotSV's Location column.

        Tags apply to split rows only; full rows always get "".

        Tag values and their meaning:
          whole_gene  — variant spans the entire transcript (txStart-txEnd)
          5prime      — variant starts at the transcript 5' boundary (txStart-*),
                        but does not span the full gene
          3prime      — variant ends at the transcript 3' boundary (*-txEnd),
                        but does not span the full gene
          intronic    — variant is entirely within one or more introns
                        (Location matches intronN-intronM pattern)
          intragenic  — all other split rows: internal to the gene, not at either
                        boundary, and not purely intronic (e.g. exon2-exon4,
                        exon2-intron3, intron3-exon4); intentionally broad —
                        may be sub-divided in future iterations
          "" (empty)  — full rows (summary rows; use split rows for location detail)
        """
    if ann_mode != "split":
        return ""
    loc = str(loc).strip()
    if loc == "txStart-txEnd":
        return "whole_gene"
    if loc.startswith("txStart-"):
        return "5prime"
    if loc.endswith("-txEnd"):
        return "3prime"
    if _INTRONIC_RE.match(loc):
        return "intronic"
    if loc:
        return "intragenic"
    return ""


def filter_repeatSV(out_df, repeatmaskerdb):
    repcolnames = ['bin', 'swScore', 'milliDiv', 'milliDel', 'milliIns', 'chrom', 'start', 'end', 'genoLeft',
                   'strand', 'repName', 'repClass', 'repFamily', 'repStart', 'repEnd', 'repLeft', 'id']
    repeatdf = pd.read_csv(repeatmaskerdb, sep='\t', header=None, names=repcolnames)
    repeatdf = repeatdf[['chrom', 'start', 'end', 'repName', 'repClass', 'repFamily']]
    # print(repeatdf[(repeatdf['chrom'] == 'chr1') & (repeatdf['start'] < 28189500) & (repeatdf['start'] > 28188584)])
    # print('No. of rows in original df: ', out_df.shape)
    out_df['chrom'] = out_df['#CHROM']
    out_df['start'] = out_df['POS'].astype(int)
    out_df['end'] = out_df['END'].astype(int)
    overlapping_intervals = bf.overlap(out_df, repeatdf, suffixes=('', '_rep'))
    # print(overlapping_intervals[overlapping_intervals['start'] == 28189163])
    out_df = overlapping_intervals[(overlapping_intervals['chrom_rep'].isna()) | (overlapping_intervals['LEN'] > 1000)]
    # print(out_df.head())
    # print('No. of rows in non-overlapping df: ', out_df.shape)
    return out_df

def extract_alt_pos(alt):
    """Return (chrom, pos) encoded in a BND ALT string."""
    m = re.search(r'[\[\]]([^\[\]:]+):(\d+)[\[\]]', str(alt))
    return (m.group(1), int(m.group(2))) if m else (None, None)

def extract_event(info_str):
    m = re.search(r'EVENT=([^;]+)', str(info_str))
    return m.group(1) if m else None

def update_info(info_str, end, svlen):
    """Rewrite INFO for a consolidated INV record."""
    s = re.sub(r'SVTYPE=BND', 'SVTYPE=INV', str(info_str))
    s = re.sub(r'MATEID=[^;]+;?', '', s)
    s = re.sub(r'HOMLEN=[^;]+;?', '', s)
    s = re.sub(r'HOMSEQ=[^;]+;?', '', s)
    s = f'END={end};SVLEN={svlen};' + s
    return re.sub(r';+', ';', s).strip(';')

def merge_nonempty(series):
    """Return first non-empty / non-null value from a pandas Series."""
    for v in series:
        if pd.notna(v) and str(v).strip() not in ('', 'nan', '.', 'None', 'NaN'):
            return v
    return series.iloc[0]

def consolidateinversions(inv_df, sample):
    inv_df['_event'] = inv_df['INFO'].apply(extract_event)
    inv_df['_mate_pos'] = inv_df['ALT'].apply(lambda x: extract_alt_pos(x)[1])
    # annot_cols = ['Gene', 'Consequence', 'Impact', 'HGVSc', 'HGVSp',
    #               'GenCC_Submission', 'Morbid_Gene', 'Gnomad_AF', 'HPRC_AF', 'CoLoRSdb_AF', 'CdLS_AF', 'ID']
    annot_cols = ['CytoBand','Gene_name','Closest_left','Closest_right','Gene_count', 'Tx', 'Overlapped_CDS_length',
                  'Overlapped_CDS_percent', 'Frameshift', 'Exon_count', 'RE_gene','P_gain_phen',
                  'P_gain_hpo','P_gain_source','P_gain_coord','P_loss_phen','P_loss_hpo','P_loss_source','P_loss_coord',
                  'P_ins_phen','P_ins_hpo','P_ins_source','P_ins_coord','po_P_gain_phen','po_P_gain_hpo',
                  'po_P_gain_source','po_P_gain_coord','po_P_gain_percent','po_P_loss_phen','po_P_loss_hpo',
                  'po_P_loss_source','po_P_loss_coord','po_P_loss_percent','P_snvindel_nb','P_snvindel_phen',
                  'B_gain_source','B_gain_coord','B_gain_AFmax','B_loss_source','B_loss_coord','B_loss_AFmax',
                  'B_ins_source','B_ins_coord','B_ins_AFmax','B_inv_source','B_inv_coord','B_inv_AFmax',
                  'po_B_gain_allG_source','po_B_gain_allG_coord','po_B_gain_someG_source','po_B_gain_someG_coord',
                  'po_B_loss_allG_source','po_B_loss_allG_coord','po_B_loss_someG_source','po_B_loss_someG_coord',
                  'GC_content_left','GC_content_right','Repeat_coord_left','Repeat_type_left','Repeat_coord_right',
                  'Repeat_type_right','Gap_left','Gap_right','SegDup_left','SegDup_right','ENCODE_blacklist_left',
                  'ENCODE_blacklist_characteristics_left','ENCODE_blacklist_right',
                  'ENCODE_blacklist_characteristics_right','ACMG','HI','TS','DDD_HI_percent','ExAC_delZ','ExAC_dupZ',
                  'ExAC_cnvZ','ExAC_synZ','ExAC_misZ','GenCC_disease','GenCC_moi','GenCC_classification','GenCC_pmid',
                  'NCBI_gene_ID','HGNC_gene_ID','OMIM_ID','OMIM_phenotype','OMIM_inheritance','OMIM_morbid',
                  'OMIM_morbid_candidate','LOEUF_bin','GnomAD_pLI','ExAC_pLI','AnnotSV_ranking_score',
                  'AnnotSV_ranking_criteria','ACMG_class', 'annotsv_location_tag', 'location_boundary_detail',
                  'cds_utr_context', 'nearest_junction', 'nearest_junction_dist']
    # annot_cols = ['Gene', 'Consequence', 'Impact', 'HGVSc', 'HGVSp',
    #               'GenCC_Submission', 'Morbid_Gene']
    segdup_cols = ['segdup_coordinates_segdup', 'segdup_OtherCoordinates_segdup',
                   'fracMatch_segdup', 'repClass_rep', 'repFamily_rep']
    sample_col = sample
    
    rows_out = []
    for event, grp in inv_df.groupby('_event', sort=False):
        grp = grp.sort_values('SV_start').reset_index(drop=True)
        row = grp.iloc[0].copy()
        # print(event, grp)

        # Inversion span: lowest POS → highest of all POS + mate positions
        all_pos = list(grp['SV_start']) + [p for p in grp['_mate_pos'] if pd.notna(p)]
        # print(all_pos)

        span_start = int(min(all_pos))
        span_end = int(max(all_pos))
        svlen = span_end - span_start

        row['SV_start'] = span_start
        row['SV_end'] = span_end
        row['LEN'] = svlen
        row['SVTYPE'] = 'INV'
        row['ALT'] = '<INV>'
        row['REF'] = 'N'
        row['INV_ID'] = event

        # Genotype: use the row with highest QUAL
        best = grp.loc[grp['QUAL'].idxmax()]
        row[sample_col] = best[sample_col]
        row['QUAL'] = best['QUAL']
        row['FORMAT'] = best['FORMAT']


        for col in annot_cols:
            if col in best:
                row[col] = merge_nonempty(grp[col])
        # Segdup: merge unique values across rows
        # for col in segdup_cols:
        #     # print(col, grp[col])
        #     vals = [str(v) for v in grp[col]]
        #             # if pd.notna(v) and str(v) not in ('[None]', 'nan', 'None')]
        #     unique_vals = list(dict.fromkeys(vals))   # deduplicate, preserve order
        #     row[col] = ' | '.join(unique_vals) if unique_vals else None


        # Update INFO
        row['INFO'] = update_info(best['INFO'], span_end, svlen)

        rows_out.append(row)

    out = pd.DataFrame(rows_out).drop(columns=['_event', '_mate_pos'])
    return out

def main():
    parser = argparse.ArgumentParser(description='Find variants overlapping differentially methylated regions', formatter_class=argparse.RawTextHelpFormatter)
    group1 = parser.add_argument_group('Mandatory inputs')
    group1.add_argument('-in_tsv', type=str, dest='in_tsv', required=True,
                        help='Input TSV file from AnnotSV output')
    group1.add_argument('-proband', type=str, dest='proband', required=True,
                        help='Output inversions file name')
    group1.add_argument('-gnomad_cutoff', type=str, dest='gnomad_cutoff', required=True,
                        help='Maximum popmax_AF in gnomad')
    group1.add_argument('-hprc_cutoff', type=str, dest='hprc_cutoff', required=True,
                        help='Maximum popmax_AF in hprc controls')
    group1.add_argument('-colorsdb_cutoff', type=str, dest='colorsdb_cutoff', required=True,
                        help='Maximum popmax_AF in CoLoRSdb')
    group1.add_argument('-cdls_cutoff', type=str, dest='cdls_cutoff', required=True,
                        help='Maximum popmax_AF in control CdLS cohort')
    # group1.add_argument('-o_vcf', type=str, dest='out_vcf', required=True,
    #                     help='Output VCF file name')
    group1.add_argument('-o_tsv', type=str, dest='out_tsv', required=True,
                        help='Output VCF file name')
    group1.add_argument('-o_inv', type=str, dest='out_inv', required=True,
                        help='Output inversions file name')

    args = parser.parse_args()

    in_tsv = args.in_tsv
    gnomad_cutoff = args.gnomad_cutoff
    hprc_cutoff = args.hprc_cutoff
    colorsdb_cutoff = args.colorsdb_cutoff
    cdls_cutoff = args.cdls_cutoff
    # out_vcf = args.out_vcf
    out_tsv = args.out_tsv
    out_inv = args.out_inv
    sample = args.proband


    # read in sv files
    vcf_df = pd.read_csv(in_tsv, sep='\t', header=0)
    vcf_df['SVTYPE'] = vcf_df['INFO'].str.split('SVTYPE=').str[1].str.split(';').str[0]
    vcf_df['LEN'] = vcf_df['INFO'].str.split('SVLEN=').str[1].str.split(';').str[0].astype('Int64').abs()
    vcf_df['Gnomad_AF'] = vcf_df['INFO'].str.split('Best_gnomAD_PopMax_AF=').str[1].str.split(';').str[0]
    vcf_df['Best_gnomAD_OFP'] = vcf_df['INFO'].str.split('Best_gnomAD_OFP=').str[1].str.split(';').str[0]
    vcf_df['HPRC_AF'] = vcf_df['INFO'].str.split('Best_HPRC_AF=').str[1].str.split(';').str[0]
    vcf_df['Best_HPRC_OFP'] = vcf_df['INFO'].str.split('Best_HPRC_OFP=').str[1].str.split(';').str[0]
    vcf_df['CoLoRSdb_AF'] = vcf_df['INFO'].str.split('Best_CoLoRSdb_AF=').str[1].str.split(';').str[0]
    vcf_df['Best_CoLoRSdb_OFP'] = vcf_df['INFO'].str.split('Best_CoLoRSdb_OFP=').str[1].str.split(';').str[0]
    vcf_df['CdLS_AF'] = vcf_df['INFO'].str.split('Best_CdLS_AF=').str[1].str.split(';').str[0]
    vcf_df['Best_CdLS_OFP'] = vcf_df['INFO'].str.split('Best_CdLS_OFP=').str[1].str.split(';').str[0]
    fillnandict = {'Gnomad_AF': 0, 'HPRC_AF': 0, 'CoLoRSdb_AF': 0, 'CdLS_AF': 0,
                   'Best_gnomAD_OFP': 0, 'Best_HPRC_OFP': 0, 'Best_CoLoRSdb_OFP': 0, 'Best_CdLS_OFP': 0}
    vcf_df.fillna(fillnandict, inplace=True)
    pd.set_option('display.max_columns', None)


    
    # inv_df = out_df[(out_df['FILTER'] == 'InvBreakpoint') & (out_df['Annotation_mode'] == 'split')]
    inv_df = vcf_df[(vcf_df['FILTER'] == 'InvBreakpoint')]
    # # print(inv_df.head(5))
    inv_df = consolidateinversions(inv_df, sample)
    conditions = [
        (((inv_df["annotsv_location_tag"] != 'whole_gene') & (inv_df["Gene_name"] != " ") & (inv_df['Overlapped_CDS_percent'] > 0)) | (inv_df['Frameshift'] == 'yes'))
    ]
    choices = ["yes"]
    # Syntax: np.select(conditions_list, choices_list, default_value)
    inv_df['gene_disrupting'] = np.select(conditions, choices, default="")
    inv_df = inv_df[['SV_chrom','SV_start','SV_end','LEN','SVTYPE','REF','ALT',sample,'Gene_name','annotsv_location_tag',
                     'GenCC_disease','GenCC_moi','GenCC_classification','GenCC_pmid','NCBI_gene_ID',
                     'OMIM_ID','OMIM_phenotype','OMIM_inheritance','OMIM_morbid','OMIM_morbid_candidate',
                     'CytoBand', 'location_boundary_detail','cds_utr_context','nearest_junction','nearest_junction_dist',
                     'Tx','Overlapped_CDS_length','Overlapped_CDS_percent','Frameshift','Exon_count', 'RE_gene',
                     'Gnomad_AF','Best_gnomAD_OFP','HPRC_AF','Best_HPRC_OFP','CoLoRSdb_AF','Best_CoLoRSdb_OFP','CdLS_AF','Best_CdLS_OFP',
                     'gene_disrupting', 'INFO', 'INV_ID']]
    inv_df.to_csv(out_inv, sep='\t', index=False) # write inversions file

    # filter variants
    af_thresholds = {
        "Gnomad_AF": gnomad_cutoff,
        "HPRC_AF": hprc_cutoff,
        "CoLoRSdb_AF": colorsdb_cutoff,
        "CdLS_AF": cdls_cutoff,
    }
    out_df = filter_common_variants(vcf_df, af_thresholds=af_thresholds, ofp_threshold=0.5)
    out_df = out_df[out_df['FILTER'] == 'PASS']
    out_df.to_csv(out_tsv, sep='\t', index=False) # write filtered SVs file


if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    main()
