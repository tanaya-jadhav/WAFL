#!/usr/bin/env python3
"""
match_bnd_ins.py - Match Sawfish BND breakends against HitSV INS insertions

Usage:
    python match_bnd_ins.py --sawfish <sawfish.vcf.gz> --hitsv <hitsv.vcf.gz> \
        [--output bnd_ins_matches.tsv] [--threshold 50]
"""
import argparse
import bisect
import csv
import re
import sys
from collections import defaultdict
import pandas as pd
import bioframe as bf
import warnings
import numpy as np

try:
    import cyvcf2
    _BACKEND = 'cyvcf2'
except ImportError:
    try:
        import pysam
        _BACKEND = 'pysam'
    except ImportError:
        raise RuntimeError("Neither cyvcf2 nor pysam is available.")

_BND_RE = re.compile(r'[\[\]]([^[\]:]+):(\d+)[\[\]]')

TSV_HEADER = [
    'sawfish_id', 'sawfish_chr', 'sawfish_pos',
    'sawfish_mate_chr', 'sawfish_mate_pos',
    'sawfish_filter', 'sawfish_homlen', 'sawfish_inslen',
    'hitsv_id', 'hitsv_chr', 'hitsv_pos', 'hitsv_svlen', 'hitsv_filter',
    'distance', 'match_end', 'Gene', 'Consequence', 'Impact',
    'GenCC_Submission', 'Morbid_Gene', 'Gnomad_AF', 'HPRC_AF', 'CoLoRSdb_AF', 'CdLS_AF'
]


def parse_bnd_alt(alt_str):
    m = _BND_RE.search(alt_str)
    # print(m.groups())
    if m:
        return m.group(1), int(m.group(2))
    return None, None


def load_hitsv_ins_cyvcf2(vcf_path):
    pos_index = defaultdict(list)
    rec_index = defaultdict(list)
    vcf = cyvcf2.VCF(vcf_path)
    for v in vcf:
        if v.INFO.get('SVTYPE') != 'INS':
            continue
        chrom = v.CHROM
        filt = v.FILTER if v.FILTER is not None else 'PASS'
        pos_index[chrom].append(v.POS)
        rec_index[chrom].append((v.ID, chrom, v.POS, v.INFO.get('SVLEN'), filt))
    vcf.close()
    return pos_index, rec_index


def load_hitsv_ins_pysam(vcf_path):
    pos_index = defaultdict(list)
    rec_index = defaultdict(list)
    vcf = pysam.VariantFile(vcf_path)
    for v in vcf.fetch():
        if v.info.get('SVTYPE') != 'INS':
            continue
        chrom = v.chrom
        filters = list(v.filter.keys())
        filt = 'PASS' if not filters or filters == ['PASS'] else ';'.join(filters)
        pos_index[chrom].append(v.pos + 1)  # pysam is 0-based
        svlen = v.info.get('SVLEN')
        if isinstance(svlen, tuple):
            svlen = svlen[0]
        rec_index[chrom].append((v.id, chrom, v.pos + 1, svlen, filt))
    vcf.close()
    return pos_index, rec_index


def load_hitsv_ins(vcf_path):
    if _BACKEND == 'cyvcf2':
        return load_hitsv_ins_cyvcf2(vcf_path)
    return load_hitsv_ins_pysam(vcf_path)


def find_nearby(pos_index, rec_index, chrom, query_pos, threshold):
    if chrom not in pos_index:
        return []
    positions = pos_index[chrom]
    records = rec_index[chrom]
    lo = bisect.bisect_left(positions, query_pos - threshold)
    hi = bisect.bisect_right(positions, query_pos + threshold)
    return [(records[i], abs(positions[i] - query_pos)) for i in range(lo, hi)]


def iter_bnd_records_cyvcf2(vcf_path):
    seen_pairs = set()
    vcf = cyvcf2.VCF(vcf_path)
    for v in vcf:
        if v.INFO.get('SVTYPE') != 'BND':
            continue
        vid = v.ID
        mate_id = v.INFO.get('MATEID')
        # Dedup: prefer :0 suffix; fallback to seen_pairs set
        if vid.endswith(':0'):
            pass
        elif vid.endswith(':1'):
            continue
        else:
            if vid in seen_pairs:
                continue
            seen_pairs.add(vid)
            if mate_id:
                seen_pairs.add(mate_id)

        filt = v.FILTER if v.FILTER is not None else 'PASS'
        alt = v.ALT[0] if v.ALT else ''
        # print(vid, v.ALT, alt)
        mate_chr, mate_pos = parse_bnd_alt(alt)
        # print(mate_chr, mate_pos)
        yield {
            'id': vid,
            'chr': v.CHROM,
            'pos': v.POS,
            'filter': filt,
            'homlen': v.INFO.get('HOMLEN') or 0,
            'inslen': v.INFO.get('INSLEN') or 0,
            'mate_chr': mate_chr,
            'mate_pos': mate_pos,
        }
    vcf.close()


def iter_bnd_records_pysam(vcf_path):
    seen_pairs = set()
    vcf = pysam.VariantFile(vcf_path)
    for v in vcf.fetch():
        if v.info.get('SVTYPE') != 'BND':
            continue
        vid = v.id
        mate_id = v.info.get('MATEID')
        if vid.endswith(':0'):
            pass
        elif vid.endswith(':1'):
            continue
        else:
            if vid in seen_pairs:
                continue
            seen_pairs.add(vid)
            if mate_id:
                seen_pairs.add(mate_id)

        filters = list(v.filter.keys())
        filt = 'PASS' if not filters or filters == ['PASS'] else ';'.join(filters)
        alts = v.alts or ('',)
        alt = alts[0]
        mate_chr, mate_pos = parse_bnd_alt(alt)
        homlen = v.info.get('HOMLEN') or 0
        inslen = v.info.get('INSLEN') or 0
        yield {
            'id': vid,
            'chr': v.chrom,
            'pos': v.pos + 1,
            'filter': filt,
            'homlen': homlen,
            'inslen': inslen,
            'mate_chr': mate_chr,
            'mate_pos': mate_pos,
        }
    vcf.close()


def iter_bnd_records(vcf_path):
    if _BACKEND == 'cyvcf2':
        yield from iter_bnd_records_cyvcf2(vcf_path)
    else:
        yield from iter_bnd_records_pysam(vcf_path)


def run_matching(args, sawfishdf):
    print(f"Loading HitSV INS records from {args.hitsv} ...", file=sys.stderr)
    pos_index, rec_index = load_hitsv_ins(args.hitsv)
    total_ins = sum(len(v) for v in rec_index.values())
    print(f"  Loaded {total_ins:,} INS records.", file=sys.stderr)

    print(f"Matching Sawfish BNDs from {args.sawfish} ...", file=sys.stderr)

    total_bnd_pairs = 0
    matched_bnd_pairs = set()
    primary_matched_bnd = set()
    matched_ins_ids = set()
    total_rows = 0
    primary_rows = 0
    mate_rows = 0

    sawfishbnd = sawfishdf[sawfishdf['SVTYPE'] == 'BND']
    # print(sawfishbnd.head(5))

    # for bnd in iter_bnd_records(args.sawfish):
    outrows = []
    for bnd_i, bnd_r in sawfishbnd.iterrows():
        total_bnd_pairs += 1

        vid = bnd_r['ID']
        # if vid == "sawfish:0:15010:0:0:1":
            # print(bnd_r)
        mate_id = bnd_r['INFO'].split('MATEID=')[1].split(';')[0]
        seen_pairs = set()
        if 'HOMLEN=' in bnd_r['INFO']:
            homlen = bnd_r['INFO'].split('HOMLEN=')[1].split(';')[0]
        else:
            homlen = 0
        if 'INSLEN=' in bnd_r['INFO']:
            inslen = bnd_r['INFO'].split('INSLEN=')[1].split(';')[0]
        else:
            inslen = 0

        alt = bnd_r['ALT'] if bnd_r['ALT'] else ''
        mate_chr, mate_pos = parse_bnd_alt(alt)

        primary_hits = find_nearby(
            pos_index, rec_index,
            bnd_r['SV_chrom'], bnd_r['SV_start'], args.threshold
        )

        mate_hits = []
        if mate_chr is not None:
            mate_hits = find_nearby(
                pos_index, rec_index,
                mate_chr, mate_pos, args.threshold
            )

        for (ins_rec, dist) in primary_hits:
            ins_id, ins_chr, ins_pos, ins_svlen, ins_filt = ins_rec
            outdict = {}
            outdict = {'sawfish_id': vid, 'sawfish_chr': bnd_r['SV_chrom'], 'sawfish_pos': bnd_r['SV_start'],
                       'REF': bnd_r['REF'], 'ALT': bnd_r['ALT'],
                       'sawfish_mate_id':mate_id, 'sawfish_mate_chr': mate_chr, 'sawfish_mate_pos': mate_pos,
                       'sawfish_filter': bnd_r['FILTER'], 'sawfish_homlen': homlen, 'sawfish_inslen': inslen,
                       'hitsv_id': ins_id, 'hitsv_chr': ins_chr, 'hitsv_pos': ins_pos, 'hitsv_svlen': ins_svlen,
                       'hitsv_filter': ins_filt, 'distance': dist, 'match_end': 'primary',
                       args.proband: bnd_r[args.proband], 'Gene_name': bnd_r['Gene_name'],
                       'Frameshift': bnd_r['Frameshift'], 'GenCC_disease': bnd_r['GenCC_disease'],
                       'GenCC_moi': bnd_r['GenCC_moi'], 'GenCC_classification': bnd_r['GenCC_classification'],
                       'GenCC_pmid': bnd_r['GenCC_pmid'], 'OMIM_ID': bnd_r['OMIM_ID'],
                       'OMIM_phenotype': bnd_r['OMIM_phenotype'], 'OMIM_inheritance': bnd_r['OMIM_inheritance'],
                       'OMIM_morbid': bnd_r['OMIM_morbid'], 'OMIM_morbid_candidate': bnd_r['OMIM_morbid_candidate'],
                       'Repeat_coord_left': bnd_r['Repeat_coord_left'], 'Repeat_type_left': bnd_r['Repeat_type_left'],
                       'Repeat_coord_right': bnd_r['Repeat_coord_right'], 'Repeat_type_right': bnd_r['Repeat_type_right'],
                       'SegDup_left': bnd_r['SegDup_left'], 'SegDup_right': bnd_r['SegDup_right'],
                       'LOEUF_bin': bnd_r['LOEUF_bin'], 'GnomAD_pLI': bnd_r['GnomAD_pLI'],
                       'ExAC_pLI': bnd_r['ExAC_pLI'],
                       'Gnomad_AF': bnd_r['Gnomad_AF'], 'Best_gnomAD_OFP': bnd_r['Best_gnomAD_OFP'],
                       'HPRC_AF': bnd_r['HPRC_AF'], 'Best_HPRC_OFP': bnd_r['Best_HPRC_OFP'],
                       'CoLoRSdb_AF': bnd_r['CoLoRSdb_AF'], 'Best_CoLoRSdb_OFP': bnd_r['Best_CoLoRSdb_OFP'],
                       'CdLS_AF': bnd_r['CdLS_AF'], 'Best_CdLS_OFP': bnd_r['Best_CdLS_OFP'],
                       'annotsv_location_tag': bnd_r['annotsv_location_tag'],
                       'location_boundary_detail': bnd_r['location_boundary_detail'],
                       'cds_utr_context': bnd_r['cds_utr_context'], 'nearest_junction': bnd_r['nearest_junction'],
                       'nearest_junction_dist': bnd_r['nearest_junction_dist']}
            # print(outdict)
            outrows.append(outdict)
            matched_bnd_pairs.add(vid)
            matched_ins_ids.add(ins_id)
            total_rows += 1
            primary_rows += 1

        for (ins_rec, dist) in mate_hits:
            ins_id, ins_chr, ins_pos, ins_svlen, ins_filt = ins_rec
            outdict = {}
            outdict = {'sawfish_id': vid, 'sawfish_chr': bnd_r['SV_chrom'], 'sawfish_pos': bnd_r['SV_start'],
                       'REF': bnd_r['REF'], 'ALT': bnd_r['ALT'],
                       'sawfish_mate_id': mate_id, 'sawfish_mate_chr': mate_chr, 'sawfish_mate_pos': mate_pos,
                       'sawfish_filter': bnd_r['FILTER'], 'sawfish_homlen': homlen, 'sawfish_inslen': inslen,
                       'hitsv_id': ins_id, 'hitsv_chr': ins_chr, 'hitsv_pos': ins_pos, 'hitsv_svlen': ins_svlen,
                       'Annotation_mode': bnd_r['Annotation_mode'], 'hitsv_filter': ins_filt, 'distance': dist,
                       'match_end': 'mate', args.proband: bnd_r[args.proband],
                       'Gene_name': '', 'Frameshift': '', 'GenCC_disease': '', 'GenCC_moi': '',
                       'GenCC_classification': '', 'GenCC_pmid': '', 'OMIM_ID': '', 'OMIM_phenotype': '',
                       'OMIM_inheritance': '', 'OMIM_morbid': '', 'OMIM_morbid_candidate': '',
                       'Repeat_coord_left': '', 'Repeat_type_left': '',
                       'Repeat_coord_right': '', 'Repeat_type_right': '',
                       'SegDup_left': '', 'SegDup_right': '',
                       'LOEUF_bin': '', 'GnomAD_pLI': '',
                       'ExAC_pLI': '',
                       'Gnomad_AF': '', 'Best_gnomAD_OFP': '',
                       'HPRC_AF': '', 'Best_HPRC_OFP': '',
                       'CoLoRSdb_AF': '', 'Best_CoLoRSdb_OFP': '',
                       'CdLS_AF': '', 'Best_CdLS_OFP': '',
                       'annotsv_location_tag': '',
                       'location_boundary_detail': '',
                       'cds_utr_context': '', 'nearest_junction': '',
                       'nearest_junction_dist': ''
                       }
            outrows.append(outdict)
            matched_bnd_pairs.add(vid)
            matched_ins_ids.add(ins_id)
            total_rows += 1
            mate_rows += 1

    out_df = pd.DataFrame(outrows)
    out_df.replace(['nan', 'nan, ', ', nan', 'NaN', 'None'], '', inplace=True)
    out_df = out_df.sort_values(by=['hitsv_id'], ascending=False)
    out_df = out_df.groupby('hitsv_id').agg(
        lambda x: ', '.join([val for val in x.dropna().unique().astype(str) if val.strip() not in ('', 'nan')])).reset_index()

    # out_df = out_df.sort_values(by=['GenCC_disease', 'hitsv_svlen'], ascending=False)
    out_df.to_csv(args.bnd_output, sep='\t', index=False)
    out_df.rename(columns={"hitsv_chr": "SV_chrom", "hitsv_pos":"SV_start", "hitsv_svlen":"LEN", "sawfish_id":"ID"}, inplace=True)
    out_cols = ['hitsv_id', 'sawfish_chr', 'sawfish_pos', 'SV_chrom', 'SV_start', 'LEN', 'REF', 'ALT', args.proband, 'Annotation_mode', 'Gene_name', 'Frameshift', 'GenCC_disease',
                     'GenCC_moi', 'GenCC_classification', 'GenCC_pmid', 'OMIM_ID', 'OMIM_phenotype', 'OMIM_inheritance',
                     'OMIM_morbid', 'OMIM_morbid_candidate', 'Repeat_coord_left', 'Repeat_type_left', 'Repeat_coord_right',
                     'Repeat_type_right', 'SegDup_left', 'SegDup_right', 'LOEUF_bin', 'GnomAD_pLI', 'ExAC_pLI',
                     'Gnomad_AF', 'Best_gnomAD_OFP', 'HPRC_AF', 'Best_HPRC_OFP', 'CoLoRSdb_AF', 'Best_CoLoRSdb_OFP',
                     'CdLS_AF', 'Best_CdLS_OFP', 'annotsv_location_tag', 'location_boundary_detail',
                     'cds_utr_context', 'nearest_junction', 'nearest_junction_dist', 'ID', 'match_end']
    out_df = out_df[out_cols]
    out_df['TYPE'] = "from_BND"

    otherinsertions = sawfishdf[sawfishdf['SVTYPE'] == 'INS']
    print(len(otherinsertions))
    otherinsertions['SV_chrom'] = 'chr'+otherinsertions['SV_chrom'].astype(str)
    ins_df = pd.concat([out_df, otherinsertions], ignore_index=True, join="outer")
    # ins_df.drop_duplicates(subset=['SV_chrom', 'SV_start', 'ID'], inplace=True)
    ins_df = ins_df[~ins_df[args.proband].str.contains('TR_OVERLAP', na=False)]
    ins_df = ins_df[['SV_chrom', 'SV_start', 'LEN','REF', 'ALT','Gene_name', 'Frameshift',args.proband, 'Annotation_mode',
                     'GenCC_disease',
                     'GenCC_moi', 'GenCC_classification', 'GenCC_pmid', 'OMIM_ID', 'OMIM_phenotype', 'OMIM_inheritance',
                     'OMIM_morbid', 'OMIM_morbid_candidate',  'Tx', 'Overlapped_CDS_length', 'Overlapped_CDS_percent',
                     'Repeat_coord_left', 'Repeat_type_left', 'Repeat_coord_right',
                     'Repeat_type_right', 'SegDup_left', 'SegDup_right', 'LOEUF_bin', 'GnomAD_pLI', 'ExAC_pLI',
                     'Gnomad_AF', 'Best_gnomAD_OFP', 'HPRC_AF', 'Best_HPRC_OFP', 'CoLoRSdb_AF', 'Best_CoLoRSdb_OFP',
                     'CdLS_AF', 'Best_CdLS_OFP', 'annotsv_location_tag', 'location_boundary_detail',
                     'cds_utr_context', 'nearest_junction', 'nearest_junction_dist', 'ID', 'TYPE', 'hitsv_id', 'match_end']]
    pd.set_option('display.max_columns', None)
    ins_df = collapse_rows(ins_df, ['SV_chrom', 'SV_start', 'LEN', 'REF', 'ALT'])
    ins_df = ins_df.replace(to_replace=['nan', 'None', '<NA>'], value='', regex=True)
    ins_df = ins_df.map(strip_comma_artifacts)
    # add tag for insertions in coding regions in disease-associated genes
    coding_context = ['CDS', 'mixed', 'UTR']
    conditions = [
        ((ins_df["GenCC_disease"] != '') | (ins_df["OMIM_morbid"] == "yes")) & (ins_df['cds_utr_context'].isin(coding_context))
    ]
    choices = ["1"]
    # Syntax: np.select(conditions_list, choices_list, default_value)
    ins_df['disease-gene-cds'] = np.select(conditions, choices, default="")
    # ins_df['nearest_junction_dist_abs'] = ins_df['nearest_junction_dist'].replace('', np.nan).astype(float).abs()
    # ins_df.sort_values(by='OMIM_morbid', inplace=True)
    ins_df = ins_df.sort_values(
        by=['GenCC_disease', 'OMIM_morbid', 'nearest_junction_dist'],
        key=lambda col: col.replace('', np.nan),
        na_position='last'
    )
    # ins_df = ins_df.sort_values(by=['GenCC_disease', 'OMIM_morbid', 'nearest_junction_dist_abs'], ascending=[False, False, True])
    ins_df.to_csv(args.ins_output, sep='\t', index=False)



    pct_bnd = 100 * len(matched_bnd_pairs) / total_bnd_pairs if total_bnd_pairs else 0
    pct_ins = 100 * len(matched_ins_ids) / total_ins if total_ins else 0

    print(file=sys.stderr)
    print("=== BND-INS Matching Summary ===", file=sys.stderr)
    print(f"Sawfish BND pairs:       {total_bnd_pairs:,}", file=sys.stderr)
    print(f"HitSV INS records:       {total_ins:,}", file=sys.stderr)
    print(f"Distance threshold:      {args.threshold} bp", file=sys.stderr)
    print("--------------------------------", file=sys.stderr)
    print(f"BND pairs matched:       {len(matched_bnd_pairs):,} / {total_bnd_pairs:,}  ({pct_bnd:.1f}%)", file=sys.stderr)
    print(f"HitSV INS matched:       {len(matched_ins_ids):,} / {total_ins:,}  ({pct_ins:.1f}%)", file=sys.stderr)
    print(f"Total output rows:       {total_rows:,}", file=sys.stderr)
    print(f"  Primary-end matches:   {primary_rows:,}", file=sys.stderr)
    print(f"  Mate-end matches:      {mate_rows:,}", file=sys.stderr)
    print(f"Output written to:       {args.bnd_output}", file=sys.stderr)

    unmatched_bnds = sawfishbnd[~sawfishbnd['ID'].isin(matched_bnd_pairs)]
    unmatched_bnds.to_csv('Unmatched_Breakends.tsv', sep='\t', index=False)


def filter_cnvs(args, sawfishdf):
    CNVlist = ['DUP', 'DEL']
    cnvdf = sawfishdf[sawfishdf['SVTYPE'].isin(CNVlist) & (sawfishdf['Annotation_mode'] == 'split')]
    # filter out CNVs that don't affect protein coding transcripts and don't overlap any CDS
    # print(df.shape)
    # print(sawfishdf[sawfishdf['Gene_name'] == 'BCL11A'])

    cnvdf = tag_segdups(cnvdf, args.segdupdb, args.proband)
    # print(list(cnvdf.columns))
    cnvdf = cnvdf[
        ['chrom', 'start', 'end', 'SV_length', 'SVTYPE', 'LEN', 'REF', 'ALT', 'Gene_name',
         'CytoBand', 'annotsv_location_tag', 'location_boundary_detail', 'cds_utr_context', 'nearest_junction',
         'nearest_junction_dist', 'segdup_overlap_bases', 'segdup_overlap_pct', 'Annotation_mode', 'Tx',
         'Overlapped_CDS_length', 'Overlapped_CDS_percent', 'Frameshift', 'Exon_count', 'SegDup_left', 'SegDup_right',
         'GenCC_disease', 'GenCC_moi', 'GenCC_classification', 'GenCC_pmid', 'OMIM_ID',
         'OMIM_phenotype', 'OMIM_inheritance', 'OMIM_morbid', 'Gnomad_AF', 'Best_gnomAD_OFP',
         'HPRC_AF', 'Best_HPRC_OFP', 'CoLoRSdb_AF', 'Best_CoLoRSdb_OFP', 'CdLS_AF', 'Best_CdLS_OFP',
         'ID', 'INFO', 'FORMAT', args.proband]]

    cnvdf['Overlapped_CDS_percent'] = pd.to_numeric(cnvdf['Overlapped_CDS_percent'], errors='coerce')
    # print(cnvdf[cnvdf['ID'] == 'sawfish:CNV:0:0:204000'])
    filteredcnv = cnvdf[(cnvdf['Tx'].str.startswith('NM')) & (cnvdf['Overlapped_CDS_percent'] > 0)]
    tempdf = cnvdf.merge(filteredcnv, how='left', indicator=True)
    filteredcnv = collapse_rows(filteredcnv, ['chrom', 'start', 'end'])
    filteredcnv.to_csv(args.cnv_output, sep='\t', index=False)
    # noncoding_df = cnvdf[((cnvdf['Tx'].str.startswith('NM')) & (cnvdf['Overlapped_CDS_percent'] > 0))]
    noncoding_df = tempdf[tempdf['_merge'] == 'left_only'].drop(columns='_merge')
    noncoding_df = collapse_rows(noncoding_df, ['chrom', 'start', 'end'])
    noncoding_df.to_csv(args.nc_cnv, sep='\t', index=False)


def tag_segdups(out_df, segdupdb, sample):
    sdcolnames = ['bin', 'chrom', 'start', 'end', 'name', 'score', 'strand', 'otherChrom', 'otherStart',
                  'otherEnd', 'otherSize', 'uid', 'posBasesHit', 'testResult', 'verdict', 'chits', 'ccov', 'alignfile',
                  'alignL', 'indelN', 'indelS', 'alignB', 'matchB', 'mismatchB', 'transitionsB', 'transversionsB',
                  'fracMatch', 'fracMatchIndel', 'jcK', 'k2K']
    segdupdf = pd.read_csv(segdupdb, sep='\t', header=None, names=sdcolnames)
    segdupdf['segdup_coordinates'] = segdupdf['chrom'] + ':' + segdupdf['start'].astype(str) + '-' + segdupdf[
        'end'].astype(str)
    segdupdf['segdup_OtherCoordinates'] = segdupdf['otherChrom'] + ':' + segdupdf['otherStart'].astype(str) + '-' + \
                                          segdupdf['otherEnd'].astype(str)
    segdupdf = segdupdf[
        ['chrom', 'start', 'end', 'segdup_coordinates', 'segdup_OtherCoordinates', 'fracMatch', 'fracMatchIndel']]
    out_df.rename(columns={'SV_chrom': 'chrom', 'SV_start': 'start', 'SV_end': 'end'}, inplace=True)
    out_df['chrom'] = 'chr' + out_df['chrom'].astype(str)
    overlapping_intervals = bf.overlap(out_df, segdupdf, suffixes=('', '_segdup'))
    # overlapping_intervals.to_csv('temp.tsv', sep='\t', index=False)
    overlap_bases = (overlapping_intervals.groupby(['chrom', 'start', 'end'])
                      .apply(total_overlap_bases, include_groups=False)
                      .rename('segdup_overlap_bases').reset_index())

    merged_df = (overlapping_intervals.groupby(['chrom', 'start', 'end', 'Gene_name']).agg(lambda x: list(set(x))).reset_index())
    merged_df = merged_df.drop(columns=['chrom_segdup', 'start_segdup', 'end_segdup'])
    merged_df.fillna('', inplace=True)
    target_cols = [x for x in merged_df.columns if x not in ['chrom', 'start', 'end', 'Gene_name']]
    merged_df[target_cols] = merged_df[target_cols].map(lambda lst: ", ".join(map(str, lst)) if isinstance(lst, list) else lst)
    merged_df = merged_df.merge(overlap_bases, on=['chrom', 'start', 'end'], how='left')
    merged_df['segdup_overlap_bases'] = merged_df['segdup_overlap_bases'].fillna(0).astype(int)
    merged_df['segdup_overlap_pct'] = (
        100 * merged_df['segdup_overlap_bases'] / (merged_df['end'] - merged_df['start'])
    ).round(2)
    merged_df = merged_df.replace(to_replace=['nan', 'None', '<NA>', 'nan, ', ', nan'], value='', regex=True)
    for col in target_cols:
        merged_df[col] = merged_df[col].str.findall(r'[^\s,](?:[^,]*[^\s,])?').str.join(', ')
    return merged_df

def strip_comma_artifacts(val):
    """Remove leading/trailing (and empty in-between) comma-separated tokens, e.g. 'no, ' -> 'no', ', no' -> 'no'."""
    if not isinstance(val, str):
        return val
    tokens = [tok.strip() for tok in val.split(',')]
    return ', '.join(tok for tok in tokens if tok != '')


def collapse_rows(df, cols):
    df = (
        df.groupby(cols).agg(lambda x: sorted({
            str(v) for v in x if pd.notna(v) and str(v).strip() != ''
        })).reset_index())
    target_cols = [x for x in df.columns if x not in cols]
    df[target_cols] = df[target_cols].map(
        lambda lst: ", ".join(lst) if isinstance(lst, list) else lst)
    for col in target_cols:
        # print(col, type(list(df[col])[0]))
        df[col] = df[col].str.findall(r'[^\s,](?:[^,]*[^\s,])?').str.join(', ')
    return df

def total_overlap_bases(group):
    # merge overlapping segdup segments first so bases covered by more than
    # one segdup entry aren't counted twice, then clip to the query interval
    _, query_start, query_end = group.name
    segments = sorted(
        (max(query_start, seg_start), min(query_end, seg_end))
        for seg_start, seg_end in zip(group['start_segdup'], group['end_segdup'])
        if pd.notna(seg_start) and pd.notna(seg_end) and min(query_end, seg_end) > max(query_start, seg_start)
    )
    if not segments:
        return 0
    merged_segments = [segments[0]]
    for seg_start, seg_end in segments[1:]:
        last_start, last_end = merged_segments[-1]
        if seg_start <= last_end:
            merged_segments[-1] = (last_start, max(last_end, seg_end))
        else:
            merged_segments.append((seg_start, seg_end))
    return sum(seg_end - seg_start for seg_start, seg_end in merged_segments)

def main():
    parser = argparse.ArgumentParser(
        description='Match Sawfish BND breakends against HitSV INS insertions'
    )
    parser.add_argument('--sawfish', required=True, help='Filtered TSV (.tsv)')
    parser.add_argument('--hitsv', required=True, help='HitSV VCF (.vcf or .vcf.gz)')
    parser.add_argument('--segdupdb', required=True, help='Reference file containing segdup coordinates')
    parser.add_argument('--proband', required=True, help='proband ID')
    parser.add_argument('--bnd_output', default='bnd_ins_matches.tsv', help='Output TSV path for BND matches')
    parser.add_argument('--ins_output', required=True, help='Output TSV path for Insertions')
    parser.add_argument('--cnv_output', required=True, help='Output TSV path for indels and CNVs')
    parser.add_argument('--nc_cnv', required=True, help='Output TSV path for indels and CNVs')
    parser.add_argument('--threshold', type=int, default=50,
                        help='Max distance (bp) to consider a match (default: 50)')
    args = parser.parse_args()

    sawfishdf = pd.read_csv(args.sawfish, sep='\t', header=0)
    warnings.filterwarnings("ignore")
    pd.set_option('display.max_columns', None)
    run_matching(args, sawfishdf)
    filter_cnvs(args, sawfishdf)


if __name__ == '__main__':
    main()