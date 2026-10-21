#!/usr/bin/env python3
"""
tag_annotsv.py

Adds location/junction annotation columns to an AnnotSV TSV. All columns are
populated for split rows only; full (summary) rows get "" / NA.

Columns added:
  annotsv_location_tag       coarse-but-refined categorical tag from Location
  location_boundary_detail   for 5prime/3prime: internal endpoint exon vs intron
  cds_utr_context            CDS / UTR / mixed, from AnnotSV's Location2
  nearest_junction           named nearest intron-exon junction (from MANE GTF)
  nearest_junction_dist      signed bp distance to that junction (transcription
                             orientation: negative = breakpoint 5' of junction)

annotsv_location_tag values (derived from AnnotSV's Location; strand-agnostic
because Location is always written left-to-right in transcription orientation):
  whole_gene      txStart-txEnd
  5prime          starts txStart- (not whole_gene)
  3prime          ends -txEnd (not whole_gene)
  single_exon     exonN-exonN (same N)
  multi_exon      exonM-exonN (M != N)
  single_intron   intronN-intronN (same N)
  multi_intron    intronM-intronN (M != N)
  exon_intron     one endpoint exon, the other intron
  intragenic_other  defensive catch-all (not expected)
  "" (empty)      full rows

Nearest-junction naming uses the MANE RefSeq GTF the pipeline already feeds to
AnnotSV for transcript selection. AnnotSV's Tx (RefSeq NM_/NR_/XM_...) joins
directly to the GTF transcript_id (minus its version). Introns are not stored
in the GTF: intron k is the gap between exon k and exon k+1 (transcription
order), bounded by the exonK/intronK donor and the intronK/exon(K+1) acceptor.
A breakpoint is anchored to whichever real splice junction is closest.
"""

import argparse
import gzip
import re
import sys

import pandas as pd


# ---- location tag (from AnnotSV Location string) ---------------------------

_EXON_RE = re.compile(r"^exon(\d+)$")
_INTRON_RE = re.compile(r"^intron(\d+)$")


def refined_location_tag(loc):
    """Refined categorical tag from an AnnotSV Location string (split rows)."""
    loc = str(loc).strip()
    if not loc:
        return ""
    if loc == "txStart-txEnd":
        return "whole_gene"
    if loc.startswith("txStart-"):
        return "5prime"
    if loc.endswith("-txEnd"):
        return "3prime"

    if "-" not in loc:
        return "intragenic_other"
    left, right = loc.split("-", 1)
    le, ri = _EXON_RE.match(left), _EXON_RE.match(right)
    li, rin = _INTRON_RE.match(left), _INTRON_RE.match(right)

    if le and ri:
        return "single_exon" if le.group(1) == ri.group(1) else "multi_exon"
    if li and rin:
        return "single_intron" if li.group(1) == rin.group(1) else "multi_intron"
    if (le and rin) or (li and ri):
        return "exon_intron"
    return "intragenic_other"


def boundary_detail(tag, loc):
    """For 5prime/3prime: is the internal (non-tx-boundary) endpoint exon or intron?"""
    loc = str(loc).strip()
    if tag == "5prime":
        internal = loc.split("-", 1)[1] if "-" in loc else ""
    elif tag == "3prime":
        internal = loc.rsplit("-", 1)[0] if "-" in loc else ""
    else:
        return ""
    if _EXON_RE.match(internal):
        return f"{tag}_exon"
    if _INTRON_RE.match(internal):
        return f"{tag}_intron"
    return ""


def cds_utr_context(loc2):
    """Coding-region context from AnnotSV's Location2 (split rows).

    Per the AnnotSV manual, Location2 tokens are exactly one of:
      CDS     between CDS start and end (may be exonic or intronic)
      5'UTR   before the CDS start (coding gene)
      3'UTR   after the CDS end (coding gene)
      UTR     the gene has NO CDS at all (non-coding transcript, e.g. NR_)
    and Location2 may be a range of two tokens (e.g. 5'UTR-CDS, CDS-3'UTR,
    5'UTR-3'UTR). We split on '-' and inspect the token set (5'UTR/3'UTR
    contain the substring 'UTR', so match exact tokens, not substrings).

    Returns:
      CDS         touches only coding sequence
      UTR         touches only 5'/3' UTR of a coding gene
      mixed       spans CDS and UTR, or spans both 5' and 3' UTR
      non_coding  gene has no CDS (bare 'UTR' token)
      ""          empty / full rows
    """
    loc2 = str(loc2).strip()
    if not loc2:
        return ""
    tokens = set(loc2.split("-"))
    has_cds = "CDS" in tokens
    has_5utr = "5'UTR" in tokens
    has_3utr = "3'UTR" in tokens
    has_bare_utr = "UTR" in tokens  # non-coding gene (no CDS)

    if has_cds and (has_5utr or has_3utr):
        return "mixed"
    if has_5utr and has_3utr:  # spans the whole coding region 5'->3'
        return "mixed"
    if has_cds:
        return "CDS"
    if has_5utr or has_3utr:
        return "UTR"
    if has_bare_utr:
        return "non_coding"
    return ""


# ---- MANE GTF exon model + nearest-junction naming -------------------------

_TXID_RE = re.compile(r'transcript_id "([^"]+)"')
_EXNO_RE = re.compile(r'exon_number "?(\d+)"?')


def load_mane_exons(gtf_path):
    """Parse the MANE RefSeq GTF into {transcript_id_no_version: model}.

    model = {"strand": "+"/"-", "exons": [(exon_number, start, end), ...]}
    Coordinates are the GTF's native 1-based inclusive.
    """
    opener = gzip.open if gtf_path.endswith(".gz") else open
    tx = {}
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            attrs = cols[8]
            m_id = _TXID_RE.search(attrs)
            m_no = _EXNO_RE.search(attrs)
            if not m_id or not m_no:
                continue
            txid = m_id.group(1).split(".", 1)[0]
            rec = tx.get(txid)
            if rec is None:
                rec = {"strand": cols[6], "exons": []}
                tx[txid] = rec
            rec["exons"].append((int(m_no.group(1)), int(cols[3]), int(cols[4])))
    for rec in tx.values():
        rec["exons"].sort(key=lambda e: e[0])  # transcription order (exon_number)
    return tx


def _junctions(model):
    """Yield (genomic_pos, name) for every splice junction of a transcript.

    Donor = 3' end of an exon (transcription); acceptor = 5' start of the next.
    On + strand transcription-3' is the larger genomic coord; on - strand it is
    the smaller. Introns are the gaps between consecutive (by exon_number) exons.
    """
    exons = model["exons"]
    plus = model["strand"] == "+"
    for i in range(len(exons) - 1):
        k_no, k_start, k_end = exons[i]
        n_no, n_start, n_end = exons[i + 1]
        if plus:
            donor_pos, acc_pos = k_end, n_start
        else:
            donor_pos, acc_pos = k_start, n_end
        yield donor_pos, f"exon{k_no}/intron{k_no} donor"
        yield acc_pos, f"intron{k_no}/exon{n_no} acceptor"


def nearest_junction(txid, sv_start, sv_end, mane):
    """Return (name, signed_dist) of the nearest splice junction to the SV.

    Considers both breakpoints; picks the globally nearest junction. Distance is
    signed in transcription orientation: negative = breakpoint lies 5' (upstream)
    of the junction, positive = 3' (downstream).
    """
    model = mane.get(txid)
    if model is None:
        return "NA_no_mane_tx", "NA"
    if len(model["exons"]) < 2:
        return "NA_single_exon_tx", "NA"

    plus = model["strand"] == "+"
    # SV_start is 0-based; convert the left breakpoint to 1-based inclusive.
    breakpoints = [sv_start + 1, sv_end]

    best = None  # (abs_dist, signed_dist, name)
    for jpos, jname in _junctions(model):
        for bp in breakpoints:
            adist = abs(bp - jpos)
            signed = (bp - jpos) if plus else (jpos - bp)
            if best is None or adist < best[0]:
                best = (adist, signed, jname)
    return best[2], best[1]


# ---- driver ----------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Add location/junction tags to AnnotSV TSV.")
    ap.add_argument("--input", required=True, help="AnnotSV TSV input")
    ap.add_argument("--mane-gtf", required=True, help="MANE RefSeq GTF (gzipped ok)")
    ap.add_argument("--output", required=True, help="Tagged TSV output")
    args = ap.parse_args()

    print(f"Loading MANE GTF {args.mane_gtf}", file=sys.stderr)
    mane = load_mane_exons(args.mane_gtf)
    print(f"Indexed {len(mane)} transcripts", file=sys.stderr)

    print(f"Loading {args.input}", file=sys.stderr)
    df = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)

    is_split = df["Annotation_mode"] == "split"

    df["annotsv_location_tag"] = ""
    df["location_boundary_detail"] = ""
    df["cds_utr_context"] = ""
    df["nearest_junction"] = ""
    df["nearest_junction_dist"] = ""

    resolved = fallback = 0
    for idx in df.index[is_split]:
        loc = df.at[idx, "Location"]
        tag = refined_location_tag(loc)
        df.at[idx, "annotsv_location_tag"] = tag
        df.at[idx, "location_boundary_detail"] = boundary_detail(tag, loc)
        df.at[idx, "cds_utr_context"] = cds_utr_context(df.at[idx, "Location2"])

        txid = str(df.at[idx, "Tx"]).strip()
        try:
            sv_start = int(df.at[idx, "SV_start"])
            sv_end = int(df.at[idx, "SV_end"])
        except ValueError:
            sv_start = sv_end = None
        if txid and sv_start is not None:
            name, dist = nearest_junction(txid, sv_start, sv_end, mane)
        else:
            name, dist = "NA_no_mane_tx", "NA"
        df.at[idx, "nearest_junction"] = name
        df.at[idx, "nearest_junction_dist"] = str(dist)
        if name.startswith("NA_"):
            fallback += 1
        else:
            resolved += 1

    print(f"Split rows: {resolved} junction-resolved, {fallback} fallback/NA",
          file=sys.stderr)

    df.to_csv(args.output, sep="\t", index=False)
    print(f"Written: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
