#!/usr/bin/env python3
"""Add TSS and TES point rows to a lifted GTF to make the master annotation."""
from __future__ import annotations

import argparse

from gtflib import format_attrs, parse_attrs, read_gtf, write_gtf

PREFIX = "KT899744delRR_"
FEATURE_RANK = {f: i for i, f in enumerate(
    ["gene", "transcript", "exon", "CDS", "start_codon", "stop_codon", "TSS", "TES"])}


def tss_tes_rows(rows: list[list[str]]) -> list[list[str]]:
    """One TSS and one TES per transcript with exons, collapsed by position and strand."""
    ends: dict[str, list] = {}
    for r in rows:
        if r[2] != "exon":
            continue
        tid = parse_attrs(r[8]).get("transcript_id")
        if tid is None:
            continue
        e = ends.setdefault(tid, [r[0], r[6], int(r[3]), int(r[4])])
        e[2] = min(e[2], int(r[3]))
        e[3] = max(e[3], int(r[4]))
    points: dict[tuple, set] = {}
    seqname: dict[tuple, str] = {}
    for tid, (seq, strand, lo, hi) in ends.items():
        five, three = (lo, hi) if strand == "+" else (hi, lo)
        for feat, pos in (("TSS", five), ("TES", three)):
            key = (feat, pos, strand)
            points.setdefault(key, set()).add(tid)
            seqname.setdefault(key, seq)
    out = []
    for (feat, pos, strand), tids in points.items():
        name = f"{PREFIX}{feat}_{pos}{'p' if strand == '+' else 'm'}"
        attrs = {"gene_id": name, "transcript_id": name, "gene_name": name,
                 "source_transcripts": ",".join(sorted(tids))}
        out.append([seqname[(feat, pos, strand)], "Genbank", feat, str(pos), str(pos),
                    ".", strand, ".", format_attrs(attrs)])
    return out


def sort_key(r: list[str]) -> tuple:
    return (int(r[3]), int(r[4]), FEATURE_RANK.get(r[2], len(FEATURE_RANK)), r[6], r[8])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-gtf", required=True)
    ap.add_argument("--out-gtf", required=True)
    ap.add_argument("--version", required=True)
    a = ap.parse_args()
    rows = read_gtf(a.in_gtf)
    rows = sorted(rows + tss_tes_rows(rows), key=sort_key)
    write_gtf(a.out_gtf, rows, [f"#!annotation-version {a.version}"])


if __name__ == "__main__":
    main()
