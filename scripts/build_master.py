#!/usr/bin/env python3
"""Add TSS and TES point rows to a lifted GTF to make the master annotation."""
from __future__ import annotations

import argparse

from gtflib import format_attrs, parse_attrs, read_fasta, read_gtf, write_gtf

PREFIX = "KT899744delRR_"
# The collinear alignment block stops 3 bp short of the genome end (135161 of
# 135164), so a clipped transcript's end lands there, not on the last base.
EDGE_SLACK = 3
FEATURE_RANK = {f: i for i, f in enumerate(
    ["gene", "transcript", "exon", "CDS", "start_codon", "stop_codon", "TSS", "TES"])}


def tss_tes_rows(rows: list[list[str]], genome_len: int | None = None) -> list[list[str]]:
    """One TSS and one TES per transcript with exons, collapsed by position and strand.

    An end is skipped when the transcript was clipped at a genome boundary
    (liftover_status contains "clipped") and that end sits at position 1 or
    within EDGE_SLACK bases of genome_len, because the real transcript end lies outside the genome.
    """
    clipped: set[str] = set()
    ends: dict[str, list] = {}
    for r in rows:
        if r[2] == "transcript":  # liftover_status is carried on transcript rows
            a = parse_attrs(r[8])
            if "clipped" in a.get("liftover_status", "") and "transcript_id" in a:
                clipped.add(a["transcript_id"])
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
            if tid in clipped and (pos == 1 or (genome_len is not None and pos >= genome_len - EDGE_SLACK)):
                continue
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
    return (int(r[3]), -int(r[4]), FEATURE_RANK.get(r[2], len(FEATURE_RANK)), r[6], r[8])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-gtf", required=True)
    ap.add_argument("--out-gtf", required=True)
    ap.add_argument("--version", required=True)
    ap.add_argument("--fasta", required=True, help="genome FASTA (for its length)")
    ap.add_argument("--mirna-gtf", help="optional miRNA GTF, appended before sorting")
    a = ap.parse_args()
    rows = read_gtf(a.in_gtf)
    extra = read_gtf(a.mirna_gtf) if a.mirna_gtf else []
    rows = sorted(rows + tss_tes_rows(rows, len(read_fasta(a.fasta))) + extra, key=sort_key)
    write_gtf(a.out_gtf, rows, [f"#!annotation-version {a.version}"])


if __name__ == "__main__":
    main()
