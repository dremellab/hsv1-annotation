#!/usr/bin/env python3
"""Strand-aware ORF containment classification from the master GTF CDS rows."""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass

from gtflib import parse_attrs, read_gtf

PREFIX = "KT899744delRR_"
Footprint = tuple[str, frozenset[int]]


@dataclass
class OrfClass:
    gene_id: str
    cls: str  # non_overlapping | partially_overlapping | embedded | retained_zero_unique
    host: str = ""
    reason: str = ""  # "" | embedded_in_single_orf | identical_footprint
    length: int = 0
    unique_len: int = 0


def footprints(rows: list[list[str]]) -> dict[str, Footprint]:
    """gene_id -> (strand, union of CDS bases across all transcripts of the gene)."""
    bases: dict[str, set[int]] = {}
    strand: dict[str, str] = {}
    for r in rows:
        if r[2] != "CDS":
            continue
        gid = parse_attrs(r[8])["gene_id"]
        strand[gid] = r[6]
        bases.setdefault(gid, set()).update(range(int(r[3]), int(r[4]) + 1))
    return {g: (strand[g], frozenset(b)) for g, b in bases.items()}


def annotated_lengths(rows: list[list[str]]) -> dict[str, int]:
    """gene_id -> sum of lengths of all its CDS rows (used only for identical-footprint ties)."""
    out: dict[str, int] = {}
    for r in rows:
        if r[2] == "CDS":
            gid = parse_attrs(r[8])["gene_id"]
            out[gid] = out.get(gid, 0) + int(r[4]) - int(r[3]) + 1
    return out


def _short(gid: str) -> str:
    return gid[len(PREFIX):] if gid.startswith(PREFIX) else gid


def classify(fp: dict[str, Footprint],
             annotated: dict[str, int] | None = None) -> dict[str, OrfClass]:
    annotated = annotated or {}
    embedded: dict[str, tuple[str, str]] = {}  # gene -> (host, reason)
    for a, (sa, ba) in fp.items():
        # strict containment in a single other same-strand gene
        hosts = [b for b, (sb, bb) in fp.items()
                 if b != a and sb == sa and ba < bb]
        if hosts:
            host = min(hosts, key=lambda b: (len(fp[b][1]), b))
            embedded[a] = (host, "embedded_in_single_orf")
            continue
        twins = sorted(b for b, (sb, bb) in fp.items() if sb == sa and bb == ba)
        if len(twins) > 1:
            def rank(g: str) -> tuple:
                plain = " " not in _short(g) and "_" not in _short(g)
                return (not plain, -annotated.get(g, 0), g)
            keeper = min(twins, key=rank)
            if keeper != a:
                embedded[a] = (keeper, "identical_footprint")

    kept = [g for g in fp if g not in embedded]
    out: dict[str, OrfClass] = {}
    for g, (s, b) in fp.items():
        others = [o for o in (kept if g not in embedded else fp) if o != g and fp[o][0] == s]
        covered: set[int] = set()
        for o in others:
            covered |= fp[o][1]
        unique = len(b - covered)
        if g in embedded:
            host, reason = embedded[g]
            out[g] = OrfClass(g, "embedded", host, reason, len(b), unique)
            continue
        if unique == 0:
            cls = "retained_zero_unique"
        elif any(fp[o][0] == s and b & fp[o][1] for o in fp if o != g):
            # shares a base with any same-strand ORF, kept or excluded
            cls = "partially_overlapping"
        else:
            cls = "non_overlapping"
        out[g] = OrfClass(g, cls, "", "", len(b), unique)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", required=True, help="master GTF")
    ap.add_argument("--out", required=True, help="classification TSV")
    a = ap.parse_args()
    rows = read_gtf(a.gtf)
    res = classify(footprints(rows), annotated_lengths(rows))
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene_id", "class", "host", "reason", "length", "unique_len"])
        for g in sorted(res):
            c = res[g]
            w.writerow([g, c.cls, c.host, c.reason, c.length, c.unique_len])


if __name__ == "__main__":
    main()
