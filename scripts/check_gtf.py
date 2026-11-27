#!/usr/bin/env python3
"""Structural check of a lifted GTF against its genome FASTA.

Exit 0 on pass; otherwise prints one error per line and exits 1.
TSS/TES rows are standalone point features: bounds, gene_id, transcript_id and
no-space checks only.
"""
from __future__ import annotations

import argparse
import collections
import re
import sys

from gtflib import read_fasta, read_gtf, revcomp

_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {x + y + z: _AA[i * 16 + j * 4 + k]
         for i, x in enumerate(_B) for j, y in enumerate(_B) for k, z in enumerate(_B)}
POINT = ("TSS", "TES")
_PAIR = re.compile(r'(\S+) "([^"]*)";')


def _seq(genome: str, segs: list[tuple[int, int]], strand: str) -> str:
    s = "".join(genome[a - 1:e] for a, e in sorted(segs))
    return revcomp(s) if strand == "-" else s


def check(rows: list[list[str]], genome: str) -> list[str]:
    """Return a list of error strings (empty when the GTF passes)."""
    errs: list[str] = []
    T: dict = collections.defaultdict(lambda: collections.defaultdict(list))
    info: dict = {}
    genes: dict = {}
    gene_biotype: dict = {}
    gene_cds: dict = collections.defaultdict(int)
    parsed = []
    for f in rows:
        if len(f) != 9:
            errs.append("row does not have 9 columns: %s" % "\t".join(f))
            continue
        s, e = int(f[3]), int(f[4])
        label = "%s %s-%s" % (f[2], f[3], f[4])
        if not 1 <= s <= e <= len(genome):
            errs.append("out of bounds: %s" % label)
            continue
        at = dict(_PAIR.findall(f[8]))
        for k, v in at.items():
            if " " in v:
                errs.append("space in attribute value %s=%r: %s" % (k, v, label))
        if "gene_id" not in at:
            errs.append("missing gene_id: %s" % label)
            continue
        if f[2] != "gene" and "transcript_id" not in at:
            errs.append("missing transcript_id: %s (%s)" % (label, at["gene_id"]))
            continue
        parsed.append((f, s, e, at))
        if f[2] == "gene":
            genes[at["gene_id"]] = (s, e, f[6])
            gene_biotype[at["gene_id"]] = at.get("gene_biotype", "")
    for f, s, e, at in parsed:
        if f[2] == "gene" or f[2] in POINT:
            continue
        t = at["transcript_id"]
        T[t][f[2]].append((s, e, f[7]))
        info[t] = (f[6], at)
        if at["gene_id"] not in genes:
            errs.append("no gene line for %s" % at["gene_id"])
            continue
        gs, ge, gst = genes[at["gene_id"]]
        if not (gs <= s and e <= ge and gst == f[6]):
            errs.append("outside gene or strand mismatch: %s %s-%s (%s)" % (f[2], s, e, t))
    for t, d in T.items():
        st, at = info[t]
        if len(d["transcript"]) != 1:
            errs.append("transcript line count != 1: %s" % t)
            continue
        ex = sorted((s, e) for s, e, _ in d["exon"])
        if not ex:
            errs.append("no exons: %s" % t)
            continue
        if not all(ex[i][1] < ex[i + 1][0] for i in range(len(ex) - 1)):
            errs.append("overlapping exons: %s" % t)
        if d["transcript"][0][:2] != (ex[0][0], ex[-1][1]):
            errs.append("transcript span != exon span: %s" % t)
        if not d["CDS"]:
            continue
        gene_cds[at["gene_id"]] += 1
        for s, e, _ in d["CDS"] + d["stop_codon"] + d["start_codon"]:
            if not any(a <= s and e <= b2 for a, b2 in ex):
                errs.append("CDS outside exon: %s" % t)
        order = sorted(d["CDS"], reverse=(st == "-"))
        done = None
        for s, e, fr in order:  # frame consistency
            if done is not None and int(fr) != (3 - done % 3) % 3:
                errs.append("frame: %s" % t)
            done = (e - s + 1 - int(fr)) if done is None else done + e - s + 1
        nt = _seq(genome, [(s, e) for s, e, _ in d["CDS"]], st)[int(order[0][2]):]
        if d["stop_codon"]:
            if _seq(genome, [(s, e) for s, e, _ in d["stop_codon"]], st) not in ("TAA", "TAG", "TGA"):
                errs.append("stop: %s" % t)
            if sum(e - s + 1 for s, e, _ in d["stop_codon"]) != 3:
                errs.append("stop codon length != 3: %s" % t)
    for gid in genes:
        if gene_cds[gid] == 0 and gene_biotype[gid] not in ("noncoding", "miRNA"):
            errs.append("coding gene without CDS: %s" % gid)
    return errs


def summary(rows: list[list[str]], genome: str) -> str:
    T = {(f[8].split('transcript_id "')[1].split('"')[0]) for f in rows
         if f[2] == "transcript" and 'transcript_id "' in f[8]}
    g = sum(1 for f in rows if f[2] == "gene")
    return "transcripts %d genes %d" % (len(T), g)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--fasta", required=True)
    a = ap.parse_args()
    genome = read_fasta(a.fasta)
    rows = read_gtf(a.gtf)
    errs = check(rows, genome)
    if errs:
        for e in errs:
            print(e, file=sys.stderr)
        sys.exit(1)
    print("PASS " + summary(rows, genome))


if __name__ == "__main__":
    main()
