#!/usr/bin/env python3
"""Place miRBase HSV-1 miRNAs on KT899744delRR and write GTF, BED and a report.

Each precursor (hairpin) is placed by the first method that works:
  1. exact_target             exact match (either strand) on KT899744delRR
  2. projected_from_BK012101.1  exact match on BK012101.1, both ends projected
                              through the collinear minimap2 block
     projected_refined       a projection that changed the length is re-aligned
                              (Smith-Waterman, as in 3) near the projected interval
  3. local_alignment          Smith-Waterman of the hairpin around its already
                              placed mature arm(s) (flag "approximate")
  4. unplaced                 report row only; reason "in_deleted_region" when the
                              sequence matches KT899744.1 exactly inside a block
                              removed from KT899744delRR, else "no_match"
Mature arms are placed by methods 1 and 2 and kept only when they lie inside a
placed precursor of the same miRNA on the same strand.
"""
from __future__ import annotations

import argparse
import os
import re
from dataclasses import dataclass, field

from build_master import sort_key
from gtflib import format_attrs, read_fasta, revcomp, write_gtf
from lift_bk012101 import build_map
from verify_delrr import KEEP

SOURCE = "miRBase_23.1"
WINDOW = 150
MATCH, MISMATCH, GAP = 2, -3, -5
EXACT, PROJECTED, LOCAL, UNPLACED = (
    "exact_target", "projected_from_BK012101.1", "local_alignment", "unplaced")
REFINED = "projected_refined"
REPORT_COLS = ["name", "mirbase_id", "type", "chrom", "start", "end", "strand",
               "placement", "identity_pct", "aln_length", "flags"]


@dataclass
class Rec:
    name: str
    acc: str
    typ: str  # precursor | mature
    start: int = 0
    end: int = 0
    strand: str = "."
    method: str = UNPLACED
    identity: float | None = None
    aln_len: int | None = None
    flags: list[str] = field(default_factory=list)
    gene_id: str = ""
    tid: str = ""


def read_named_fasta(path: str) -> list[tuple[str, str, str]]:
    """(name, accession, sequence) per record; accession is the 2nd header field."""
    out: list[tuple[str, str, str]] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                f = line[1:].split()
                out.append((f[0], f[1] if len(f) > 1 else "NA", ""))
            elif line:
                n, a, s = out[-1]
                out[-1] = (n, a, s + line)
    return out


def exact_hits(genome: str, seq: str) -> list[tuple[int, int, str]]:
    """1-based (start, end, strand) of exact matches on either strand."""
    hits = []
    for strand, pat in (("+", seq), ("-", revcomp(seq))):
        i = genome.find(pat)
        while i != -1:
            hits.append((i + 1, i + len(pat), strand))
            i = genome.find(pat, i + 1)
    return sorted(hits)


def hit_flags(hits: list[tuple[int, int, str]]) -> list[str]:
    if len(hits) <= 1:
        return []
    if (len(hits) == 2 and hits[0][:2] == hits[1][:2]
            and {hits[0][2], hits[1][2]} == {"+", "-"}):
        return ["strand_ambiguous"]
    return ["multi_hit"]


def project(hits, pmap):
    """Project source hits through the collinear block; keep hits with both ends mapped."""
    m, _, _ = pmap
    out = []
    for s, e, strand in hits:
        ts, te = m[s], m[e]
        if ts and te and ts <= te:
            out.append((ts, te, strand, e - s))
    return out


def place_hits(seq, tgt, src, pmap):
    """Rules 1 and 2. Returns (hits, method, flags); hits are (start, end, strand)."""
    hits = exact_hits(tgt, seq)
    if hits:
        return hits, EXACT, hit_flags(hits)
    if src and pmap:
        sh = exact_hits(src, seq)
        proj = project(sh, pmap)
        if proj:
            hits = [(s, e, st) for s, e, st, _ in proj]
            flags = hit_flags(hits)
            if any(e - s != d for s, e, _, d in proj):
                flags.append("length_changed")
            return hits, PROJECTED, flags
    return [], UNPLACED, []


def smith_waterman(a: str, b: str):
    """Local alignment of a (query) in b. Returns (score, b_start0, b_end0excl, columns, matches)."""
    n, m = len(a), len(b)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    best = (0, 0, 0)
    for i in range(1, n + 1):
        ai, Hi, Hp = a[i - 1], H[i], H[i - 1]
        for j in range(1, m + 1):
            v = max(0, Hp[j - 1] + (MATCH if ai == b[j - 1] else MISMATCH),
                    Hp[j] + GAP, Hi[j - 1] + GAP)
            Hi[j] = v
            if v > best[0]:
                best = (v, i, j)
    score, i, j = best
    j_end, cols, matches = j, 0, 0
    while i > 0 and j > 0 and H[i][j] > 0:
        h = H[i][j]
        if h == H[i - 1][j - 1] + (MATCH if a[i - 1] == b[j - 1] else MISMATCH):
            matches += a[i - 1] == b[j - 1]
            i, j = i - 1, j - 1
        elif h == H[i - 1][j] + GAP:
            i -= 1
        else:
            j -= 1
        cols += 1
    return score, j, j_end, cols, matches


def align_near(seq, tgt, start, end, strand):
    """Smith-Waterman of seq within WINDOW nt of [start, end] on strand.

    Returns (score, start, end, columns, matches) in 1-based target coordinates, or None.
    """
    lo, hi = max(1, start - WINDOW), min(len(tgt), end + WINDOW)
    win = tgt[lo - 1:hi]
    if strand == "-":
        win = revcomp(win)
    score, b0, b1, cols, matches = smith_waterman(seq, win)
    if cols == 0:
        return None
    if strand == "+":
        return score, lo + b0, lo + b1 - 1, cols, matches
    return score, hi - b1 + 1, hi - b0, cols, matches


def precursor_of(mature_name: str) -> str:
    return re.sub(r"-[35]p$", "", mature_name).replace("-miR-", "-mir-")


def unplaced_reason(seq: str, kt: str, keep) -> str:
    """"in_deleted_region" if seq matches KT899744.1 entirely inside a removed block."""
    if kt:
        deleted, prev = [], 0
        for a, b in sorted(keep):
            if a - 1 > prev:
                deleted.append((prev + 1, a - 1))
            prev = b
        if len(kt) > prev:
            deleted.append((prev + 1, len(kt)))
        for s, e, _ in exact_hits(kt, seq):
            if any(a <= s and e <= b for a, b in deleted):
                return "in_deleted_region"
    return "no_match"


def place(hairpins, matures, tgt, src, pmap, kt="", keep=KEEP):
    """Place all miRNAs; returns a list of Rec (precursors first, then matures)."""
    hairpins = [(n, a, s.upper().replace("U", "T")) for n, a, s in hairpins]
    matures = [(n, a, s.upper().replace("U", "T")) for n, a, s in matures]
    mat_hits = {n: place_hits(s, tgt, src, pmap) for n, _, s in matures}
    pre_recs: dict[str, list[Rec]] = {}
    out: list[Rec] = []
    for name, acc, seq in hairpins:
        hits, method, flags = place_hits(seq, tgt, src, pmap)
        recs = []
        if hits and method == PROJECTED and "length_changed" in flags:
            # the projected interval is unreliable (e.g. a homopolymer indel at an edge)
            refined = [(h, align_near(seq, tgt, h[0], h[1], h[2])) for h in hits]
            if all(a and a[3] >= 0.8 * len(seq) and a[4] >= 0.8 * a[3] for _, a in refined):
                for (s0, e0, st), (_, s, e, cols, matches) in refined:
                    recs.append(Rec(name, acc, "precursor", s, e, st, REFINED,
                                    round(100.0 * matches / cols, 1), cols, list(flags)))
                hits = []
                method = REFINED
        if hits:
            for s, e, st in hits:
                recs.append(Rec(name, acc, "precursor", s, e, st, method,
                                100.0 if method == EXACT else None, e - s + 1, list(flags)))
        elif method != REFINED:
            best = None
            for mn, _, _ in matures:
                if precursor_of(mn) != name:
                    continue
                for ms, me, mst in mat_hits[mn][0]:
                    aln = align_near(seq, tgt, ms, me, mst)
                    if aln is None:
                        continue
                    score, s, e, cols, matches = aln
                    cand = (-score, mst, s, e, cols, matches)
                    if best is None or cand < best:
                        best = cand
            if best:
                _, mst, s, e, cols, matches = best
                recs.append(Rec(name, acc, "precursor", s, e, mst, LOCAL,
                                round(100.0 * matches / cols, 1), cols, ["approximate"]))
        if not recs:
            recs = [Rec(name, acc, "precursor",
                        flags=[unplaced_reason(seq, kt, keep)])]
        for r in recs:
            if r.method != UNPLACED:
                suffix = {"+": "_plus", "-": "_minus"}[r.strand] if "strand_ambiguous" in r.flags else ""
                r.gene_id = r.tid = name + suffix
        pre_recs[name] = [r for r in recs if r.method != UNPLACED]
        out += recs
    for name, acc, seq in matures:
        hits, mmethod, mflags = mat_hits[name]
        recs = []
        for s, e, st in hits:
            for p in pre_recs.get(precursor_of(name), []):
                if p.strand == st and p.start <= s and e <= p.end:
                    suffix = {"+": "_plus", "-": "_minus"}[st] if "strand_ambiguous" in p.flags else ""
                    flags = [f for f in mflags if f == "length_changed"]
                    flags += [f for f in p.flags if f == "strand_ambiguous"]
                    recs.append(Rec(name, acc, "mature", s, e, st, mmethod,
                                    100.0 if mmethod == EXACT else None, e - s + 1, flags,
                                    p.gene_id, name + suffix))
                    break
        reason = "outside_precursor" if hits else unplaced_reason(seq, kt, keep)
        out += recs or [Rec(name, acc, "mature", flags=[reason])]
    return out


def flag_str(r: Rec) -> str:
    return ",".join(r.flags) or "none"


def gtf_rows(recs: list[Rec], chrom: str) -> list[list[str]]:
    rows = []
    for r in recs:
        if r.method == UNPLACED:
            continue
        s, e = str(r.start), str(r.end)
        gattrs = {"gene_id": r.gene_id, "gene_name": r.gene_id}
        if r.typ == "precursor":
            rows.append([chrom, SOURCE, "gene", s, e, ".", r.strand, ".",
                         format_attrs({**gattrs, "gene_biotype": "miRNA"})])
        t = {**gattrs, "transcript_id": r.tid,
             "transcript_biotype": "pre_miRNA" if r.typ == "precursor" else "miRNA",
             "mirbase_id": r.acc, "placement": r.method, "flags": flag_str(r)}
        for feat in ("transcript", "exon"):
            rows.append([chrom, SOURCE, feat, s, e, ".", r.strand, ".", format_attrs(t)])
    return sorted(rows, key=sort_key)


def bed_lines(recs: list[Rec], chrom: str) -> list[str]:
    placed = [r for r in recs if r.method != UNPLACED]
    placed.sort(key=lambda r: (r.start, r.end, r.tid, r.strand))
    return [f"{chrom}\t{r.start - 1}\t{r.end}\t{r.tid}\t0\t{r.strand}" for r in placed]


def report_lines(recs: list[Rec], chrom: str) -> list[str]:
    lines = ["\t".join(REPORT_COLS)]
    for r in recs:
        if r.method == UNPLACED:
            cols = [r.name, r.acc, r.typ, "NA", "NA", "NA", "NA", UNPLACED, "NA", "NA", flag_str(r)]
        else:
            ident = "NA" if r.identity is None else f"{r.identity:.1f}"
            cols = [r.tid, r.acc, r.typ, chrom, r.start, r.end, r.strand, r.method, ident,
                    r.aln_len, flag_str(r)]
        lines.append("\t".join(str(c) for c in cols))
    return lines


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hairpin", required=True, help="HSV-1 miRBase hairpin FASTA")
    ap.add_argument("--mature", required=True, help="HSV-1 miRBase mature FASTA")
    ap.add_argument("--tgt-fa", required=True, help="KT899744delRR FASTA")
    ap.add_argument("--src-fa", required=True, help="BK012101.1 FASTA")
    ap.add_argument("--kt-fa", help="KT899744.1 FASTA (to explain unplaced miRNAs)")
    ap.add_argument("--version", required=True, help="annotation version for the GTF header")
    ap.add_argument("--paf", required=True, help="BK012101.1 -> KT899744delRR minimap2 PAF")
    ap.add_argument("--out-gtf", required=True)
    ap.add_argument("--out-bed", required=True)
    ap.add_argument("--report", required=True)
    a = ap.parse_args()
    with open(a.tgt_fa) as fh:
        chrom = fh.readline()[1:].split()[0]
    tgt, src = read_fasta(a.tgt_fa), read_fasta(a.src_fa)
    pmap = build_map(a.paf, len(src)) if os.path.getsize(a.paf) else None
    kt = read_fasta(a.kt_fa) if a.kt_fa else ""
    recs = place(read_named_fasta(a.hairpin), read_named_fasta(a.mature), tgt, src, pmap, kt)
    write_gtf(a.out_gtf, gtf_rows(recs, chrom), [f"#!annotation-version {a.version}"])
    with open(a.out_bed, "w") as fh:
        fh.writelines(line + "\n" for line in bed_lines(recs, chrom))
    with open(a.report, "w") as fh:
        fh.writelines(line + "\n" for line in report_lines(recs, chrom))


if __name__ == "__main__":
    main()
