#!/usr/bin/env python3
"""Derive quantification and visualization GTFs plus overlap reports from the master GTF."""
from __future__ import annotations

import argparse
import csv
import os

from classify import annotated_lengths, classify, footprints
from gtflib import parse_attrs, read_gtf, write_gtf

METRICS = ["total_cds", "non_overlapping", "partially_overlapping", "embedded", "excluded",
           "retained", "retained_zero_unique", "mirna", "tss", "tes"]
EXCL_COLS = ["cds_id", "gene_name", "chrom", "start", "end", "strand", "cds_length",
             "containing_cds_id", "overlapping_bases", "unique_bases", "reason", "cds_check"]
OVL_COLS = ["gene_a", "gene_b", "strand_a", "strand_b", "same_strand", "overlap_bases",
            "class_a", "class_b", "length_a", "length_b", "unique_len_a", "unique_len_b"]


def _header(path: str) -> list[str]:
    with open(path) as fh:
        first = fh.readline().rstrip("\n")
    return [first] if first.startswith("#!annotation-version") else []


def _write_tsv(path: str, cols: list[str], rows: list[list]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        w.writerows(rows)


def derive(master: str, outdir: str, name: str) -> None:
    rows = read_gtf(master)
    header = _header(master)
    fp = footprints(rows)
    res = classify(fp, annotated_lengths(rows))
    embedded = {g for g, c in res.items() if c.cls == "embedded"}

    gene_name: dict[str, str] = {}
    chrom: dict[str, str] = {}
    noncoding: set[str] = set()
    cds_check: dict[str, str] = {}
    for r in rows:
        a = parse_attrs(r[8])
        gid = a.get("gene_id", "")
        if r[2] == "gene":
            gene_name[gid] = a.get("gene_name", gid)
            chrom[gid] = r[0]
            if a.get("gene_biotype") == "noncoding":
                noncoding.add(gid)
        elif r[2] == "transcript":
            v = a.get("cds_check", "ok")
            if cds_check.get(gid, "ok") == "ok":
                cds_check[gid] = v

    os.makedirs(outdir, exist_ok=True)
    pre = os.path.join(outdir, name)
    quant = [r for r in rows if r[2] not in ("TSS", "TES")
             and (g := parse_attrs(r[8]).get("gene_id", "")) not in embedded
             and g not in noncoding]
    write_gtf(f"{pre}.quantification.gtf", quant, header)
    write_gtf(f"{pre}.visualization.gtf", rows, header)

    excl = []
    for g in embedded:
        c = res[g]
        strand, bases = fp[g]
        excl.append([int(min(bases)), g, [
            g, gene_name.get(g, g), chrom.get(g, ""), min(bases), max(bases), strand,
            len(bases), c.host, len(bases & fp[c.host][1]), c.unique_len, c.reason,
            cds_check.get(g, "ok")]])
    excl.sort(key=lambda x: (x[0], x[1]))
    _write_tsv(f"{pre}.excluded_CDSs.tsv", EXCL_COLS, [x[2] for x in excl])

    ovl = []
    genes = sorted(fp)
    for i, a in enumerate(genes):
        for b in genes[i + 1:]:
            n = len(fp[a][1] & fp[b][1])
            if n:
                ovl.append([a, b, fp[a][0], fp[b][0], "yes" if fp[a][0] == fp[b][0] else "no",
                            n, res[a].cls, res[b].cls, res[a].length, res[b].length,
                            res[a].unique_len, res[b].unique_len])
    _write_tsv(f"{pre}.CDS_overlaps.tsv", OVL_COLS, ovl)

    count = {k: sum(c.cls == k for c in res.values())
             for k in ("non_overlapping", "partially_overlapping", "embedded",
                       "retained_zero_unique")}
    summary = {"total_cds": len(res), **count, "excluded": len(embedded),
               "retained": len(res) - len(embedded), "mirna": 0,
               "tss": sum(r[2] == "TSS" for r in rows), "tes": sum(r[2] == "TES" for r in rows)}
    _write_tsv(f"{pre}.annotation_summary.tsv", ["metric", "value"],
               [[m, summary[m]] for m in METRICS])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--master", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--name", required=True)
    a = ap.parse_args()
    derive(a.master, a.outdir, a.name)


if __name__ == "__main__":
    main()
