#!/usr/bin/env python3
"""Lift the BK012101.1 (HSV-1 strain 17) annotation onto KT899744delRR.

Coordinates are projected through the collinear minimap2 alignment block
(PAF with cg:Z CIGAR, --eqx). Features lying in the terminal repeats, which
are deleted in KT899744delRR, are dropped (their internal-repeat copies are
annotated separately in BK012101.1); features spanning a terminal repeat
boundary are clipped. Every lifted CDS is re-translated on the target and
compared with the source protein.
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import OrderedDict, defaultdict
from urllib.parse import unquote

from gtflib import read_fasta, revcomp

CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            CODON[a + b + c] = _aa[i * 16 + j * 4 + k]


def translate(nt: str) -> str:
    return "".join(CODON.get(nt[i : i + 3], "X") for i in range(0, len(nt) - 2, 3))


def attrs(s: str) -> dict:
    d = {}
    for kv in s.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = unquote(v)
    return d


def clean(v: str) -> str:
    # no spaces: _fix_gtf.py (RIBOLD and HAROLD) splits attributes on whitespace
    return v.replace('"', "'").replace(";", ",").replace(" ", "_")


def build_map(paf: str, qlen: int):
    """1-based source pos -> 1-based target pos for the primary + strand block."""
    best = None
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        if f[4] != "+":
            continue
        if best is None or int(f[9]) > int(best[9]):
            best = f
    cg = [x for x in best[12:] if x.startswith("cg:Z:")][0][5:]
    q, t = int(best[2]), int(best[7])
    m = [0] * (qlen + 2)
    for n, op in re.findall(r"(\d+)([=XIDM])", cg):
        n = int(n)
        if op in "=XM":
            for i in range(n):
                m[q + i + 1] = t + i + 1
            q += n
            t += n
        elif op == "I":
            q += n
        elif op == "D":
            t += n
    return m, int(best[2]) + 1, int(best[3])


def seg_seq(genome: str, segs, strand: str) -> str:
    s = "".join(genome[a - 1 : b] for a, b in sorted(segs))
    return revcomp(s) if strand == "-" else s


def take(segs, strand, n, from5=True):
    """Return sub-segments covering n nt at the 5' (or 3') end, and the rest."""
    order = sorted(segs, reverse=(strand == "-"))
    if not from5:
        order = order[::-1]
    got, rest = [], []
    for a, b in order:
        if n <= 0:
            rest.append((a, b))
            continue
        L = b - a + 1
        if L <= n:
            got.append((a, b))
            n -= L
        else:
            lowside = (strand == "+") == from5
            if lowside:
                got.append((a, a + n - 1))
                rest.append((a + n, b))
            else:
                got.append((b - n + 1, b))
                rest.append((a, b - n))
            n = 0
    return sorted(got), sorted(rest)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--src-fa", required=True)
    ap.add_argument("--tgt-fa", required=True)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--out-gtf", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--chrom", default="KT899744")
    ap.add_argument("--prefix", default="KT899744delRR_")
    ap.add_argument("--date", required=True)
    ap.add_argument("--keep-noncoding", action="store_true",
                    help="write transcripts without a CDS as noncoding genes (<prefix><transcript_id>_nc)")
    a = ap.parse_args()

    src, tgt = read_fasta(a.src_fa), read_fasta(a.tgt_fa)
    pmap, qlo, qhi = build_map(a.paf, len(src))

    genes: OrderedDict = OrderedDict()
    mrnas: OrderedDict = OrderedDict()
    cds: OrderedDict = OrderedDict()
    for line in open(a.gff):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) != 9:
            continue
        typ, s, e, strand, ph, at = f[2], int(f[3]), int(f[4]), f[6], f[7], attrs(f[8])
        if typ == "gene":
            genes[at["ID"]] = dict(name=at["Name"], strand=strand, s=s, e=e, biotype=at.get("gene_biotype", "protein_coding"))
        elif typ == "mRNA":
            mrnas[at["ID"]] = dict(gene=at["Parent"], at=at, exons=[], strand=strand)
        elif typ == "exon":
            mrnas[at["Parent"]]["exons"].append((s, e))
        elif typ == "CDS":
            key = (at["ID"], at["Parent"])
            c = cds.setdefault(key, dict(parent=at["Parent"], at=at, segs=[], strand=strand))
            c["segs"].append((s, e, int(ph)))

    # assemble transcripts
    txs = []
    cds_by_rna = defaultdict(list)
    for (cid, par), c in cds.items():
        if par.startswith("rna-"):
            cds_by_rna[par].append(c)
    for rid, m in mrnas.items():
        cl = cds_by_rna.get(rid, [])
        assert len(cl) <= 1, rid
        txs.append(dict(gene=m["gene"], tid=rid[4:], name=m["at"].get("product", rid[4:]), note=m["at"].get("Note", ""),
                        strand=m["strand"], exons=sorted(m["exons"]), cds=cl[0] if cl else None, kind="mRNA"))
    for (cid, par), c in cds.items():
        if par.startswith("gene-"):
            txs.append(dict(gene=par, tid=c["at"]["protein_id"], name=c["at"].get("product", ""), note="",
                            strand=c["strand"], exons=sorted((s, e) for s, e, _ in c["segs"]), cds=c, kind="CDS_only"))
    # one gene per ORF: every CDS record becomes its own gene, named by the
    # paper's CDS product; the paper's transcription-unit gene is kept as "locus".
    # Transcripts without a CDS keep the locus as gene and are dropped below.
    for t in txs:
        t["locus"] = genes[t["gene"]]["name"]
        if t["cds"]:
            orf = t["cds"]["at"]["product"]
            key = "orf-%s|%s" % (orf, t["cds"]["at"]["protein_id"])
            genes[key] = dict(name=orf, strand=t["strand"], biotype="protein_coding")
            t["gene"] = key
    ids = [t["tid"] for t in txs]
    assert len(ids) == len(set(ids)), "duplicate transcript ids"

    def lift_seg(s, e):
        pos = [pmap[p] for p in range(s, e + 1) if pmap[p]]
        outside = sum(1 for p in range(s, e + 1) if p < qlo or p > qhi)
        if not pos:
            return None, outside, 0
        lo, hi = min(pos), max(pos)
        # net length change inside the chain (target span minus source bases inside chain)
        delta = (hi - lo + 1) - ((e - s + 1) - outside)
        return (lo, hi), outside, delta

    rep = []
    out_by_gene = defaultdict(list)
    for t in txs:
        strand = t["strand"]
        g = genes[t["gene"]]
        ex_l, outside, delta, lost = [], 0, 0, 0
        for s, e in t["exons"]:
            seg, o, d = lift_seg(s, e)
            outside += o
            if seg is None:
                lost += 1
                continue
            delta += d
            ex_l.append(seg)
        row = dict(gkey=t["gene"], locus=t["locus"], gene=g["name"], tid=t["tid"], product=t["name"], kind=t["kind"], strand=strand,
                   src="%d-%d" % (t["exons"][0][0], t["exons"][-1][1]), tgt="", status="", cds="", pid="",
                   aa_src="", aa_tgt="", aa_diff="", detail="")
        if t["cds"]:
            row["pid"] = t["cds"]["at"]["protein_id"]
        if not ex_l:
            row["status"] = "dropped_terminal_repeat"
            rep.append(row)
            continue
        status = "exact" if (outside == 0 and delta == 0) else ("clipped_terminal_repeat" if outside else "indel(%+d)" % delta)
        if outside and delta:
            status += ";indel(%+d)" % delta
        row["status"] = status
        row["tgt"] = "%d-%d" % (min(x[0] for x in ex_l), max(x[1] for x in ex_l))

        cds_l, stop_l, start_l, frame0, cds_flag = [], [], [], 0, ""
        if t["cds"]:
            c = t["cds"]
            at = c["at"]
            segs = sorted((s, e) for s, e, _ in c["segs"])
            order = sorted(c["segs"], reverse=(strand == "-"))
            ph0 = order[0][2]
            src_nt = seg_seq(src, segs, strand)[ph0:]
            src_aa = translate(src_nt)
            lost5 = 0
            done5 = False
            for s, e, _ in order:
                seg, o, d = lift_seg(s, e)
                if seg is None:
                    if not done5:
                        lost5 += e - s + 1
                    continue
                if not done5:
                    # bases clipped off the 5' end of the first surviving segment
                    rng = range(s, e + 1) if strand == "+" else range(e, s - 1, -1)
                    for p in rng:
                        if pmap[p]:
                            break
                        lost5 += 1
                    done5 = True
                cds_l.append(seg)
            src_len = sum(e - s + 1 for s, e in segs)
            n_out = sum(1 for s, e in segs for p in range(s, e + 1) if p < qlo or p > qhi)
            if not cds_l:
                cds_flag = "CDS_lost"
            else:
                cds_l = sorted(cds_l)
                frame0 = (ph0 - lost5) % 3 if n_out else ph0
                tnt = seg_seq(tgt, cds_l, strand)[frame0:]
                taa = translate(tnt)
                p5 = "start_range" in at if strand == "+" else "end_range" in at
                p3 = "end_range" in at if strand == "+" else "start_range" in at
                tlen = sum(e - s + 1 for s, e in cds_l)
                flags = []
                if n_out:
                    flags.append("CDS_clipped(%dnt)" % n_out)
                elif (tlen - src_len) % 3:
                    flags.append("frameshift(%+dnt)" % (tlen - src_len))
                elif tlen != src_len:
                    flags.append("inframe_indel(%+dnt)" % (tlen - src_len))
                body = taa[:-1] if taa.endswith("*") else taa
                if "*" in body and not (tlen - src_len) % 3:
                    flags.append("internal_stop(aa%d)" % (body.index("*") + 1))
                elif "*" in body:
                    flags.append("premature_stop(aa%d)" % (body.index("*") + 1))
                if not p3 and not n_out and not taa.endswith("*"):
                    flags.append("no_stop")
                src_body = src_aa.rstrip("*")
                if len(body) == len(src_body):
                    nd = sum(1 for x, y in zip(body, src_body) if x != y)
                    row["aa_diff"] = str(nd)
                    if not flags:
                        flags.append("identical_aa" if nd == 0 else "aa_subs")
                if not p5 and not n_out and frame0 == 0 and tnt[:3] != src_nt[:3]:
                    flags.append("start_codon_changed(%s>%s)" % (src_nt[:3], tnt[:3]))
                cds_flag = ";".join(flags)
                row["aa_src"], row["aa_tgt"] = str(len(src_body)), str(len(body))
                # GTF convention: stop codon is not part of the CDS
                clip5 = bool(n_out) and lost5 > 0
                if taa.endswith("*") and (len(tnt) % 3 == 0):
                    stop_l, cds_l = take(cds_l, strand, 3, from5=False)
                if not p5 and not clip5 and frame0 == 0 and cds_l:
                    start_l, _ = take(cds_l, strand, 3, from5=True)
            row["cds"] = cds_flag
        rep.append(row)
        out_by_gene[t["gene"]].append(dict(t=t, exons=sorted(ex_l), cds=cds_l, stop=stop_l, start=start_l,
                                          frame0=frame0, status=status, cds_flag=cds_flag))

    # every gene must carry at least one CDS: drop genes with no coding transcript
    # (with --keep-noncoding each such transcript becomes its own noncoding gene)
    for gid in list(out_by_gene):
        if not any(x["cds"] for x in out_by_gene[gid]):
            members = out_by_gene.pop(gid)
            if a.keep_noncoding:
                for x in members:
                    tid = x["t"]["tid"]
                    nkey = "nc|" + tid
                    genes[nkey] = dict(name=tid + "_nc", strand=x["t"]["strand"], biotype="noncoding")
                    x["status"] = ("clipped_noncoding" if x["status"].startswith("clipped_terminal_repeat")
                                   else "lifted_noncoding")
                    out_by_gene[nkey].append(x)
                    for r in rep:
                        if r["tid"] == tid and r["tgt"]:
                            r["gkey"], r["gene"], r["status"] = nkey, tid + "_nc", x["status"]
                continue
            for r in rep:
                if r["gkey"] == gid and r["tgt"]:
                    r["status"] = "dropped_no_CDS(" + r["status"].replace("(", "[").replace(")", "]") + ")"

    # Fold the paper's variant records ("X *1", "X #1", "X iso1", optionally
    # "X RNA *1") into the gene of their base ORF "X" when their CDS overlap on
    # the same strand. Otherwise a base ORF nested in its variant would share
    # every CDS base with another gene and receive no uniquely assigned reads.
    def base_name(n):
        return re.sub(r"( RNA)? ([*#]\d+|iso\d+)$", "", n)

    groups = defaultdict(list)
    for gid in out_by_gene:
        groups[(base_name(genes[gid]["name"]), genes[gid]["strand"])].append(gid)
    for (base, strand), members in groups.items():
        if len(members) == 1 and genes[members[0]]["name"] == base:
            continue
        comp = {m: m for m in members}

        def find(x):
            while comp[x] != x:
                x = comp[x]
            return x

        def segs(m):
            return [s for x in out_by_gene[m] for s in x["cds"] + x["stop"]]

        for i, m1 in enumerate(members):
            for m2 in members[i + 1:]:
                if any(a1 <= b2 and a2 <= b1 for a1, b1 in segs(m1) for a2, b2 in segs(m2)):
                    comp[find(m1)] = find(m2)
        comps = defaultdict(list)
        for m in members:
            comps[find(m)].append(m)
        cl = sorted(comps.values(), key=lambda c: (not any(genes[m]["name"] == base for m in c), -len(c)))
        keep = cl[0]  # the component holding the base ORF (or the largest one) takes the base name
        if len(keep) == 1 and genes[keep[0]]["name"] != base and len(cl) == 1:
            pass  # lone variant with no base ORF: named after the base
        head = keep[0]
        for m in keep[1:]:
            out_by_gene[head].extend(out_by_gene.pop(m))
        genes[head] = dict(genes[head], name=base)
        for r in rep:
            if r["gkey"] in keep:
                r["gene"] = base

    names = [genes[gid]["name"] for gid in out_by_gene]
    dup = sorted({n for n in names if names.count(n) > 1})
    assert not dup, "ORF names not unique after lifting: %s" % dup

    lines = []
    glist = []
    for gid, g in genes.items():
        ts = out_by_gene.get(gid)
        if not ts:
            continue
        gs = min(x["exons"][0][0] for x in ts)
        ge = max(x["exons"][-1][1] for x in ts)
        glist.append((gs, ge, gid))
    for gs, ge, gid in sorted(glist):
        g = genes[gid]
        strand = g["strand"]
        gtype = "protein_coding" if g["biotype"] == "protein_coding" else "ncRNA"
        gname = g["name"].replace(" ", "_")
        gattr = 'gene_id "%s%s"; gene_name "%s%s"; gene_type "%s";' % (a.prefix, gname, a.prefix, gname, gtype)
        if g["biotype"] == "noncoding":
            gattr += ' gene_biotype "noncoding";'
        lines.append([a.chrom, "Genbank", "gene", gs, ge, ".", strand, ".", gattr])
        for x in sorted(out_by_gene[gid], key=lambda x: (x["exons"][0][0], x["exons"][-1][1], x["t"]["tid"])):
            t = x["t"]
            ttype = "protein_coding" if x["cds"] else ("ncRNA" if gtype == "ncRNA" else "processed_transcript")
            base = 'gene_id "%s%s"; transcript_id "%s"; gene_name "%s%s"; transcript_name "%s";' % (
                a.prefix, gname, clean(t["tid"]), a.prefix, gname, clean(t["name"]))
            tail = ' gene_type "%s"; transcript_type "%s"; locus "%s";' % (gtype, ttype, clean(t["locus"]))
            if x["cds"]:
                tail += ' transcript_biotype "protein_coding";'
            elif g["biotype"] == "noncoding":
                tail += ' transcript_biotype "noncoding"; gene_biotype "noncoding";'
            if t["cds"]:
                tail += ' protein_id "%s";' % t["cds"]["at"]["protein_id"]
            extra = ' transcript_support "%s"; liftover_status "%s";' % (
                "BK012101.1_mRNA" if t["kind"] == "mRNA" else "BK012101.1_CDS_only", x["status"])
            if x["cds_flag"]:
                extra += ' cds_check "%s";' % x["cds_flag"]
            note = t["note"] or (t["cds"]["at"].get("Note", "") if t["cds"] else "")
            if note:
                extra += ' note "%s";' % clean(note)
            lines.append([a.chrom, "Genbank", "transcript", x["exons"][0][0], x["exons"][-1][1], ".", strand, ".", base + tail + extra])
            order = sorted(x["exons"], reverse=(strand == "-"))
            for n, (s, e) in enumerate(order, 1):
                lines.append([a.chrom, "Genbank", "exon", s, e, ".", strand, ".", base + ' exon_number "%d";' % n + tail])
            if x["cds"]:
                done = x["frame0"]
                first = True
                for s, e in sorted(x["cds"], reverse=(strand == "-")):
                    fr = x["frame0"] if first else (3 - (done % 3)) % 3
                    if first:
                        done = (e - s + 1) - x["frame0"]
                    else:
                        done += e - s + 1
                    first = False
                    lines.append([a.chrom, "Genbank", "CDS", s, e, ".", strand, str(fr), base + tail])
            for s, e in x["start"]:
                lines.append([a.chrom, "Genbank", "start_codon", s, e, ".", strand, "0", base + tail])
            for s, e in x["stop"]:
                lines.append([a.chrom, "Genbank", "stop_codon", s, e, ".", strand, "0", base + tail])

    with open(a.out_gtf, "w") as fh:
        fh.write("##description: strain KOS HSV-1 genome annotation (KT899744.1 modified so only one copy of the joint region "
                 "was present: Δ1-9603, Δ125845-126977, Δ145361-151974); gene models lifted from the strain 17 annotation "
                 "GenBank BK012101.1 (Whisnant et al. 2020, PMID 32341360) through a minimap2 whole-genome alignment; "
                 "one gene per ORF (gene_id = paper CDS product, locus = paper transcription-unit gene); "
                 "terminal-repeat copies and transcripts without a CDS dropped; CDS excludes the stop codon\n")
        fh.write("##provider: Genbank BK012101.1, lifted to KT899744delRR\n##format: gtf\n##date: %s\n" % a.date)
        for l in lines:
            fh.write("\t".join(map(str, l)) + "\n")
    cols = ["gene", "locus", "tid", "product", "kind", "strand", "src", "tgt", "status", "pid", "cds", "aa_src", "aa_tgt", "aa_diff"]
    with open(a.report, "w") as fh:
        fh.write("\t".join(["orf_gene", "locus", "transcript_id", "product", "kind", "strand", "BK012101.1_coords", "KT899744delRR_coords",
                            "liftover_status", "protein_id", "cds_check", "aa_len_src", "aa_len_tgt", "aa_diffs"]) + "\n")
        for r in rep:
            fh.write("\t".join(r[c] for c in cols) + "\n")
    print("chain: source %d-%d" % (qlo, qhi), file=sys.stderr)


if __name__ == "__main__":
    main()
