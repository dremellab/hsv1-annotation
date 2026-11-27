import hashlib
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from gtflib import parse_attrs, read_fasta, read_gtf, revcomp
from place_mirnas import EXACT, place

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "scripts" / "place_mirnas.py"
MB = ROOT / "sources" / "miRBase_23.1"
TGT = ROOT / "sources" / "KT899744delRR" / "KT899744delRR.fa"
SRC = ROOT / "sources" / "BK012101.1" / "BK012101.1.fa"
KT = ROOT / "sources" / "KT899744.1" / "KT899744.1.fa"
QUANT_MD5 = "dae588e5069a3fd0efd3b776085611af"  # committed v0.1.0 quantification GTF


def rnd(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def identity_map(n, off=0):
    """Toy build_map result: source pos i -> target pos i + off, whole block."""
    return [0] + [i + off for i in range(1, n + 1)] + [0], 1, n


def recs_of(res, name):
    return [r for r in res if r.name == name]


def test_palindrome_gives_both_strands_flagged():
    x = rnd(30, 1)
    hp = x + revcomp(x)                      # equals its own reverse complement
    assert hp == revcomp(hp)
    tgt = rnd(100, 2) + hp + rnd(100, 3)
    mat = hp[36:58]                          # 22 nt arm inside the precursor
    res = place([("hsv1-mir-H11", "MI1", hp.replace("T", "U"))],
                [("hsv1-miR-H11-3p", "MIMAT1", mat.replace("T", "U"))],
                tgt, "", None)
    pre = recs_of(res, "hsv1-mir-H11")
    assert sorted(r.strand for r in pre) == ["+", "-"]
    assert all(r.start == 101 and r.end == 160 for r in pre)
    assert all(r.method == "exact_target" and "strand_ambiguous" in r.flags for r in pre)
    mat_recs = recs_of(res, "hsv1-miR-H11-3p")
    assert sorted(r.strand for r in mat_recs) == ["+", "-"]
    plus = [r for r in mat_recs if r.strand == "+"][0]
    minus = [r for r in mat_recs if r.strand == "-"][0]
    assert plus.start == 101 + 36 and minus.start == 101 + 60 - 58  # mirrored arm


def test_antisense_pair_not_ambiguous():
    core = rnd(80, 4)
    tgt = rnd(50, 5) + core + rnd(50, 6)
    h1 = core[:50]
    h6 = revcomp(core[:60])                  # same start, opposite strand, different seq
    res = place([("hsv1-mir-H1", "MI1", h1), ("hsv1-mir-H6", "MI6", h6)],
                [("hsv1-miR-H1-5p", "MIMAT1", h1[5:27]), ("hsv1-miR-H6-3p", "MIMAT6", h6[10:32])],
                tgt, "", None)
    a, b = recs_of(res, "hsv1-mir-H1"), recs_of(res, "hsv1-mir-H6")
    assert len(a) == len(b) == 1
    assert (a[0].strand, a[0].start) == ("+", 51)
    assert (b[0].strand, b[0].start) == ("-", 51)
    assert not a[0].flags and not b[0].flags


def test_projection_path():
    src = rnd(200, 7)
    hp = src[50:100]
    tgt = rnd(190, 8)                        # hairpin absent from target
    assert hp not in tgt and revcomp(hp) not in tgt
    res = place([("hsv1-mir-H9", "MI9", hp)], [("hsv1-miR-H9-5p", "MIMAT9", hp[3:25])],
                tgt, src, identity_map(200, -10))
    p = recs_of(res, "hsv1-mir-H9")[0]
    assert p.method == "projected_from_BK012101.1"
    assert (p.start, p.end, p.strand) == (41, 90, "+")
    assert not p.flags
    m = recs_of(res, "hsv1-miR-H9-5p")
    assert len(m) == 1 and m[0].method == "projected_from_BK012101.1" and m[0].start == 44


def test_projection_length_change_recorded():
    src = rnd(200, 9)
    hp = src[50:100]
    m = [0] * 202
    for i in range(1, 201):
        m[i] = i - 10 if i <= 80 else i - 12  # 2 nt deletion inside the hairpin
    res = place([("hsv1-mir-H9", "MI9", hp)], [], rnd(190, 10), src, (m, 1, 200))
    p = recs_of(res, "hsv1-mir-H9")[0]
    assert p.end - p.start + 1 == 48


def test_local_alignment_path():
    tgt = rnd(400, 11)
    hp = list(tgt[150:230])
    for i in (5, 70):                        # 2 mismatches, outside the mature arm
        hp[i] = "A" if hp[i] != "A" else "C"
    hp = "".join(hp)
    mat = hp[20:42]
    res = place([("hsv1-mir-H8", "MI8", hp)], [("hsv1-miR-H8-5p", "MIMAT8", mat)],
                tgt, "", None)
    p = recs_of(res, "hsv1-mir-H8")[0]
    assert p.method == "local_alignment" and "approximate" in p.flags
    assert p.identity < 100 and p.aln_len == 80
    assert (p.start, p.end, p.strand) == (151, 230, "+")
    m = recs_of(res, "hsv1-miR-H8-5p")
    assert len(m) == 1 and m[0].start == 151 + 20


def test_unplaced_no_match():
    res = place([("hsv1-mir-H99", "MI99", rnd(60, 12))], [], rnd(300, 13), "", None)
    assert [(r.name, r.method, r.flags) for r in res] == [
        ("hsv1-mir-H99", "unplaced", ["no_match"])]


def test_unplaced_in_deleted_region():
    kt = rnd(500, 16)
    keep = [(1, 200), (301, 500)]            # 201-300 deleted
    tgt = kt[:200] + kt[300:]
    inside_del, straddle = kt[220:280], kt[180:240]  # straddling hit is not "inside"
    res = place([("hsv1-mir-H5", "MI5", inside_del), ("hsv1-mir-H6", "MI6", straddle)],
                [("hsv1-miR-H5-5p", "MIMAT5", inside_del[5:27])], tgt, "", None, kt, keep)
    assert {r.name: (r.method, r.flags) for r in res} == {
        "hsv1-mir-H5": ("unplaced", ["in_deleted_region"]),
        "hsv1-mir-H6": ("unplaced", ["no_match"]),
        "hsv1-miR-H5-5p": ("unplaced", ["in_deleted_region"])}


def test_mature_own_method_inherits_only_strand_ambiguous():
    src = rnd(200, 17)
    hp = src[50:100]
    tgt = rnd(190, 18)
    tgt = tgt[:60] + hp[5:27] + tgt[82:]     # mature exact on target, hairpin only via source
    res = place([("hsv1-mir-H9", "MI9", hp)], [("hsv1-miR-H9-5p", "MIMAT9", hp[5:27])],
                tgt, src, identity_map(200, -10))
    m = recs_of(res, "hsv1-miR-H9-5p")[0]
    assert (m.method, m.identity, m.flags) == (EXACT, 100.0, [])
    p = recs_of(res, "hsv1-mir-H9")[0]
    assert (p.method, p.identity) == ("projected_from_BK012101.1", None)


def test_mature_kept_only_inside_precursor_same_strand():
    tgt = rnd(400, 14)
    hp = tgt[100:160]
    mat = hp[10:32]
    tgt = tgt[:300] + mat + tgt[322:]        # extra copy outside the precursor
    tgt = tgt[:350] + revcomp(mat) + tgt[372:]  # and one on the other strand
    res = place([("hsv1-mir-H3", "MI3", hp)], [("hsv1-miR-H3-5p", "MIMAT3", mat)], tgt, "", None)
    m = recs_of(res, "hsv1-miR-H3-5p")
    assert [(r.start, r.strand) for r in m] == [(111, "+")]


def write_fa(path, name, seq):
    Path(path).write_text(f">{name} x\n{seq}\n")


def test_cli_outputs(tmp_path):
    x = rnd(30, 1)
    hp = x + revcomp(x)
    tgt = rnd(100, 2) + hp + rnd(100, 3)
    write_fa(tmp_path / "t.fa", "KT899744", tgt)
    write_fa(tmp_path / "s.fa", "BK", rnd(300, 15))
    (tmp_path / "hp.fa").write_text(f">hsv1-mir-H11 MI0014 Herpes Simplex Virus 1 miR-H11 stem-loop\n{hp.replace('T','U')}\n")
    (tmp_path / "m.fa").write_text(f">hsv1-miR-H11-3p MIMAT0014688 Herpes Simplex Virus 1 miR-H11-3p\n{hp[40:62].replace('T','U')}\n")
    (tmp_path / "x.paf").write_text("")
    outs = []
    for k in (1, 2):
        o = tmp_path / f"o{k}"
        o.mkdir()
        subprocess.run([sys.executable, str(SCRIPT), "--hairpin", str(tmp_path / "hp.fa"),
                        "--mature", str(tmp_path / "m.fa"), "--tgt-fa", str(tmp_path / "t.fa"),
                        "--src-fa", str(tmp_path / "s.fa"), "--paf", str(tmp_path / "x.paf"),
                        "--version", "v9", "--out-gtf", str(o / "a.gtf"), "--out-bed", str(o / "a.bed"),
                        "--report", str(o / "a.tsv")], check=True, capture_output=True)
        outs.append(o)
    for f in ("a.gtf", "a.bed", "a.tsv"):
        assert (outs[0] / f).read_bytes() == (outs[1] / f).read_bytes()
    assert (outs[0] / "a.gtf").read_text().splitlines()[0] == "#!annotation-version v9"
    rows = read_gtf(str(outs[0] / "a.gtf"))
    assert {r[1] for r in rows} == {"miRBase_23.1"} and {r[0] for r in rows} == {"KT899744"}
    genes = {parse_attrs(r[8])["gene_id"]: parse_attrs(r[8]) for r in rows if r[2] == "gene"}
    assert set(genes) == {"hsv1-mir-H11_plus", "hsv1-mir-H11_minus"}
    assert all(a["gene_biotype"] == "miRNA" for a in genes.values())
    tx = {parse_attrs(r[8])["transcript_id"]: parse_attrs(r[8]) for r in rows if r[2] == "transcript"}
    assert set(tx) == {"hsv1-mir-H11_plus", "hsv1-mir-H11_minus",
                       "hsv1-miR-H11-3p_plus", "hsv1-miR-H11-3p_minus"}
    a = tx["hsv1-miR-H11-3p_plus"]
    assert a["transcript_biotype"] == "miRNA" and a["mirbase_id"] == "MIMAT0014688"
    assert a["placement"] == "exact_target" and a["flags"] == "strand_ambiguous"
    assert a["gene_id"] == "hsv1-mir-H11_plus"
    assert tx["hsv1-mir-H11_minus"]["transcript_biotype"] == "pre_miRNA"
    assert not any(" " in v for at in tx.values() for v in at.values())
    bed = (outs[0] / "a.bed").read_text().splitlines()
    assert len(bed) == 4 and all(len(b.split("\t")) == 6 for b in bed)
    rep = (outs[0] / "a.tsv").read_text().splitlines()
    assert rep[0].split("\t") == ["name", "mirbase_id", "type", "chrom", "start", "end",
                                  "strand", "placement", "identity_pct", "aln_length", "flags"]
    assert len(rep) == 5


needs = pytest.mark.skipif(not MB.exists() or not SRC.exists() or shutil.which("minimap2") is None,
                           reason="needs sources and minimap2 on PATH")


def fa_headers(path):
    return [l[1:].split()[0] for l in Path(path).read_text().splitlines() if l.startswith(">")]


@needs
def test_real_data(paf, tmp_path):
    o = tmp_path
    subprocess.run([sys.executable, str(SCRIPT), "--hairpin", str(MB / "hsv1_hairpin.fa"),
                    "--mature", str(MB / "hsv1_mature.fa"), "--tgt-fa", str(TGT),
                    "--src-fa", str(SRC), "--kt-fa", str(KT), "--paf", str(paf), "--version", "v0",
                    "--out-gtf", str(o / "m.gtf"),
                    "--out-bed", str(o / "m.bed"), "--report", str(o / "m.tsv")],
                   check=True, capture_output=True)
    rep = [l.split("\t") for l in (o / "m.tsv").read_text().splitlines()[1:]]
    pre = [r for r in rep if r[2] == "precursor"]
    mat = [r for r in rep if r[2] == "mature"]
    unpl = {r[0]: r[10] for r in rep if r[7] == "unplaced"}
    assert unpl == {"hsv1-mir-H5": "in_deleted_region", "hsv1-miR-H5-5p": "in_deleted_region",
                    "hsv1-miR-H5-3p": "in_deleted_region"}
    assert {r[0].replace("_plus", "").replace("_minus", "") for r in pre} == set(fa_headers(MB / "hsv1_hairpin.fa")) and len(pre) == 19
    assert {r[0].replace("_plus", "").replace("_minus", "") for r in mat} == set(fa_headers(MB / "hsv1_mature.fa")) and len(mat) == 28
    h11 = [r for r in pre if r[0].startswith("hsv1-mir-H11")]
    assert sorted(r[6] for r in h11) == ["+", "-"] and {r[4] for r in h11} == {"53208"}
    assert {r[7] for r in pre if r[0] in ("hsv1-mir-H8", "hsv1-mir-H17")} == {"local_alignment"}
    h1 = [r for r in pre if r[0] == "hsv1-mir-H1"][0]
    h6 = [r for r in pre if r[0] == "hsv1-mir-H6"][0]
    assert (h1[4], h1[6], h6[4], h6[6]) == ("109078", "+", "109078", "-")
    assert not h1[10].count("strand_ambiguous") and not h6[10].count("strand_ambiguous")
    arms = {(r[6], r[4]) for r in mat if r[0].startswith("hsv1-miR-H11-3p")}
    assert arms == {("+", "53280"), ("-", "53246")}


@needs
def test_real_quant_gtf_unchanged(paf, tmp_path):
    tgt = str(TGT)
    scripts = ROOT / "scripts"
    run = lambda *a: subprocess.run([sys.executable, *map(str, a)], check=True, capture_output=True)
    run(scripts / "lift_bk012101.py", "--gff", SRC.parent / "BK012101.1.gff3", "--src-fa", SRC,
        "--tgt-fa", tgt, "--paf", paf, "--out-gtf", tmp_path / "l.gtf",
        "--report", tmp_path / "r.tsv", "--date", "2026.11.10", "--keep-noncoding")
    # same --date as build.sh's VERSION is not needed: quant rows carry no date
    run(scripts / "place_mirnas.py", "--hairpin", MB / "hsv1_hairpin.fa", "--mature",
        MB / "hsv1_mature.fa", "--tgt-fa", tgt, "--src-fa", SRC, "--kt-fa", KT, "--paf", paf, "--version", "v0.1.0",
        "--out-gtf", tmp_path / "m.gtf", "--out-bed", tmp_path / "m.bed", "--report", tmp_path / "m.tsv")
    run(scripts / "build_master.py", "--in-gtf", tmp_path / "l.gtf", "--mirna-gtf", tmp_path / "m.gtf",
        "--out-gtf", tmp_path / "master.gtf", "--version", "v0.1.0", "--fasta", tgt)
    run(scripts / "derive_gtfs.py", "--master", tmp_path / "master.gtf", "--outdir", tmp_path / "o",
        "--name", "KT899744delRR", "--mirna-report", tmp_path / "m.tsv")
    summ = dict(l.split("\t") for l in (tmp_path / "o" / "KT899744delRR.annotation_summary.tsv"
                                        ).read_text().splitlines()[1:])
    assert (summ["mirna"], summ["mirna_precursors"], summ["mirna_unplaced"]) == ("26", "18", "3")
    assert hashlib.md5((tmp_path / "o" / "KT899744delRR.quantification.gtf").read_bytes()
                       ).hexdigest() == QUANT_MD5
