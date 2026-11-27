import csv
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from check_gtf import check
from gtflib import parse_attrs, read_fasta, read_gtf

P = "KT899744delRR_"
ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "sources" / "BK012101.1"
SCRIPT = ROOT / "scripts" / "derive_gtfs.py"
NAME = "KT899744delRR"


def row(feat, s, e, attrs, strand="+", frame="."):
    return ["KT899744", "Genbank", feat, str(s), str(e), ".", strand, frame, attrs]


def coding(name, s, e, strand="+", cds_check="ok"):
    g = P + name
    ga = f'gene_id "{g}"; gene_name "{name}";'
    ta = f'gene_id "{g}"; transcript_id "{g}_t"; gene_name "{name}"; cds_check "{cds_check}";'
    return [row("gene", s, e, ga, strand), row("transcript", s, e, ta, strand),
            row("exon", s, e, ta, strand), row("CDS", s, e, ta, strand, "0")]


def write_master(path):
    rows = []
    rows += coding("A", 100, 300)
    rows += coding("B", 400, 600)
    rows += coding("C", 150, 250, cds_check="stop_missing")  # embedded in A
    rows += coding("D", 200, 350, strand="-")                # opposite strand over A
    rows += [row("gene", 700, 800, f'gene_id "{P}N"; gene_name "N"; gene_biotype "noncoding";'),
             row("transcript", 700, 800, f'gene_id "{P}N"; transcript_id "{P}N_t";'),
             row("exon", 700, 800, f'gene_id "{P}N"; transcript_id "{P}N_t";')]
    rows += [row("TSS", 100, 100, f'gene_id "{P}A"; transcript_id "{P}A_t";'),
             row("TES", 300, 300, f'gene_id "{P}A"; transcript_id "{P}A_t";')]
    with open(path, "w") as fh:
        fh.write("#!annotation-version v9\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return rows


@pytest.fixture
def out(tmp_path):
    m = tmp_path / "master.gtf"
    write_master(m)
    o = tmp_path / "out"
    subprocess.run([sys.executable, str(SCRIPT), "--master", str(m), "--outdir", str(o),
                    "--name", NAME], check=True, capture_output=True)
    return o


def tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_quant_excludes_embedded_noncoding_tss(out):
    text = (out / f"{NAME}.quantification.gtf").read_text()
    assert text.splitlines()[0] == "#!annotation-version v9"
    rows = read_gtf(str(out / f"{NAME}.quantification.gtf"))
    genes = {r[8].split('"')[1] for r in rows}
    assert genes == {P + "A", P + "B", P + "D"}
    assert not {r[2] for r in rows} & {"TSS", "TES"}
    assert "noncoding" not in text


def test_vis_equals_master(out, tmp_path):
    assert (out / f"{NAME}.visualization.gtf").read_text() == (tmp_path / "master.gtf").read_text()


def test_excluded_columns_and_values(out):
    path = out / f"{NAME}.excluded_CDSs.tsv"
    assert path.read_text().splitlines()[0].split("\t") == [
        "cds_id", "gene_name", "chrom", "start", "end", "strand", "cds_length",
        "containing_cds_id", "overlapping_bases", "unique_bases", "reason", "cds_check"]
    rows = tsv(path)
    assert rows == [{
        "cds_id": P + "C", "gene_name": "C", "chrom": "KT899744", "start": "150",
        "end": "250", "strand": "+", "cds_length": "101", "containing_cds_id": P + "A",
        "overlapping_bases": "101", "unique_bases": "0",
        "reason": "embedded_in_single_orf", "cds_check": "stop_missing"}]


def test_summary_arithmetic(out):
    s = {r["metric"]: int(r["value"]) for r in tsv(out / f"{NAME}.annotation_summary.tsv")}
    assert list(s) == ["total_cds", "non_overlapping", "partially_overlapping", "embedded",
                       "excluded", "retained", "retained_zero_unique", "mirna", "tss", "tes"]
    assert s["total_cds"] == 4 and s["embedded"] == 1 and s["tss"] == 1 and s["tes"] == 1
    assert s["total_cds"] == (s["non_overlapping"] + s["partially_overlapping"]
                              + s["embedded"] + s["retained_zero_unique"])
    assert s["excluded"] == s["embedded"]
    assert s["retained"] == s["total_cds"] - s["excluded"]
    assert s["mirna"] == 0


def test_overlaps_includes_opposite_strand(out):
    rows = tsv(out / f"{NAME}.CDS_overlaps.tsv")
    pairs = [(r["gene_a"], r["gene_b"]) for r in rows]
    assert pairs == sorted(pairs)
    d = [r for r in rows if r["gene_a"] == P + "A" and r["gene_b"] == P + "D"][0]
    assert d["same_strand"] == "no" and d["strand_a"] == "+" and d["strand_b"] == "-"
    assert d["overlap_bases"] == "101" and d["length_a"] == "201"
    ac = [r for r in rows if r["gene_a"] == P + "A" and r["gene_b"] == P + "C"][0]
    assert ac["same_strand"] == "yes" and ac["class_b"] == "embedded"
    assert not [r for r in rows if P + "B" in (r["gene_a"], r["gene_b"])]


@pytest.mark.skipif(not SRC.exists() or shutil.which("minimap2") is None,
                    reason="needs sources/BK012101.1 and minimap2 on PATH")
def test_real_derive(paf, tmp_path):
    tgt = ROOT / "sources" / "KT899744delRR" / "KT899744delRR.fa"
    lifted, master = tmp_path / "lift.gtf", tmp_path / "master.gtf"
    subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "lift_bk012101.py"),
         "--gff", str(SRC / "BK012101.1.gff3"), "--src-fa", str(SRC / "BK012101.1.fa"),
         "--tgt-fa", str(tgt), "--paf", str(paf), "--out-gtf", str(lifted),
         "--report", str(tmp_path / "rep.tsv"), "--date", "2026.11.10",
         "--keep-noncoding"], check=True, capture_output=True)
    subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "build_master.py"), "--in-gtf", str(lifted),
         "--out-gtf", str(master), "--version", "v0", "--fasta", str(tgt)], check=True, capture_output=True)
    names = {parse_attrs(r[8])["gene_id"] for r in read_gtf(str(master)) if r[2] in ("TSS", "TES")}
    assert not names & {P + "TSS_1p", P + "TSS_135161m", P + "TES_135161p"}
    o = tmp_path / "out"
    subprocess.run([sys.executable, str(SCRIPT), "--master", str(master), "--outdir", str(o),
                    "--name", NAME], check=True, capture_output=True)
    s = {r["metric"]: int(r["value"]) for r in tsv(o / f"{NAME}.annotation_summary.tsv")}
    assert (s["total_cds"], s["embedded"], s["excluded"], s["retained"],
            s["retained_zero_unique"]) == (258, 40, 40, 218, 0)
    assert s["tss"] > 0 and s["tes"] > 0
    genome = read_fasta(str(tgt))
    for kind in ("quantification", "visualization"):
        assert check(read_gtf(str(o / f"{NAME}.{kind}.gtf")), genome) == []
