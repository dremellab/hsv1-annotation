import random
import subprocess
import sys
from pathlib import Path

from build_master import tss_tes_rows
from gtflib import parse_attrs, read_gtf

BUILD = Path(__file__).resolve().parent.parent / "scripts" / "build_master.py"


def exon(tid, s, e, strand):
    return ["KT899744", "Genbank", "exon", str(s), str(e), ".", strand, ".",
            f'gene_id "g_{tid}"; transcript_id "{tid}";']


def by_feature(rows, feat):
    return [r for r in rows if r[2] == feat]


def test_plus_strand_ends():
    rows = tss_tes_rows([exon("t1", 100, 200, "+"), exon("t1", 300, 400, "+")])
    tss, tes = by_feature(rows, "TSS")[0], by_feature(rows, "TES")[0]
    assert (tss[3], tss[4]) == ("100", "100")
    assert (tes[3], tes[4]) == ("400", "400")
    assert tss[0] == "KT899744" and tss[1] == "Genbank" and tss[5] == "." and tss[7] == "."
    a = parse_attrs(tss[8])
    assert a["gene_id"] == a["transcript_id"] == a["gene_name"] == "KT899744delRR_TSS_100p"
    assert a["source_transcripts"] == "t1"
    assert parse_attrs(tes[8])["gene_id"] == "KT899744delRR_TES_400p"


def test_minus_strand_ends():
    rows = tss_tes_rows([exon("t1", 100, 200, "-"), exon("t1", 300, 400, "-")])
    tss, tes = by_feature(rows, "TSS")[0], by_feature(rows, "TES")[0]
    assert tss[3] == "400" and tes[3] == "100"
    assert parse_attrs(tss[8])["gene_id"] == "KT899744delRR_TSS_400m"
    assert parse_attrs(tes[8])["gene_id"] == "KT899744delRR_TES_100m"


def test_collapse_shared_tss():
    rows = tss_tes_rows([exon("t2", 100, 200, "+"), exon("t1", 100, 300, "+")])
    tss = by_feature(rows, "TSS")
    assert len(tss) == 1
    assert parse_attrs(tss[0][8])["source_transcripts"] == "t1,t2"
    assert len(by_feature(rows, "TES")) == 2


def test_deterministic_order(tmp_path):
    rows = [exon("t1", 100, 200, "+"), exon("t1", 300, 400, "+"),
            exon("t2", 100, 250, "-"), exon("t3", 100, 200, "+")]
    rows.append(["KT899744", "Genbank", "gene", "100", "400", ".", "+", ".",
                 'gene_id "g_t1";'])
    fa = tmp_path / "g.fa"
    fa.write_text(">g\n" + "A" * 1000 + "\n")
    outs = []
    for i in range(3):
        shuffled = rows[:]
        random.Random(i).shuffle(shuffled)
        inp = tmp_path / f"in{i}.gtf"
        inp.write_text("".join("\t".join(r) + "\n" for r in shuffled))
        out = tmp_path / f"out{i}.gtf"
        subprocess.run([sys.executable, str(BUILD), "--in-gtf", str(inp),
                        "--out-gtf", str(out), "--version", "v1.2.3",
                        "--fasta", str(fa)], check=True)
        outs.append(out.read_bytes())
    assert outs[0] == outs[1] == outs[2]
    assert outs[0].startswith(b"#!annotation-version v1.2.3\n")
    body = read_gtf(str(tmp_path / "out0.gtf"))
    assert body == sorted(body, key=lambda r: (int(r[3]), -int(r[4])))
    at100 = [r[2] for r in body if r[3] == "100" and r[6] == "+"]
    assert at100.index("gene") < at100.index("TSS")
    assert at100.index("gene") < at100.index("exon")
    assert len(body) == len(rows) + 5


def transcript(tid, strand, status):
    return ["KT899744", "Genbank", "transcript", "1", "1", ".", strand, ".",
            f'gene_id "g_{tid}"; transcript_id "{tid}"; liftover_status "{status}";']


def test_clipped_ends_at_genome_boundary_skipped():
    rows = [transcript("p", "+", "clipped_terminal_repeat"), exon("p", 1, 500, "+"),
            transcript("m", "-", "clipped_terminal_repeat;indel(-36)"), exon("m", 600, 997, "-"),
            transcript("q", "+", "clipped_terminal_repeat"), exon("q", 700, 997, "+")]
    out = tss_tes_rows(rows, 1000)
    tss = {r[3] for r in by_feature(out, "TSS")}
    tes = {r[3] for r in by_feature(out, "TES")}
    assert tss == {"700"} and tes == {"500", "600"}  # p TSS@1, m TSS@997, q TES@997 skipped (alignment ends 3 bp short of genome end)


def test_unclipped_or_interior_ends_kept():
    rows = [transcript("p", "+", "exact"), exon("p", 1, 1000, "+"),
            transcript("q", "+", "clipped_terminal_repeat"), exon("q", 10, 900, "+")]
    out = tss_tes_rows(rows, 1000)
    assert {r[3] for r in by_feature(out, "TSS")} == {"1", "10"}
    assert {r[3] for r in by_feature(out, "TES")} == {"1000", "900"}


def test_cli_mirna_gtf_merged_and_sorted(tmp_path):
    base = tmp_path / "in.gtf"
    mi = tmp_path / "mi.gtf"
    (tmp_path / "g.fa").write_text(">KT899744\n" + "A" * 1000 + "\n")
    base.write_text("\t".join(exon("t1", 100, 200, "+")) + "\n")
    mi.write_text("#x\nKT899744\tmiRBase_23.1\tgene\t50\t70\t.\t+\t.\tgene_id \"hsv1-mir-H1\";\n")
    out = tmp_path / "out.gtf"
    subprocess.run([sys.executable, str(BUILD), "--in-gtf", str(base), "--mirna-gtf", str(mi),
                    "--out-gtf", str(out), "--version", "v0", "--fasta", str(tmp_path / "g.fa")],
                   check=True, capture_output=True)
    body = read_gtf(str(out))
    assert body[0][8].startswith('gene_id "hsv1-mir-H1"') and body[0][3] == "50"
