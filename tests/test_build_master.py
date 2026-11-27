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
    outs = []
    for i in range(3):
        shuffled = rows[:]
        random.Random(i).shuffle(shuffled)
        inp = tmp_path / f"in{i}.gtf"
        inp.write_text("".join("\t".join(r) + "\n" for r in shuffled))
        out = tmp_path / f"out{i}.gtf"
        subprocess.run([sys.executable, str(BUILD), "--in-gtf", str(inp),
                        "--out-gtf", str(out), "--version", "v1.2.3"], check=True)
        outs.append(out.read_bytes())
    assert outs[0] == outs[1] == outs[2]
    assert outs[0].startswith(b"#!annotation-version v1.2.3\n")
    body = read_gtf(str(tmp_path / "out0.gtf"))
    assert body == sorted(body, key=lambda r: (int(r[3]), int(r[4])))
    assert len(body) == len(rows) + 5
