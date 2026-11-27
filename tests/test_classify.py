import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from classify import classify, footprints
from gtflib import read_gtf

P = "KT899744delRR_"
ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "sources" / "BK012101.1"


def fp(*items):
    """items: (name, strand, [(start, end), ...]) -> footprint dict."""
    return {P + n: (s, frozenset(b for a, e in iv for b in range(a, e + 1)))
            for n, s, iv in items}


def test_non_overlapping():
    r = classify(fp(("A", "+", [(1, 10)]), ("B", "+", [(20, 30)])))
    assert r[P + "A"].cls == "non_overlapping" and r[P + "A"].unique_len == 10


def test_partial_overlap_kept():
    r = classify(fp(("A", "+", [(1, 20)]), ("B", "+", [(10, 40)])))
    assert r[P + "A"].cls == r[P + "B"].cls == "partially_overlapping"
    assert r[P + "A"].unique_len == 9


def test_embedded_single_host():
    r = classify(fp(("A", "+", [(10, 20)]), ("B", "+", [(1, 100)])))
    a = r[P + "A"]
    assert (a.cls, a.host, a.reason) == ("embedded", P + "B", "embedded_in_single_orf")
    assert r[P + "B"].cls == "non_overlapping"


def test_opposite_strand_not_embedded():
    r = classify(fp(("A", "+", [(10, 20)]), ("B", "-", [(1, 100)])))
    assert r[P + "A"].cls == "non_overlapping"


def test_union_cover_not_embedded():
    r = classify(fp(("A", "+", [(40, 60)]), ("B", "+", [(1, 50)]), ("C", "+", [(50, 100)])))
    assert r[P + "A"].cls == "retained_zero_unique"
    assert r[P + "A"].unique_len == 0


def test_identical_footprint_tiebreak():
    r = classify(fp(("UL26", "+", [(1, 30)]), ("UL26_iORF", "+", [(1, 30)])))
    assert r[P + "UL26"].cls == "non_overlapping"
    x = r[P + "UL26_iORF"]
    assert (x.cls, x.reason, x.host) == ("embedded", "identical_footprint", P + "UL26")


def test_identical_footprint_no_suffix_variants():
    # neither plain: longer annotated CDS wins, then alphabetical
    f = fp(("A_x", "+", [(1, 30)]), ("A_y", "+", [(1, 30)]))
    r = classify(f)
    assert r[P + "A_x"].cls == "non_overlapping"
    assert r[P + "A_y"].reason == "identical_footprint"
    # annotated length overrides alphabetical order
    r = classify(f, {P + "A_x": 30, P + "A_y": 60})
    assert r[P + "A_y"].cls == "non_overlapping"
    assert (r[P + "A_x"].cls, r[P + "A_x"].reason, r[P + "A_x"].host) == (
        "embedded", "identical_footprint", P + "A_y")


def test_chain():
    r = classify(fp(("A", "+", [(10, 20)]), ("B", "+", [(5, 50)]), ("C", "+", [(1, 100)])))
    assert r[P + "A"].cls == r[P + "B"].cls == "embedded"
    assert r[P + "A"].host == P + "B"  # smallest containing footprint
    assert r[P + "C"].cls == "non_overlapping"


def test_variant_transcripts_same_gene():
    rows = [["KT899744", "Genbank", "CDS", str(s), str(e), ".", "+", "0",
             f'gene_id "{P}G"; transcript_id "{t}";']
            for t, s, e in (("t1", 1, 100), ("t2", 10, 20))]
    f = footprints(rows)
    assert f[P + "G"] == ("+", frozenset(range(1, 101)))
    assert classify(f)[P + "G"].cls == "non_overlapping"


def test_unique_len_after_exclusion():
    r = classify(fp(("A", "+", [(10, 20)]), ("B", "+", [(1, 100)])))
    assert r[P + "B"].unique_len == 100
    assert r[P + "A"].unique_len == 0 and r[P + "A"].length == 11


EXPECTED_EMBEDDED = """ORF-P, UL5 uORF, UL5.5 dORF, UL5.5 iORF 1, UL5.5 iORF 2, UL5.5 iORF 3, UL5.6, UL5.6 dORF, UL6 iORF, UL8 iORF, UL8.5, UL9.5 iORF 1, UL9.5 iORF 2, UL12.5, UL15 iORF 3, UL15.5, UL17 uORF RNA *1, UL20 iORF, UL20 uORF, UL23 uORF 1, UL23 uORF 2, UL24.5, UL26.5, UL29.5, UL30.5 sORF, UL32.5, UL32.5 uORF, UL39.4 sORF, UL40.5 iORF, UL41 iORF, UL41 uORF 2, UL43.5 iORF, UL44 iORF, UL49.5 uORF, UL53 iORF, UL53 uORF, US1.5, US3.5, US3.6 sORF, US3.6 uORF"""


@pytest.mark.skipif(not SRC.exists() or shutil.which("minimap2") is None,
                    reason="needs sources/BK012101.1 and minimap2 on PATH")
def test_real_master_counts(paf, tmp_path):
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
         "--out-gtf", str(master), "--version", "v0"], check=True, capture_output=True)
    res = classify(footprints(read_gtf(str(master))))
    emb = {g[len(P):] for g, c in res.items() if c.cls == "embedded"}
    expected = {n.strip().replace(" ", "_") for n in EXPECTED_EMBEDDED.split(",")}
    assert len(res) == 258
    assert sorted(expected - emb) == [] and sorted(emb - expected) == []
    assert len(emb) == 40
    assert sum(c.cls != "embedded" for c in res.values()) == 218
    assert not [g for g, c in res.items() if c.reason == "identical_footprint"]
    assert not [g for g, c in res.items() if c.cls == "retained_zero_unique"]
