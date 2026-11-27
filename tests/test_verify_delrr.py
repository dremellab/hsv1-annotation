from pathlib import Path

import pytest

from gtflib import format_attrs, parse_attrs, read_fasta
from verify_delrr import KEEP, rebuild

ROOT = Path(__file__).resolve().parent.parent
KT = ROOT / "sources" / "KT899744.1" / "KT899744.1.fa"
DELRR = ROOT / "sources" / "KT899744delRR" / "KT899744delRR.fa"


def test_rebuild_matches():
    kt = "ACGTACGTAC" "GGGGCCCCAA" "TTTTAAAAGG"  # 30 bp
    keep = [(3, 8), (13, 18), (25, 30)]
    expected = kt[2:8] + kt[12:18] + kt[24:30]
    assert rebuild(kt, keep) == expected
    assert len(expected) == 18


def test_real_delrr():
    if not KT.exists():
        pytest.skip("KT899744.1 not downloaded (run scripts/fetch_sources.sh)")
    out = rebuild(read_fasta(str(KT)), KEEP)
    assert out == read_fasta(str(DELRR))
    assert len(out) == 135164


def test_parse_format_roundtrip():
    s = 'gene_id "a"; x "b_c";'
    assert format_attrs(parse_attrs(s)) == s
