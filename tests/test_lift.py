import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from check_gtf import check
from gtflib import read_fasta, read_gtf

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "sources" / "BK012101.1"
TGT = ROOT / "sources" / "KT899744delRR" / "KT899744delRR.fa"
EXPECTED = Path("/project/dremel_lab/workflows/reference_data/fasta_gtf/KT899744delRR.BK012101.gtf")
LIFT = ROOT / "scripts" / "lift_bk012101.py"

pytestmark = pytest.mark.skipif(
    not SRC.exists() or shutil.which("minimap2") is None,
    reason="needs sources/BK012101.1 and minimap2 on PATH",
)


@pytest.fixture(scope="module")
def paf(tmp_path_factory):
    out = tmp_path_factory.mktemp("paf") / "bk_to_delrr.paf"
    with open(out, "w") as fh:
        subprocess.run(
            ["minimap2", "-cx", "asm20", "--cs", "--eqx", "-N", "50", "-p", "0.3",
             str(TGT), str(SRC / "BK012101.1.fa")],
            stdout=fh, stderr=subprocess.DEVNULL, check=True)
    return out


def run_lift(paf, tmp_path, *extra):
    out = tmp_path / "out.gtf"
    subprocess.run(
        [sys.executable, str(LIFT), "--gff", str(SRC / "BK012101.1.gff3"),
         "--src-fa", str(SRC / "BK012101.1.fa"), "--tgt-fa", str(TGT),
         "--paf", str(paf), "--out-gtf", str(out), "--report", str(tmp_path / "rep.tsv"),
         "--date", "2026.11.10", *extra],
        check=True, capture_output=True)
    return out


def body(path):
    return [l for l in Path(path).read_text().splitlines() if l and not l.startswith("#")]


def test_regression_without_flag(paf, tmp_path):
    out = run_lift(paf, tmp_path)
    got = [l.replace(' transcript_biotype "protein_coding";', "") for l in body(out)]
    assert got == body(EXPECTED)


def test_keep_noncoding_adds_18(paf, tmp_path):
    rows = read_gtf(str(run_lift(paf, tmp_path, "--keep-noncoding")))
    nc = [r for r in rows if r[2] == "transcript" and 'transcript_biotype "noncoding"' in r[8]]
    assert len(nc) == 18
    assert sum(1 for r in rows if r[2] == "gene") == 258 + 18


def test_checker_passes(paf, tmp_path):
    out = run_lift(paf, tmp_path, "--keep-noncoding")
    assert check(read_gtf(str(out)), read_fasta(str(TGT))) == []
