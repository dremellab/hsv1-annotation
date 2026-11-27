import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

_ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="session")
def paf(tmp_path_factory):
    src = _ROOT / "sources" / "BK012101.1" / "BK012101.1.fa"
    tgt = _ROOT / "sources" / "KT899744delRR" / "KT899744delRR.fa"
    out = tmp_path_factory.mktemp("paf") / "bk_to_delrr.paf"
    with open(out, "w") as fh:
        subprocess.run(
            ["minimap2", "-cx", "asm20", "--cs", "--eqx", "-N", "50", "-p", "0.3",
             str(tgt), str(src)],
            stdout=fh, stderr=subprocess.DEVNULL, check=True)
    return out
