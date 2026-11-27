#!/usr/bin/env python3
"""Check that KT899744delRR is KT899744.1 with three blocks removed."""
from __future__ import annotations

import argparse
import sys

from gtflib import read_fasta

# 1-based, inclusive KT899744.1 blocks that are kept.
KEEP = [(9064, 125844), (126978, 145360)]


def rebuild(kt: str, keep: list[tuple[int, int]]) -> str:
    return "".join(kt[a - 1 : b] for a, b in keep)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--kt", required=True, help="KT899744.1.fa")
    ap.add_argument("--delrr", required=True, help="KT899744delRR.fa")
    args = ap.parse_args()
    built = rebuild(read_fasta(args.kt), KEEP)
    ref = read_fasta(args.delrr)
    if built == ref:
        print(f"OK: identical ({len(ref)} bp)")
        return 0
    n = min(len(built), len(ref))
    mism = sum(1 for i in range(n) if built[i] != ref[i])
    print(f"FAIL: rebuilt {len(built)} bp vs {len(ref)} bp, {mism} mismatches")
    return 1


if __name__ == "__main__":
    sys.exit(main())
