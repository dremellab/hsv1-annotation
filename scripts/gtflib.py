#!/usr/bin/env python3
"""Small helpers for FASTA and GTF handling (standard library only)."""
from __future__ import annotations

import re

_COMP = str.maketrans("ACGTUNacgtunRYKMBDHVryk mbdhv".replace(" ", ""),
                      "TGCAANtgcaanYRMKVHDByrmkvhdb")
_ATTR = re.compile(r'\s*([^\s"]+)\s+"([^"]*)"\s*;?')


def read_fasta(path: str) -> str:
    """Return all sequence lines concatenated, uppercase."""
    parts: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            parts.append(line.strip().upper())
    return "".join(parts)


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def parse_attrs(field: str) -> dict[str, str]:
    return {m.group(1): m.group(2) for m in _ATTR.finditer(field)}


def format_attrs(d: dict[str, str]) -> str:
    return " ".join(f'{k} "{v}";' for k, v in d.items())


def read_gtf(path: str) -> list[list[str]]:
    rows: list[list[str]] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            rows.append(line.split("\t"))
    return rows


def write_gtf(path: str, rows: list[list[str]], header: list[str]) -> None:
    with open(path, "w") as fh:
        for h in header:
            fh.write(h if h.startswith("#") else f"#{h}")
            fh.write("\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
