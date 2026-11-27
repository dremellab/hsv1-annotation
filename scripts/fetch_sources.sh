#!/usr/bin/env bash
# Download source records from NCBI into <outdir>/sources and verify checksums.
# Usage: scripts/fetch_sources.sh <outdir>   (outdir = repository root)
set -euo pipefail

OUT="${1:?usage: fetch_sources.sh <outdir>}"
EU="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

mkdir -p "$OUT/sources/BK012101.1" "$OUT/sources/KT899744.1"
BK="$OUT/sources/BK012101.1"
curl -sS -o "$BK/BK012101.1.gb" "$EU?db=nuccore&id=BK012101.1&rettype=gbwithparts&retmode=text"
curl -sS -o "$BK/BK012101.1.fa" "$EU?db=nuccore&id=BK012101.1&rettype=fasta&retmode=text"
curl -sS -o "$BK/BK012101.1.gff3" \
  "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=BK012101.1"
curl -sS -o "$OUT/sources/KT899744.1/KT899744.1.fa" \
  "$EU?db=nuccore&id=KT899744.1&rettype=fasta&retmode=text"

(cd "$OUT" && md5sum -c sources/checksums.md5)
