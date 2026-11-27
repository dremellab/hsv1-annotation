#!/usr/bin/env bash
# Download source records from NCBI and miRBase into <outdir>/sources and verify checksums.
# The miRBase URLs are unversioned and serve the current release. The pinned files
# are release 23.1, so once miRBase publishes a newer release the md5 check below
# fails by design; set MIRBASE_DIR to a directory holding 23.1 hairpin.fa and mature.fa.
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

# miRBase 23.1: keep the HSV-1 entries (header starts ">hsv1-") of the full files.
# MIRBASE_DIR (optional) points at already-downloaded hairpin.fa and mature.fa.
MB="$OUT/sources/miRBase_23.1"
mkdir -p "$MB"
MBTMP="$(mktemp -d)"
trap 'rm -rf "$MBTMP"' EXIT
for kind in hairpin mature; do
  if [[ -n "${MIRBASE_DIR:-}" ]]; then
    cp "$MIRBASE_DIR/$kind.fa" "$MBTMP/$kind.fa"
  else
    curl -sS -o "$MBTMP/$kind.fa" "https://www.mirbase.org/download/$kind.fa"
  fi
  awk '/^>/ {keep = ($1 ~ /^>hsv1-/)} keep' "$MBTMP/$kind.fa" > "$MB/hsv1_$kind.fa"
done

(cd "$OUT" && md5sum -c sources/checksums.md5)
