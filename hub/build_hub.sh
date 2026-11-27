#!/usr/bin/env bash
# Build a UCSC assembly hub for KT899744delRR from the committed outputs.
# Tracks: quantification GTF, visualization GTF, TSS/TES points, excluded ORFs, miRNAs.
# Needs UCSC tools (faToTwoBit, twoBitInfo, gtfToGenePred, genePredToBed,
# bedToBigBed); on Rivanna they are in the biocondatools env (BIN below).
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN=${UCSC_BIN:-/project/dremel_lab/miniforge3/envs/biocondatools/bin}
OUT=$ROOT/output/KT899744delRR
FA=$ROOT/sources/KT899744delRR/KT899744delRR.fa
HUB=$ROOT/hub/KT899744delRR
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

for t in faToTwoBit twoBitInfo gtfToGenePred genePredToBed bedToBigBed; do
  [[ -x $BIN/$t ]] || { echo "missing $BIN/$t (set UCSC_BIN)" >&2; exit 1; }
done
mkdir -p "$HUB"

"$BIN/faToTwoBit" "$FA" "$HUB/KT899744delRR.2bit"
"$BIN/twoBitInfo" "$HUB/KT899744delRR.2bit" "$HUB/chrom.sizes"

# GTF -> bigBed 12, one item per transcript, labelled by gene name
# (prefix KT899744delRR_ dropped). Only transcript structure rows are used.
gtf_to_bb() {
  local gtf=$1 out=$2
  awk -F'\t' '$3=="exon"||$3=="CDS"||$3=="start_codon"||$3=="stop_codon"' "$gtf" > "$TMP/in.gtf"
  "$BIN/gtfToGenePred" -genePredExt "$TMP/in.gtf" "$TMP/gp"
  # genePredExt col 12 = gene name; use it as the item name
  awk -F'\t' 'BEGIN{OFS="\t"} {n=$12; sub(/^KT899744delRR_/, "", n); $1=n; print}' "$TMP/gp" > "$TMP/gp.named"
  "$BIN/genePredToBed" "$TMP/gp.named" stdout | LC_ALL=C sort -k1,1 -k2,2n -k4,4 > "$TMP/bed"
  "$BIN/bedToBigBed" -type=bed12 "$TMP/bed" "$HUB/chrom.sizes" "$out"
}
gtf_to_bb "$OUT/KT899744delRR.quantification.gtf" "$HUB/quantification.bb"
gtf_to_bb "$OUT/KT899744delRR.visualization.gtf" "$HUB/visualization.bb"

# TSS (green) and TES (red) points, bed9 with itemRgb
awk -F'\t' 'BEGIN{OFS="\t"} $3=="TSS"||$3=="TES" {
  c = ($3=="TSS") ? "0,150,0" : "200,0,0"
  print $1, $4-1, $5, $3, 0, $7, $4-1, $5, c }' "$OUT/KT899744delRR.visualization.gtf" \
  | LC_ALL=C sort -k1,1 -k2,2n -k4,4 -k6,6 > "$TMP/tss_tes.bed"
"$BIN/bedToBigBed" -type=bed9 "$TMP/tss_tes.bed" "$HUB/chrom.sizes" "$HUB/tss_tes.bb"

# ORFs excluded from quantification (CDS span), named "ORF in HOST"
awk -F'\t' 'BEGIN{OFS="\t"} NR>1 {
  n=$1; h=$8; sub(/^KT899744delRR_/, "", n); sub(/^KT899744delRR_/, "", h)
  print $3, $4-1, $5, n "_in_" h, 0, $6 }' "$OUT/KT899744delRR.excluded_CDSs.tsv" \
  | LC_ALL=C sort -k1,1 -k2,2n -k4,4 > "$TMP/excluded.bed"
"$BIN/bedToBigBed" -type=bed6 "$TMP/excluded.bed" "$HUB/chrom.sizes" "$HUB/excluded.bb"

# miRNA precursors and mature arms (BED6 written by scripts/place_mirnas.py)
LC_ALL=C sort -k1,1 -k2,2n -k4,4 -k6,6 "$OUT/KT899744delRR.miRNA.bed" > "$TMP/mirna.bed"
"$BIN/bedToBigBed" -type=bed6 "$TMP/mirna.bed" "$HUB/chrom.sizes" "$HUB/mirna.bb"

echo "hub files written to $HUB"
