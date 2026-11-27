#!/usr/bin/env bash
# Copy the quantification and visualization GTFs into the shared reference folder.
# Usage: ./publish.sh [--apply] [--dest DIR]
# Without --apply it only prints what it would do and changes nothing.
# The existing KT899744delRR.BK012101.gtf is never touched.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NAME="KT899744delRR"
REF="$NAME.BK012101"
DEST="/project/dremel_lab/workflows/reference_data/fasta_gtf"
APPLY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --apply) APPLY=1; shift ;;
    --dest) DEST="${2:?--dest needs a value}"; shift 2 ;;
    -h|--help) sed -n 2,5p "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

OUT="$ROOT/output/$NAME"
SRC_FA="$ROOT/sources/$NAME/$NAME.fa"
QUANT_SRC="$OUT/$NAME.quantification.gtf"
VIS_SRC="$OUT/$NAME.visualization.gtf"
QUANT="$DEST/$REF.quant.gtf"
VIS="$DEST/$REF.vis.gtf"

for f in "$SRC_FA" "$QUANT_SRC" "$VIS_SRC"; do
  [[ -f "$f" ]] || { echo "error: missing $f (run ./build.sh first)" >&2; exit 1; }
done
[[ -d "$DEST" ]] || { echo "error: destination $DEST does not exist" >&2; exit 1; }
[[ -f "$DEST/$NAME.fa" && -f "$DEST/$NAME.fa.fai" ]] \
  || { echo "error: $DEST/$NAME.fa or .fai not found" >&2; exit 1; }

want="$(md5sum < "$SRC_FA" | cut -d' ' -f1)"
have="$(md5sum < "$DEST/$NAME.fa" | cut -d' ' -f1)"
if [[ "$want" != "$have" ]]; then
  echo "error: $DEST/$NAME.fa (md5 $have) differs from sources/$NAME/$NAME.fa (md5 $want); aborting" >&2
  exit 1
fi
echo "FASTA md5 matches: $want"

if [[ "$APPLY" -ne 1 ]]; then
  echo "dry run (no changes). With --apply it would, in $DEST:"
  echo "  copy  $QUANT_SRC -> $REF.quant.gtf"
  echo "  copy  $VIS_SRC -> $REF.vis.gtf"
  echo "  link  $REF.quant.fa -> $NAME.fa"
  echo "  link  $REF.quant.fa.fai -> $NAME.fa.fai"
  echo "  write $REF.quant.fa.regions"
  echo "  leave $REF.gtf unchanged"
  exit 0
fi

cp "$QUANT_SRC" "$QUANT"
cp "$VIS_SRC" "$VIS"
ln -sfn "$NAME.fa" "$DEST/$REF.quant.fa"
ln -sfn "$NAME.fa.fai" "$DEST/$REF.quant.fa.fai"
printf '%s.quant\tKT899744\n' "$REF" > "$DEST/$REF.quant.fa.regions"
echo "published to $DEST:"
ls -l "$DEST"/"$REF".quant.* "$DEST/$REF.vis.gtf"
