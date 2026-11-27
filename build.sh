#!/usr/bin/env bash
# Build the KT899744delRR annotation from the pinned sources.
# Usage: ./build.sh [--outdir output/KT899744delRR]
# Requires minimap2 on PATH:  module load gcc/14.2.0 ; module load minimap2/2.30
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-/project/dremel_lab/miniforge3/envs/pipelines/bin/python}"
HAROLD_SCRIPTS="${HAROLD_SCRIPTS:-/project/dremel_lab/workflows/pipelines/HAROLD/dev/workflow/scripts}"
RIBOLD_SCRIPTS="${RIBOLD_SCRIPTS:-/sfs/gpfs/tardis/project/dremel_lab/workflows/pipelines/RIBOLD/ribold/workflow/scripts}"
NAME="KT899744delRR"
OUTDIR="output/$NAME"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="${2:?--outdir needs a value}"; shift 2 ;;
    -h|--help) sed -n 2,4p "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

if ! command -v minimap2 >/dev/null 2>&1; then
  echo "error: minimap2 not found on PATH. Run:" >&2
  echo "  module load gcc/14.2.0" >&2
  echo "  module load minimap2/2.30" >&2
  exit 1
fi

cd "$ROOT"
[[ "$OUTDIR" = /* ]] || OUTDIR="$ROOT/$OUTDIR"
VERSION="$(tr -d '[:space:]' < VERSION)"
BK=sources/BK012101.1
DELRR=sources/$NAME/$NAME.fa
WORK="$ROOT/work"
mkdir -p "$OUTDIR" "$WORK"

echo "== 1. sources"
if ! md5sum --quiet -c sources/checksums.md5 >/dev/null 2>&1; then
  scripts/fetch_sources.sh "$ROOT"
fi
md5sum -c sources/checksums.md5

echo "== 2. verify delRR genome"
"$PYTHON" scripts/verify_delrr.py --kt sources/KT899744.1/KT899744.1.fa --delrr "$DELRR"

echo "== 3. align BK012101.1 to $NAME"
minimap2 -cx asm20 --cs --eqx -N 50 -p 0.3 "$DELRR" "$BK/BK012101.1.fa" \
  > "$OUTDIR/bk_to_delrr.paf" 2> "$WORK/minimap2.log"

echo "== 4. lift annotation"
"$PYTHON" scripts/lift_bk012101.py --gff "$BK/BK012101.1.gff3" --src-fa "$BK/BK012101.1.fa" \
  --tgt-fa "$DELRR" --paf "$OUTDIR/bk_to_delrr.paf" --out-gtf "$WORK/lifted.gtf" \
  --report "$OUTDIR/liftover_report.tsv" --date "$VERSION" --keep-noncoding

echo "== 5. master annotation"
"$PYTHON" scripts/build_master.py --in-gtf "$WORK/lifted.gtf" \
  --out-gtf "$OUTDIR/$NAME.master.gtf" --version "$VERSION" --fasta "$DELRR"

echo "== 6. derive quantification and visualization GTFs"
"$PYTHON" scripts/derive_gtfs.py --master "$OUTDIR/$NAME.master.gtf" --outdir "$OUTDIR" --name "$NAME"

echo "== 7. check GTFs"
for kind in quantification visualization; do
  "$PYTHON" scripts/check_gtf.py --gtf "$OUTDIR/$NAME.$kind.gtf" --fasta "$DELRR"
done

echo "== 8. pipeline-helper check"
printf 'KT899744delRR.BK012101.quant\tKT899744\n' > "sources/$NAME/$NAME.fa.regions"
QUANT="$OUTDIR/$NAME.quantification.gtf"
# _fix_gtf.py writes gene_id_2_gene_name.tsv into the current directory
(cd "$WORK" \
  && "$PYTHON" "$HAROLD_SCRIPTS/_fix_gtf.py" --ingtf "$QUANT" --outgtf "$WORK/fixed.gtf" \
  && "$PYTHON" "$HAROLD_SCRIPTS/_build_cds_gtf.py" --ingtf "$WORK/fixed.gtf" \
       --regions "$ROOT/sources/$NAME/$NAME.fa.regions" --host "" --outgtf "$WORK/cds.gtf" \
  && "$PYTHON" "$RIBOLD_SCRIPTS/_fix_gtf.py" --ingtf "$QUANT" --outgtf "$WORK/fixed_ribold.gtf")
ncds() { awk -F'\t' '$3=="CDS"' "$1" | wc -l; }
N_Q="$(ncds "$QUANT")"; N_F="$(ncds "$WORK/fixed.gtf")"; N_C="$(ncds "$WORK/cds.gtf")"
N_R="$(ncds "$WORK/fixed_ribold.gtf")"
echo "CDS rows: quant=$N_Q harold_fix=$N_F build_cds=$N_C ribold_fix=$N_R"
if [[ "$N_Q" != "$N_F" || "$N_F" != "$N_C" || "$N_Q" != "$N_R" ]]; then
  echo "error: CDS row counts differ across pipeline helpers" >&2
  exit 1
fi

echo "== 9. summary arithmetic"
SUM="$OUTDIR/$NAME.annotation_summary.tsv"
"$PYTHON" - "$SUM" <<'PY'
import csv, sys
s = {r["metric"]: int(r["value"]) for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
assert s["total_cds"] == (s["non_overlapping"] + s["partially_overlapping"]
                          + s["embedded"] + s["retained_zero_unique"]), s
assert s["excluded"] == s["embedded"], s
assert s["retained"] == s["total_cds"] - s["excluded"], s
assert s["tss"] > 0 and s["tes"] > 0, s
if s["retained_zero_unique"] > 0:
    print("WARNING: %d retained_zero_unique ORF(s); expect zero counts for them" % s["retained_zero_unique"])
print("summary arithmetic OK")
PY
cat "$SUM"
echo "build complete: $OUTDIR"
