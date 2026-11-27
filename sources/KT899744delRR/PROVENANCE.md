# KT899744delRR provenance

There is no publication for KT899744delRR. Aligning it back to KT899744.1
(Colgrove et al. 2015, Virology 487:215-221) shows it is KT899744.1 with three
blocks removed and no other change. Keeping KT899744.1 bases 9064-125844 and
126978-145360 reproduces `KT899744delRR.fa` exactly (135,164 bp, zero
mismatches).

| Deleted from KT899744.1 | Length | Interpretation |
|---|---|---|
| 1-9063 | 9,063 bp | TRL; inverted copy retained from 126,012 (IRL) |
| 125845-126977 | 1,133 bp | Extra copy of terminal sequence at the IRL/IRS junction; its start matches, inverted, the sequence from base 168. Consistent with the "a" sequence; not confirmed. |
| 145361-151974 | 6,614 bp | TRS; inverted copy retained from 132,227 (IRS) |

## Coordinate conversion

- delRR position = KT899744.1 position - 9063, for positions 9064-125844
- delRR position = KT899744.1 position - 10196, for positions 126978-145360

## Header typo

The FASTA header reads `Δ1-9603`; the sequence shows the deletion is 1-9063.
The shared reference file is not edited; pipelines read only the first word
(`KT899744`).

## Verification

Download the sources, then check the rebuild:

```bash
scripts/fetch_sources.sh .
python scripts/verify_delrr.py \
  --kt sources/KT899744.1/KT899744.1.fa \
  --delrr sources/KT899744delRR/KT899744delRR.fa
```

The script exits 0 if the rebuilt sequence is identical, 1 otherwise.
