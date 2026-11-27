# hsv1-annotation

GTF annotations of the herpes simplex virus 1 (HSV-1) genome KT899744delRR,
built for read counting (quantification) and for genome browsers
(visualization). Version: see `VERSION`.

The annotation is the BK012101.1 annotation (Whisnant et al. 2020) lifted
onto KT899744delRR. One gene is defined per ORF, ORFs that lie entirely
inside another ORF are removed from the quantification GTF, and TSS and TES
features are added to the visualization GTF. Host-virus RNA-seq and Ribo-seq
pipelines (HAROLD, RIBOLD) can use the quantification GTF without code changes.

miRNAs are out of scope for this version and will be added to the
visualization GTF later.

## Quick start

Use the published quantification GTF in a pipeline config:

```yaml
viruses: "KT899744delRR.BK012101.quant"
```

Do not add `-O` or `--fraction` to featureCounts with this GTF (see
"Ambiguous reads").

Rebuild everything from the pinned sources:

```bash
module load gcc/14.2.0
module load minimap2/2.30
./build.sh
```

`./build.sh` is the entry point. It runs with the Python 3.11 of the shared
`pipelines` mamba env (override with `PYTHON=...`); the scripts use only the
standard library. minimap2 is loaded with the two `module load` commands above.

## Reference genome

KT899744delRR is KT899744.1 (HSV-1 KOS, Colgrove et al. 2015) with three
blocks removed and no other change; there is no publication for it. Keeping
KT899744.1 bases 9064-125844 and 126978-145360 reproduces
`sources/KT899744delRR/KT899744delRR.fa` exactly (135,164 bp, zero
mismatches). `build.sh` re-checks this on every run.

| Deleted from KT899744.1 | Length | Interpretation |
|---|---|---|
| 1-9063 | 9,063 bp | TRL; inverted copy retained from 126,012 (IRL) |
| 125845-126977 | 1,133 bp | Extra copy of terminal sequence at the IRL/IRS junction; consistent with the "a" sequence, not confirmed |
| 145361-151974 | 6,614 bp | TRS; inverted copy retained from 132,227 (IRS) |

Details are in `sources/KT899744delRR/PROVENANCE.md`. The shared FASTA is not
edited. Its header reads `Δ1-9603`, but the sequence shows the deletion is
1-9063; pipelines read only the first word of the header (`KT899744`).

## Method

1. **Fetch and verify.** BK012101.1 (GenBank record, FASTA, GFF3) and
   KT899744.1 (FASTA) are checked against `sources/checksums.md5`, and
   KT899744delRR is rebuilt from KT899744.1 and compared with the committed
   FASTA.
2. **Align.** `minimap2 -cx asm20 --cs --eqx -N 50 -p 0.3` of BK012101.1 onto
   KT899744delRR. Only the collinear block is used for lifting.
3. **Lift.** One gene per ORF, named by the paper's CDS product. Variant
   records (`*1`, `#1`, `iso1`) are folded into their base ORF. Copies in the
   terminal repeats are dropped, and transcripts crossing a deleted boundary
   are clipped. Each CDS is re-translated and flagged in `cds_check` (for
   example `identical_aa`, `aa_subs`, `frameshift`). The 18 transcripts
   without a CDS are kept as `transcript_biotype "noncoding"` under their own
   gene.
4. **Master.** One 1 bp `TSS` row (at the 5' end) and one 1 bp `TES` row (at
   the 3' end) per BK012101.1 mRNA, collapsing exact duplicates. Each row
   carries the `source_transcripts` it came from. A TSS or TES is omitted
   when its transcript was clipped at a genome boundary and that end lies at
   the clip point (position 1 or the end of the alignment at 135,161); the
   true end is outside the genome.
5. **Derive.** ORFs are classified (below), then the quantification GTF, the
   visualization GTF and the reports are written.
6. **Validate.** Coordinates in bounds, gene/transcript/exon structure,
   frames, stop codons, `gene_id` and `transcript_id` on every row, and no
   spaces in attribute values (`_fix_gtf.py` splits on whitespace). The
   quantification GTF is passed through HAROLD `_fix_gtf.py` and
   `_build_cds_gtf.py` and the CDS row counts must match. Summary arithmetic
   is checked. Outputs carry no timestamps (the date in the header comes from
   `VERSION`), so two builds give identical files.

### ORF rules

An ORF is a `gene_id` in the master annotation. Its CDS footprint is the union
of the genomic bases of all CDS rows in all its transcripts; stop codons are
excluded, as in the GTF. Only same-strand overlaps count. Each ORF gets one
class:

| Class | Definition | In quantification GTF |
|---|---|---|
| `non_overlapping` | No CDS base shared with any other same-strand ORF, kept or excluded | Kept |
| `partially_overlapping` | At least one CDS base shared with another same-strand ORF (including an excluded embedded ORF) and not itself embedded | Kept |
| `embedded` | All CDS bases inside one other ORF | Excluded; host in `containing_cds_id` |
| `retained_zero_unique` | Kept, but every base is shared with other kept ORFs | Kept and flagged |

Rules for special cases:
- If two ORFs have identical footprints, the one whose name has no suffix is
  kept, then the one with the longer annotated CDS, then the first
  alphabetically. The other is excluded with reason `identical_footprint`.
  The current build has none.
- Chains (A inside B inside C): all classes are decided against the original
  master in one pass, so A and B are both excluded.
- CDSs that are frameshifted in KOS relative to BK012101.1 are classified
  normally, flagged `cds_check "frameshift"`, and kept.

Result of the current build (`annotation_summary.tsv`): 258 ORFs, of which 143
are non-overlapping, 75 partially overlapping and 40 embedded. The 40 embedded
ORFs are excluded, leaving 218 in the quantification GTF; 0 are
`retained_zero_unique`. Each excluded ORF has exactly one host, listed in
`excluded_CDSs.tsv` (for example UL26.5 in UL26, UL12.5 in UL12); a host of an
excluded ORF is therefore itself `partially_overlapping` (UL26 is one). The
master and visualization GTFs have 290 TSS and 157 TES rows.

### Strand and ambiguous reads

Opposite-strand overlaps are listed in `CDS_overlaps.tsv` for information only
and do not affect classes. They do not cause lost reads only when counting is
strand-aware (stranded library, featureCounts `-s 1` or `-s 2`); with `-s 0`
they are also ambiguous. Both pipelines run featureCounts with
`-t CDS -g gene_id` and without `-O`, so a read overlapping more than one kept
ORF is discarded as ambiguous. Excluding embedded ORFs makes reads in them
count for the host ORF instead of being discarded. `-O` would assign reads in
overlapping ORFs to every ORF they touch (with `--fraction`, as split
fractional counts), which changes counts relative to earlier runs. Do not add
either when using this GTF.

### `retained_zero_unique`

An ORF is `retained_zero_unique` if it is kept but every one of its CDS bases
is shared with another kept same-strand ORF. Because ambiguous reads are
discarded, such an ORF is expected to get zero counts. The class exists to make
these cases visible; there are none in the current build.

### Effective unique length

For each kept ORF, the number of its CDS bases not shared with any other kept
same-strand ORF, computed on the quantification GTF. `CDS_overlaps.tsv` reports
it (`unique_len_a`, `unique_len_b`) next to the full CDS length. Neither
pipeline uses it. If you want length-normalized values that reflect only the
bases that can yield uniquely assigned reads, divide counts by this length
instead of the full CDS length (for example in TPM or RPKM). It is only an
approximation of the region that can yield uniquely assigned reads: a read that
starts in a unique region and runs into a shared one is still discarded, so the
usable region is shorter by about one read length at each boundary (this
matters for short ORFs and for Ribo-seq). ORFs with
effective length 0 are the `retained_zero_unique` ORFs and cannot be
normalized this way.

## Outputs

Committed in `output/KT899744delRR/`:

| File | Contents |
|---|---|
| `KT899744delRR.master.gtf` | 258 ORF genes plus 18 noncoding genes, all transcripts, TSS and TES rows |
| `KT899744delRR.quantification.gtf` | 218 ORFs, CDS-carrying transcripts only; no TSS, TES or noncoding transcripts |
| `KT899744delRR.visualization.gtf` | Same content as the master, sorted, attributes tidied for genome browsers |
| `KT899744delRR.excluded_CDSs.tsv` | One row per excluded ORF: CDS ID, gene name, chromosome, start, end, strand, CDS length, containing CDS ID, overlapping bases, unique bases, reason, `cds_check` |
| `KT899744delRR.CDS_overlaps.tsv` | One row per overlapping ORF pair (either strand): each ORF's class, full length and effective unique length |
| `KT899744delRR.annotation_summary.tsv` | Totals: CDSs, non-overlapping, partially overlapping, embedded, excluded, retained, `retained_zero_unique`, miRNA (0 until they are added), TSS, TES |
| `liftover_report.tsv` | Per BK012101.1 transcript: lift status, coordinates, protein ID, `cds_check`, amino-acid differences |
| `bk_to_delrr.paf` | The minimap2 alignment used for lifting |

All GTFs use sequence name `KT899744`, source `Genbank`, and gene IDs prefixed
`KT899744delRR_`.

## Publishing to the shared reference folder

Pipelines find references as `<name>.fa` and `<name>.gtf` in
`/project/dremel_lab/workflows/reference_data/fasta_gtf/`, selected through
`viruses:` in the config. `publish.sh` puts the new files there. Run it only
when asked; it writes outside the repository.

```bash
./publish.sh            # dry run: checks the FASTA and lists what it would do
./publish.sh --apply    # copy files into the shared folder
./publish.sh --apply --dest DIR   # use another folder (for testing)
```

It first checks that `KT899744delRR.fa` in the destination has the same MD5 as
`sources/KT899744delRR/KT899744delRR.fa` and aborts otherwise. Then it writes:

- `KT899744delRR.BK012101.quant.gtf` (quantification GTF)
- `KT899744delRR.BK012101.quant.fa` and `.fa.fai`, relative symlinks to
  `KT899744delRR.fa` and `KT899744delRR.fa.fai`
- `KT899744delRR.BK012101.quant.fa.regions` (`KT899744delRR.BK012101.quant`,
  tab, `KT899744`)
- `KT899744delRR.BK012101.vis.gtf` (visualization GTF)

### Why the original GTF is unchanged

`KT899744delRR.BK012101.gtf` was used for earlier RIBOLD and HAROLD runs.
Replacing it would change what those runs mean when repeated. It is left
as it is, and `publish.sh` never writes to it. New analyses choose the
`.quant` name explicitly.

## Visualization GTF

`KT899744delRR.BK012101.vis.gtf` is meant for IGV and UCSC: it keeps all
ORFs (including embedded ones such as UL26.5), the noncoding transcripts and
the TSS/TES features. Do not use it for counting. Checks in IGV and UCSC
have not been recorded yet; any feature types a browser rejects will be listed
here.

## Tests

```bash
/project/dremel_lab/miniforge3/envs/pipelines/bin/python -m pytest tests
```

The tests use small synthetic GTFs: each class, strand-aware overlap,
containment across several transcripts of one gene, the identical-footprint
tie-break, chains, and effective unique length after exclusion.

## Methods text for publication

A draft is in `docs/methods.md`.

## Citing

If you use these annotations, please cite the sources they derive from:

- Whisnant AW et al. 2020. Integrative functional genomics decodes herpes
  simplex virus 1. Nat Commun 11:2038. doi:10.1038/s41467-020-15992-5
  (annotation deposited as GenBank BK012101.1).
- Colgrove RC et al. 2015. History and genomic sequence analysis of the herpes
  simplex virus 1 KOS and KOS1.1 sub-strains. Virology 487:215-221
  (GenBank KT899744.1).
- Li H. 2018. Minimap2: pairwise alignment for nucleotide sequences.
  Bioinformatics 34:3094-3100.
