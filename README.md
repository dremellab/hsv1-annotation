# hsv1-annotation

GTF annotations of the herpes simplex virus 1 (HSV-1) genome KT899744delRR,
built for read counting (quantification) and for genome browsers
(visualization). Version: see `VERSION`.

The annotation is the BK012101.1 annotation (Whisnant et al. 2020) lifted
onto KT899744delRR. One gene is defined per ORF, ORFs that lie entirely
inside another ORF are removed from the quantification GTF, and TSS and TES
features are added to the visualization GTF. Host-virus RNA-seq and Ribo-seq
pipelines (HAROLD, RIBOLD) can use the quantification GTF without code changes.

HSV-1 miRNAs (miRBase 23.1) are placed on KT899744delRR and added to the
master and visualization GTFs; they have their own GTF, BED and report and are
never in the quantification GTF (see "miRNAs").

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
   KT899744.1 (FASTA) and the miRBase HSV-1 files are checked against
   `sources/checksums.md5` (`build.sh` runs `scripts/fetch_sources.sh` to
   download them only when this check fails), and
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
   visualization GTF and the reports are written. miRNAs (placed by
   `scripts/place_mirnas.py`, see "miRNAs") are merged into the master before
   this step and excluded from the quantification GTF.
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
| `KT899744delRR.master.gtf` | 258 ORF genes plus 18 noncoding genes and 18 miRNA precursor genes, all transcripts, TSS and TES rows |
| `KT899744delRR.quantification.gtf` | 218 ORFs, CDS-carrying transcripts only; no TSS, TES or noncoding transcripts |
| `KT899744delRR.visualization.gtf` | Same content as the master (including the miRNAs), sorted, attributes tidied for genome browsers |
| `KT899744delRR.miRNA.gtf` | miRNA precursor and mature records only (source `miRBase_23.1`); use this for small RNA-seq counting |
| `KT899744delRR.miRNA.bed` | The same records as BED6, name = `transcript_id` |
| `KT899744delRR.miRNA_report.tsv` | One row per placed miRNA record and one per unplaced miRBase entry: name, miRBase ID, type, coordinates, strand, placement, identity, alignment length, flags |
| `KT899744delRR.excluded_CDSs.tsv` | One row per excluded ORF: CDS ID, gene name, chromosome, start, end, strand, CDS length, containing CDS ID, overlapping bases, unique bases, reason, `cds_check` |
| `KT899744delRR.CDS_overlaps.tsv` | One row per overlapping ORF pair (either strand): each ORF's class, full length and effective unique length |
| `KT899744delRR.annotation_summary.tsv` | Totals: CDSs, non-overlapping, partially overlapping, embedded, excluded, retained, `retained_zero_unique`, miRNA mature records, miRNA precursors, miRNA entries unplaced, TSS, TES |
| `liftover_report.tsv` | Per BK012101.1 transcript: lift status, coordinates, protein ID, `cds_check`, amino-acid differences |
| `bk_to_delrr.paf` | The minimap2 alignment used for lifting |

All GTFs use sequence name `KT899744`. Lifted rows use source `Genbank` and gene
IDs prefixed `KT899744delRR_`; miRNA rows use source `miRBase_23.1` and ids
`hsv1-mir-*` (precursors) and `hsv1-miR-*` (mature arms), with no prefix.

## miRNAs

**Source.** The HSV-1 entries (`hsv1-`) of miRBase release 23.1 (Kozomara et al.
2019): 18 hairpins and 27 mature miRNAs, stored as
`sources/miRBase_23.1/hsv1_hairpin.fa` and `hsv1_mature.fa` (sequences as
published, U kept). `scripts/fetch_sources.sh` downloads the full
`hairpin.fa` and `mature.fa` from miRBase and filters them, so the committed
files are reproducible; their MD5s are in `sources/checksums.md5`.

**Result.** 17 of 18 precursors and 25 of 27 mature miRNAs are placed on
KT899744delRR. With H11 annotated on both strands this gives 18 precursor
records and 26 mature records (the H11-3p arm has one record per strand).
Of the 17 placed precursors, 12 are exact matches on KT899744delRR, 3 are
projected from BK012101.1 (H3, H18, and H15 after refinement) and 2 come from local alignment
(H8, H17). The per-miRNA table is `output/KT899744delRR/KT899744delRR.miRNA_report.tsv`
(columns: `name mirbase_id type chrom start end strand placement identity_pct
aln_length flags`; coordinates are 1-based on `KT899744`).

| miRNA | Placement | Coordinates (strand) | Identity | Flags |
|---|---|---|---|---|
| H1 | exact_target | 109078-109163 (+) | 100.0 | none |
| H2 | exact_target | 112632-112706 (+) | 100.0 | none |
| H3 | projected_from_BK012101.1 | 116519-116599 (+) | NA | none |
| H4 | exact_target | 116664-116744 (+) | 100.0 | none |
| H5 | unplaced | not in KT899744delRR | NA | in_deleted_region |
| H6 | exact_target | 109078-109163 (-) | 100.0 | none |
| H7 | exact_target | 114202-114284 (+) | 100.0 | none |
| H8 | local_alignment | 114597-114670 (+) | 97.3 | approximate |
| H11 | exact_target | 53208-53337 (+ and -) | 100.0 | strand_ambiguous |
| H12 | exact_target | 121584-121671 (-) | 100.0 | none |
| H13 | exact_target | 121595-121658 (+) | 100.0 | none |
| H14 | exact_target | 112625-112713 (-) | 100.0 | none |
| H15 | projected_refined | 109955-110052 (-) | 100.0 | length_changed |
| H16 | exact_target | 59819-59918 (+) | 100.0 | none |
| H17 | local_alignment | 117183-117277 (-) | 97.9 | approximate |
| H18 | projected_from_BK012101.1 | 121460-121554 (+) | NA | none |
| H26 | exact_target | 83681-83743 (+) | 100.0 | none |
| H27 | exact_target | 114201-114283 (-) | 100.0 | none |

### How the miRNA annotation was made

1. **Download and filter.** `scripts/fetch_sources.sh` (run by `build.sh` only
   when the checksums fail) downloads miRBase 23.1
   `hairpin.fa` and `mature.fa`, keeps the `>hsv1-` entries and writes
   `sources/miRBase_23.1/hsv1_{hairpin,mature}.fa`; the MD5s are checked.
2. **Place precursors.** `scripts/place_mirnas.py` converts U to T and
   searches both strands. (Report names of H11 records carry `_plus` or
   `_minus`, the same suffix as the GTF ids.) For each hairpin it tries, in order:
   1. `exact_target`: exact match on KT899744delRR. One hit is placed. Two
      hits at the same start on opposite strands (a hairpin equal to its own
      reverse complement, H11) are both kept and flagged `strand_ambiguous`.
      Any other multiple hit is kept and flagged `multi_hit`.
   2. `projected_from_BK012101.1`: otherwise, exact match on BK012101.1, with
      start and end projected through the collinear minimap2 block (the
      `build_map` function of `scripts/lift_bk012101.py`). A hit is kept if
      both ends project. If the projection changes the length, the hit is
      kept and flagged `length_changed`. Such a projection can be stretched by
      an indel at a homopolymer edge (H15: an A9 run in BK012101.1 is A28 in
      KT899744delRR), so it is refined: the hairpin is aligned
      (Smith-Waterman, as in 3) within 150 nt of the projected interval on the
      same strand, the placement becomes `projected_refined` and identity and
      `aln_length` come from that alignment (the refinement is used only if the
      alignment covers at least 80% of the hairpin at 80% identity or better).
      Identity is `NA` for other projected hits.
   3. `local_alignment`: otherwise (H8, H17), the mature arm(s) are placed
      by methods 1 and 2, then the hairpin is aligned (Smith-Waterman,
      match +2, mismatch -3, gap -5) to KT899744delRR within 150 nt of the
      placed arm on the arm's strand. Identity (%) and alignment length are
      reported and the flag `approximate` is set.
   4. `unplaced`: no GTF rows, only a report row. The flag is
      `in_deleted_region` if the sequence matches KT899744.1 exactly inside
      a block removed from KT899744delRR (blocks from `scripts/verify_delrr.py`),
      otherwise `no_match`. A mature arm that matches the genome but not inside
      its precursor gets `outside_precursor`.
3. **Place mature arms.** Each mature miRNA is placed by methods 1 and 2 on
   its own and kept only where it lies inside a placed precursor of the same
   miRNA (`hsv1-miR-H1-5p` belongs to `hsv1-mir-H1`) on the same strand. A
   mature record reports its own placement method and identity, and carries
   `strand_ambiguous` from its precursor.
4. **H11.** See below. H11 is annotated on both strands.
5. **Outputs.** `place_mirnas.py` writes `KT899744delRR.miRNA.gtf`,
   `.miRNA.bed` and `.miRNA_report.tsv`. Each precursor has a `gene`,
   `transcript` (`transcript_biotype "pre_miRNA"`) and `exon` row, and each
   mature arm a `transcript` and `exon` row (`transcript_biotype "miRNA"`);
   all carry `mirbase_id`, `placement` and `flags`. `build_master.py
   --mirna-gtf` merges them into the master; `derive_gtfs.py` keeps genes with
   `gene_biotype "miRNA"` out of the quantification GTF and reports the
   `mirna`, `mirna_precursors` and `mirna_unplaced` summary rows;
   `check_gtf.py` accepts miRNA genes without a CDS. The UCSC hub gets an
   `mirna` track (`hub/build_hub.sh`).
6. **Rebuild.** `./build.sh` runs steps 1 to 5 (placement is step 4b, after
   the minimap2 alignment because it needs the PAF); `./hub/build_hub.sh`
   rebuilds the hub. Both are reproducible: two builds give identical files.

### H11 strand decision

The H11 precursor is a perfect 130 bp inverted repeat in OriL, so the hairpin
equals its own reverse complement and matches KT899744delRR at 53208-53337 on
both strands. In BK012101.1, transcripts from both strands span it (the
UL27.5-UL30 family on the + strand, the UL29.5-UL30.6 family on the - strand),
so the data do not say which strand makes the miRNA. It is annotated on both
strands, flagged `strand_ambiguous`, and the ids get `_plus` and `_minus`
suffixes (`hsv1-mir-H11_plus`, `hsv1-miR-H11-3p_minus`). The 3p arm is at
53280-53299 (+) and 53246-53265 (-). See Zou et al. 2020 and Duan et al.
2023.

### H5 is not in KT899744delRR

H5 (`hsv1-mir-H5` and both arms) matches KT899744.1 at 126577-126656 (+) and
151154-151233 (-), and BK012101.1 at 126720-126799 (+) and 151398-151477 (-).
Both KT899744.1 positions lie in blocks removed from KT899744delRR
(125845-126977 and 145361-151974), so H5 has no sequence on this genome. It is
reported as `unplaced` with flag `in_deleted_region`, and reads from the H5
locus cannot map to KT899744delRR. Other entries are not affected.

### H1 and H6

H1 and H6 are an antisense pair: both precursors start at 109078 and span
109078-109163, H1 on the + strand and H6 on the - strand. They have different
sequences, so they are two ordinary records, not `strand_ambiguous`.

### Use

The miRNAs are in the master and visualization GTFs and in the miRNA GTF and
BED, never in the quantification GTF (it is byte-identical to the previous
version). For small RNA-seq counting use `KT899744delRR.miRNA.gtf`; count the
mature records (`transcript_biotype "miRNA"`), not the precursors.

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
the TSS/TES features and the miRNA precursor and mature records. Do not use it
for counting. Checks in IGV and UCSC
have not been recorded yet; any feature types a browser rejects will be listed
here.

## UCSC Genome Browser hub

`hub/` is a UCSC assembly hub for KT899744delRR (KT899744delRR is not a UCSC
genome, so the hub carries its own sequence).

**Open in UCSC:**
[KT899744delRR hub at UL26/UL26.5](https://genome.ucsc.edu/cgi-bin/hgTracks?genome=KT899744delRR&hubUrl=https://raw.githubusercontent.com/dremellab/hsv1-annotation/main/hub/hub.txt&position=KT899744:41000-44000).
This link is safe to share: it loads the hub from this repo into the
viewer's own browser session and always shows the hub on `main`.

Alternatively, go to
<https://genome.ucsc.edu/cgi-bin/hgHubConnect>, choose "Connected Hubs" and
paste:

```
https://raw.githubusercontent.com/dremellab/hsv1-annotation/main/hub/hub.txt
```

Tracks:

| Track | Source | Items |
|---|---|---|
| Visualization GTF | `KT899744delRR.visualization.gtf` | 347 transcripts |
| Quantification GTF | `KT899744delRR.quantification.gtf` | 245 transcripts |
| Excluded ORFs | `KT899744delRR.excluded_CDSs.tsv`, labelled `ORF_in_HOST` | 40 |
| TSS / TES | TSS (green) and TES (red) rows of the visualization GTF | 290 + 157 |
| miRNAs (miRBase 23.1) | `KT899744delRR.miRNA.bed`: 18 precursor and 26 mature records, H11 on both strands | 44 |

Items are labelled by ORF name. The default view (KT899744:41,000-44,000)
shows UL26 with the embedded UL26.5, which is present in the visualization
track and absent from the quantification track. Rebuild after any change to
`output/`:

```bash
./hub/build_hub.sh   # UCSC tools from the biocondatools env; override with UCSC_BIN
```

## Tests

```bash
/project/dremel_lab/miniforge3/envs/pipelines/bin/python -m pytest tests
```

The tests use small synthetic GTFs: each class, strand-aware overlap,
containment across several transcripts of one gene, the identical-footprint
tie-break, chains, and effective unique length after exclusion. The miRNA
tests (`tests/test_place_mirnas.py`) use toy genomes for each placement path:
a palindromic hairpin on both strands, an antisense pair, projection, refinement
of a projection that changes length, local alignment, mature arms kept only
inside their precursor, and the unplaced reasons. They also check, on the real
sources (needs minimap2), the placement counts and that the quantification GTF
is unchanged.

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
- Kozomara A, Birgaoanu M, Griffiths-Jones S. 2019. miRBase: from microRNA
  sequences to function. Nucleic Acids Res 47:D155-D162.
- Zou W, Zhou X, Wang L, Zhou GG, Chen X. 2020. Degradation of herpes simplex
  virus-1 viral miRNA H11 by vaccinia virus protein VP55 attenuates viral
  replication. Front Microbiol 11:717. doi:10.3389/fmicb.2020.00717
- Duan Y, Sun L, Li Q. 2023. Herpes simplex virus 1 microRNAs: an update.
  Intervirology 66:97-110. doi:10.1159/000531348
