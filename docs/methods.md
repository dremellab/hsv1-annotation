# Draft methods paragraph

Draft text for a publication. Check the numbers against
`output/KT899744delRR/KT899744delRR.annotation_summary.tsv` for the release used.

---

HSV-1 reads were mapped to KT899744delRR, a version of the KOS genome
KT899744.1 (Colgrove et al. 2015) in which the terminal repeats TRL and TRS
and a 1,133 bp terminal sequence copy at the IRL/IRS junction were removed
(135,164 bp). The annotation of strain KOS published by Whisnant et al. (2020;
GenBank BK012101.1) was transferred to this genome. BK012101.1 was aligned to
KT899744delRR with minimap2 v2.30 (Li 2018; preset asm20 with `--cs --eqx`,
`-N 50 -p 0.3`), and features in the collinear alignment block were lifted by
coordinate. Features in the removed terminal repeats were dropped, features
crossing a deleted boundary were clipped, and each lifted CDS was
re-translated and compared with the source protein. This gave 258 ORFs, with
one gene per ORF, plus 18 noncoding transcripts. For read counting, 40 ORFs
whose coding bases all lie inside one other same-strand ORF were excluded, so
that reads in these regions are assigned to the containing kept ORF; the remaining
218 ORFs (143 non-overlapping and 75 partially overlapping) were used with
featureCounts (strand-aware counting, `-t CDS -g gene_id`, without counting multi-overlap
reads), which discards reads overlapping more than one ORF. Transcription
start and end sites were taken from the 5' and 3' ends of the BK012101.1
mRNAs and are included in the annotation used for visualization only.
The annotation files and build scripts are available at
https://github.com/dremellab/hsv1-annotation (release v0.1.0).
