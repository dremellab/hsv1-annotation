# Changelog

All notable changes to this project are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- UCSC assembly hub (`hub/`) with visualization, quantification, excluded-ORF
  and TSS/TES tracks, and `hub/build_hub.sh` to rebuild it.
- HSV-1 miRNAs from miRBase 23.1: `scripts/place_mirnas.py` places 17 of 18
  precursors and 25 of 27 mature miRNAs on KT899744delRR (exact match,
  projection from BK012101.1, or local alignment; H11 on both strands, H5
  unplaced because it lies in the deleted blocks). New outputs
  `KT899744delRR.miRNA.{gtf,bed}` and `miRNA_report.tsv`, miRNAs in the master
  and visualization GTFs, a hub track, and `mirna`, `mirna_precursors` and
  `mirna_unplaced` summary rows. The quantification GTF is unchanged.
- miRBase sources in `sources/miRBase_23.1/`, fetched and filtered by
  `scripts/fetch_sources.sh`.

## [v0.1.0] - 2026-11-27

### Added

- Source records (BK012101.1, KT899744.1), KT899744delRR provenance and
  verification script.
- Liftover of the BK012101.1 annotation onto KT899744delRR, keeping the 18
  noncoding transcripts.
- Master GTF with ORF genes, transcripts, TSS and TES rows.
- Class definitions: `non_overlapping` and `partially_overlapping` are decided
  against all same-strand ORFs (kept or excluded), so the host of an embedded
  ORF is `partially_overlapping`; effective unique length is computed against
  kept ORFs only.
- TSS and TES rows are skipped for transcripts clipped at a genome boundary
  when the end lies at the clip point (position 1 or the end of the genome).
- ORF classification and the quantification and visualization GTFs, with
  excluded-CDS, overlap and summary reports.
- GTF validator, pipeline-helper check and unit tests.
- `build.sh` to rebuild all outputs and `publish.sh` to copy the GTFs into the
  shared reference folder.
- README and draft publication methods (`docs/methods.md`).
