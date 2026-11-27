# Changelog

All notable changes to this project are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

## [v0.1.0]

### Added

- Source records (BK012101.1, KT899744.1), KT899744delRR provenance and
  verification script.
- Liftover of the BK012101.1 annotation onto KT899744delRR, keeping the 18
  noncoding transcripts.
- Master GTF with ORF genes, transcripts, TSS and TES rows.
- ORF classification and the quantification and visualization GTFs, with
  excluded-CDS, overlap and summary reports.
- GTF validator, pipeline-helper check and unit tests.
- `build.sh` to rebuild all outputs and `publish.sh` to copy the GTFs into the
  shared reference folder.
- README and draft publication methods (`docs/methods.md`).
