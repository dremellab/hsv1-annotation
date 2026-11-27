from check_gtf import check

GENOME = "ACGT" * 25  # 100 bp


def row(feat, s, e, attrs, strand="+", frame="."):
    return ["KT899744", "Genbank", feat, str(s), str(e), ".", strand, frame, attrs]


def gene_rows(gene_attr="", cds=True, tx_attr=""):
    g = 'gene_id "G1";' + gene_attr
    t = 'gene_id "G1"; transcript_id "T1";' + tx_attr
    rows = [
        row("gene", 1, 30, g),
        row("transcript", 1, 30, t),
        row("exon", 1, 30, t),
    ]
    if cds:
        rows.append(row("CDS", 1, 30, t, frame="0"))
    return rows


def test_coding_gene_with_cds_passes():
    assert check(gene_rows(), GENOME) == []


def test_space_in_value_fails():
    rows = gene_rows(tx_attr=' transcript_name "a b";')
    assert any("space" in e for e in check(rows, GENOME))


def test_missing_transcript_id_fails():
    rows = gene_rows()
    rows[2][8] = 'gene_id "G1";'
    assert any("transcript_id" in e for e in check(rows, GENOME))


def test_noncoding_gene_without_cds_passes():
    rows = gene_rows(gene_attr=' gene_biotype "noncoding";', cds=False)
    assert check(rows, GENOME) == []


def test_coding_gene_without_cds_fails():
    rows = gene_rows(cds=False)
    assert any("without CDS" in e for e in check(rows, GENOME))


def test_out_of_bounds_fails():
    rows = gene_rows()
    rows[0][4] = "101"
    assert any("out of bounds" in e for e in check(rows, GENOME))


def test_tss_tes_rows_standalone_pass():
    rows = [
        row("TSS", 5, 5, 'gene_id "G9"; transcript_id "T9";'),
        row("TES", 50, 50, 'gene_id "G9"; transcript_id "T9";'),
    ]
    assert check(rows, GENOME) == []
    assert any("out of bounds" in e for e in check([row("TSS", 101, 101, 'gene_id "G9"; transcript_id "T9";')], GENOME))
    assert any("transcript_id" in e for e in check([row("TES", 5, 5, 'gene_id "G9";')], GENOME))
    assert any("space" in e for e in check([row("TSS", 5, 5, 'gene_id "G 9"; transcript_id "T9";')], GENOME))


def test_mirna_gene_without_cds_passes():
    rows = gene_rows(gene_attr=' gene_biotype "miRNA";', cds=False)
    assert check(rows, GENOME) == []


def test_untyped_gene_without_cds_fails():
    rows = gene_rows(gene_attr=' gene_biotype "other";', cds=False)
    assert any("without CDS" in e for e in check(rows, GENOME))
