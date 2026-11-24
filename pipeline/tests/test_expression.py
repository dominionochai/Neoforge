"""Expression join tests: crosswalk ambiguity, log2(TPM+1) inversion, honesty."""
from __future__ import annotations

import math

from neoforge.expression import crosswalk, gene_tpm


def test_crosswalk_maps_and_drops_ambiguity(tmp_path):
    p = tmp_path / "map.tsv"
    p.write_text("ENSGA\t111\nENSGA\t222\nENSG\t111\n")
    ensg_by_entrez, dropped = crosswalk(p)
    assert ensg_by_entrez == {"222": "ENSGA"}
    assert dropped == 1


def test_gene_tpm_inverts_the_documented_transform(tmp_path):
    # source stores log2(TPM+1); tpm = 2^v - 1 must round-trip exactly
    p = tmp_path / "row.csv"
    p.write_text(
        "depmap_id,TSPAN6 (7105),CFTR (1080),NOENTREZ,DUP (9999),GAP基因 (7777),EMPTY (555)\n"
        "ACH-TEST," + f"{math.log2(22.0)},{math.log2(33.0)},{math.log2(4.0)},{math.log2(11.0)},{math.log2(1000.0)},\n")
    mapping, dropped = crosswalk(_write_map(tmp_path))
    depmap_id, tpm, counts = gene_tpm(p, mapping)
    assert depmap_id == "ACH-TEST"
    assert tpm["GH"] == 21.0
    assert tpm["GH2"] == 32.0
    assert tpm["GH4"] == 10.0
    assert tpm["GH5"] == 999.0
    assert counts["columns"] == 6
    assert counts["value_gaps"] == 0 and counts["unmapped_entrez"] == 2
    assert dropped == 0


def _write_map(tmp_path):
    p = tmp_path / "map.tsv"
    p.write_text("GH\t7105\nGH2\t1080\nGH3\t1234\nGH4\t9999\nGH5\t7777\n")
    return p


def test_ingest_expression_block(tmp_path):
    from neoforge.expression import crosswalk, gene_tpm

    mapping_path = _write_map(tmp_path)
    mapping, _ = crosswalk(mapping_path)
    row_path = tmp_path / "row.csv"
    row_path.write_text(
        "depmap_id,TSPAN6 (7105),CFTR (1080),DUP (9999),GAP (7777),EMPTY (555)\n"
        "ACH-TEST," + f"{math.log2(22.0)},{math.log2(33.0)},{math.log2(4.0)},{math.log2(1000.0)},\n")
    depmap_id, tpm, counts = gene_tpm(row_path, mapping)
    assert depmap_id == "ACH-TEST"
    assert set(tpm) == {"GH", "GH2", "GH4", "GH5"}
    assert counts["unmapped_entrez"] == 1 and counts["value_gaps"] == 0
