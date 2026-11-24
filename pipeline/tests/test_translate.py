"""Codon-exact derivation tests on a synthetic mini-genome.

The fixture is a 400nt chromosome with one + strand transcript (T1) and one
− strand transcript (T2) sharing two CDS blocks split by a 49nt gap. Codon 17
(T1) and codon 15 (T2) span the exon junction, so passing these tests means the
machinery assembles codons across junctions and handles the minus strand
correctly — not merely that a string was sliced."""
from __future__ import annotations

from pathlib import Path

import pytest

from neoforge.artifacts import ArtifactStore, StageRun
from neoforge.gtf import GtfTranscriptIndex
from neoforge.sequence import InMemoryGenome
from neoforge.stages.translate import CODON_TABLE, complement, run_translate

PROTEIN1 = "MSKFEKTGKLVVQKRFLGSPPGCVAVAATIS"  # 31 aa, codon-exact (93 nt CDS)
# Exclude codons whose reverse complement is a stop, so both strand orientations
# of the fixture translate stop-free.
_MINUS_UNSAFE = {"TAA", "TAG", "TGA", "TTA", "CTA", "TCA"}
FIRST_CODON = {
    aa: next(c for c, a in CODON_TABLE.items() if a == aa and c not in _MINUS_UNSAFE)
    for aa in set(PROTEIN1)
}
CDS = "".join(FIRST_CODON[aa] for aa in PROTEIN1)

# exon 1 (101-150) carries 50 nt; exon 2 (200-242) carries 43 nt; the CDS is
# 93 nt, so codons tile exactly on both strands, and the junction codon spans
# the intron on each strand (plus: offsets 48-50, minus: offsets 42-44).
GENOME = "N" * 100 + CDS[:50] + "N" * 49 + CDS[50:] + "N" * 158

_CDS_ON_MINUS = complement(CDS)[::-1]
PROTEIN2 = "".join(CODON_TABLE[_CDS_ON_MINUS[i:i + 3]] for i in range(0, len(_CDS_ON_MINUS) - 2, 3))
assert "*" not in PROTEIN1 and "*" not in PROTEIN2
assert len(CDS) == 93


def alternate(base: str) -> str:
    return next(b for b in "ACGT" if b != base)


def alt_for(codon: str, codon_pos: int) -> str:
    """Alt base giving a non-stop, non-synonymous substitution in one codon."""
    ref_aa = CODON_TABLE[codon]
    for b in "ACGT":
        if b == codon[codon_pos]:
            continue
        aa = CODON_TABLE[codon[:codon_pos] + b + codon[codon_pos + 1:]]
        if aa != "*" and aa != ref_aa:
            return b
    raise AssertionError("fixture codon admits no clean substitution")


def alt_for_both(plus_codon: str, plus_pos: int, minus_codon: str, minus_pos: int) -> str:
    """Alt base giving a clean substitution on both strands simultaneously."""
    return next(
        b for b in "ACGT"
        if b != plus_codon[plus_pos]
        and CODON_TABLE[plus_codon[:plus_pos] + b + plus_codon[plus_pos + 1:]]
        not in ("*", CODON_TABLE[plus_codon])
        and CODON_TABLE[minus_codon[:minus_pos] + complement(b) + minus_codon[minus_pos + 1:]]
        not in ("*", CODON_TABLE[minus_codon])
    )


def find_single_base_change(cds: str, length: int, predicate) -> tuple[int, str] | None:
    """Find (cds_index, alt_base) with predicate(ref_aa, alt_aa), interior only."""
    for i in range(1, length - 3):
        codon_pos, codon_idx = i % 3, i // 3
        codon = cds[codon_idx * 3:codon_idx * 3 + 3]
        for b in "ACGT":
            if b == codon[codon_pos]:
                continue
            if predicate(CODON_TABLE[codon], CODON_TABLE[codon[:codon_pos] + b + codon[codon_pos + 1:]]):
                return i, b
    return None


def make_gtf(tmp_path: Path, tids: tuple[str, ...] = ("T1", "T2")) -> Path:
    path = tmp_path / ("mini-" + "".join(tids) + ".gtf")
    lines = ["#test fixture"]
    for tid in tids:
        strand = "+" if tid == "T1" else "-"
        for start, end in ((101, 150), (200, 242)):
            lines.append(
                f'chrT\tsrc\tCDS\t{start}\t{end}\t.\t{strand}\t0\t'
                f'gene_id "G"; transcript_id "{tid}"; transcript_type "protein_coding";'
            )
    path.write_text("\n".join(lines) + "\n")
    return path


def write_variants(store: ArtifactStore, rows: list[tuple[str, int, str, str]], case: str = "mini") -> None:
    payload = {
        "selection": {"ingested_snv_records": len(rows), "selected": len(rows),
                      "excluded": 0, "exclusion_reasons": {}, "indels_seen": 0},
        "selected": {
            "columns": ["variant_key", "chrom", "pos", "ref", "alt", "tvaf", "nvaf", "tumor_depth"],
            "variant_key": [f"{c}:{p}{r}>{a}" for c, p, r, a in rows],
            "chrom": [c for c, _, _, _ in rows],
            "pos": [p for _, p, _, _ in rows],
            "ref": [r for _, _, r, _ in rows],
            "alt": [a for _, _, _, a in rows],
            "tvaf": [None] * len(rows),
            "nvaf": [None] * len(rows),
            "tumor_depth": [None] * len(rows),
        },
    }
    store.write(StageRun(stage="variants", case=case), "variants", payload)


def run(mini, tmp_path: Path, rows, lengths=(8,), tids=("T1", "T2")) -> dict:
    _index, genome = mini
    store = ArtifactStore(tmp_path / f"work-{rows[0][1]}-{rows[0][3]}-{len(tids)}")
    write_variants(store, rows)
    return run_translate(store, "mini", make_gtf(tmp_path, tids), genome, peptide_lengths=lengths)


def as_rows(derived: dict) -> list[dict]:
    cols = derived["columns"]
    return [dict(zip(cols, vals)) for vals in zip(*[derived[c] for c in cols])]


@pytest.fixture()
def mini(tmp_path: Path) -> tuple[GtfTranscriptIndex, InMemoryGenome]:
    return GtfTranscriptIndex.from_gtf(make_gtf(tmp_path)), InMemoryGenome({"chrT": GENOME})


def test_codon_table_is_complete_standard_code() -> None:
    assert len(CODON_TABLE) == 64
    assert sum(1 for v in CODON_TABLE.values() if v == "*") == 3
    assert CODON_TABLE["ATG"] == "M" and CODON_TABLE["TTT"] == "F"


def test_index_coverage_and_translation_order(mini) -> None:
    index, _ = mini
    t1, t2 = index.transcripts["T1"], index.transcripts["T2"]
    assert set(index.covering_transcripts("chrT", 150)) == {"T1", "T2"}
    assert index.covering_transcripts("chrT", 151) == []
    assert index.covering_transcripts("chrT", 199) == []
    assert t1.offset_of(150) == 49 and t1.offset_of(200) == 50
    assert t2.offset_of(150) == 43 and t2.offset_of(200) == 42
    assert t1.positions_for_offsets([48, 49, 50]) == [149, 150, 200]
    # minus strand: translation order walks the genome right-to-left
    assert t2.positions_for_offsets([42, 43, 44]) == [200, 150, 149]
    assert t1.n_codons == 31 and t2.n_codons == 31


def test_junction_codon_plus_strand(mini, tmp_path: Path) -> None:
    alt = alt_for(CDS[48:51], 1)
    out = run(mini, tmp_path, [("chrT", 150, CDS[49], alt)])
    rows = [r for r in as_rows(out["payload"]["derived"]) if r["transcript_id"] == "T1"]
    assert rows, "plus-strand transcript must yield windows spanning the junction"
    expected_ref_aa = CODON_TABLE[CDS[48:51]]
    expected_alt_aa = CODON_TABLE[CDS[48] + alt + CDS[50]]
    for row in rows:
        assert row["protein_pos"] == 17
        assert row["ref_aa"] == expected_ref_aa
        assert row["alt_aa"] == expected_alt_aa
        assert expected_alt_aa not in (expected_ref_aa, "*")
        s, k = row["window_start"], len(row["window_seq"])
        assert row["ref_window_seq"] == PROTEIN1[s - 1:s + k - 1]
        assert row["window_seq"] == expected_alt_aa.join([PROTEIN1[s - 1:16], PROTEIN1[17:s + k - 1]])
        assert row["mutant_offset"] == 17 - s + 1
    assert [r["window_start"] for r in rows] == list(range(max(1, 17 - 8 + 1), min(31 - 8 + 1, 17) + 1))


def test_junction_codon_minus_strand(mini, tmp_path: Path) -> None:
    minus_codon = _CDS_ON_MINUS[42:45]
    alt = alt_for(minus_codon, 1)
    out = run(mini, tmp_path, [("chrT", 150, CDS[49], alt)])
    rows = [r for r in as_rows(out["payload"]["derived"]) if r["transcript_id"] == "T2"]
    assert rows, "minus-strand transcript must yield windows spanning the junction"
    expected_ref_aa = CODON_TABLE[minus_codon]
    assert expected_ref_aa == PROTEIN2[14]
    for row in rows:
        assert row["protein_pos"] == 15
        assert row["ref_aa"] == expected_ref_aa
        expected_alt_aa = CODON_TABLE[_CDS_ON_MINUS[42] + complement(alt) + _CDS_ON_MINUS[44]]
        assert row["alt_aa"] == expected_alt_aa
        assert expected_alt_aa not in (expected_ref_aa, "*")
        s, k = row["window_start"], len(row["window_seq"])
        assert row["ref_window_seq"] == PROTEIN2[s - 1:s + k - 1]
        assert row["window_seq"] == expected_alt_aa.join([PROTEIN2[s - 1:14], PROTEIN2[15:s + k - 1]])


def test_window_count_formula_across_lengths(mini, tmp_path: Path) -> None:
    plus_codon, minus_codon = CDS[18:21], _CDS_ON_MINUS[72:75]
    alt = alt_for_both(plus_codon, 0, minus_codon, 2)
    out = run(mini, tmp_path, [("chrT", 119, plus_codon[0], alt)], lengths=(8, 9, 10, 11))
    rows = as_rows(out["payload"]["derived"])
    assert {r["transcript_id"] for r in rows} == {"T1", "T2"}
    n = len(rows)
    interior = sum(min(k, 7, 31 - 7 + 1) for k in (8, 9, 10, 11))
    assert n == 2 * interior  # both transcripts share the same geometry


def test_first_residue_window_edges(mini, tmp_path: Path) -> None:
    plus_codon, minus_codon = CDS[0:3], _CDS_ON_MINUS[90:93]
    alt = alt_for_both(plus_codon, 0, minus_codon, 2)
    out = run(mini, tmp_path, [("chrT", 101, CDS[0], alt)])
    rows = as_rows(out["payload"]["derived"])
    t1 = [r for r in rows if r["transcript_id"] == "T1"]
    t2 = [r for r in rows if r["transcript_id"] == "T2"]
    assert t1 and t2
    assert [r["window_start"] for r in t1] == [1] and [r["mutant_offset"] for r in t1] == [1]
    assert [r["protein_pos"] for r in t2] == [31]
    assert [r["window_start"] for r in t2] == [24] and [r["mutant_offset"] for r in t2] == [8]


def test_exclusions_are_counted_not_silent(mini, tmp_path: Path) -> None:
    syn = find_single_base_change(CDS, 93, lambda a, b: a == b)
    stop = find_single_base_change(CDS, 93, lambda a, b: b == "*")
    assert syn is not None and stop is not None
    rows = [
        ("chrT", 150, "N", "A"),
        ("chrT", 101 + syn[0], CDS[syn[0]], syn[1]),
        ("chrT", 101 + stop[0], CDS[stop[0]], stop[1]),
    ]
    out = run(mini, tmp_path, rows, tids=("T1",))
    summary = out["payload"]["summary"]
    assert summary["ref_mismatch"] == 1
    assert summary["synonymous"] == 1
    assert summary["stop_gained"] == 1
    assert summary["windows_derived"] == 0
    assert summary["variants_processed"] == 3


def test_ref_windows_match_independent_translation(mini, tmp_path: Path) -> None:
    plus_codon, minus_codon = CDS[9:12], _CDS_ON_MINUS[81:84]
    alt = alt_for_both(plus_codon, 1, minus_codon, 1)
    out = run(mini, tmp_path, [("chrT", 111, CDS[10], alt)])
    rows = as_rows(out["payload"]["derived"])
    by_transcript = {r["transcript_id"]: r for r in rows}
    # each transcript's ref windows must be slices of a protein derived
    # independently here from the same fixture strings
    for tid, protein in (("T1", PROTEIN1), ("T2", PROTEIN2)):
        row = by_transcript[tid]
        expected = protein[row["window_start"] - 1:row["window_start"] + 7]
        assert row["ref_window_seq"] == expected
