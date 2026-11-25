"""Stage 7 construct tests.

The construct stage is a design layer, so these tests pin the discipline, not
a model: mechanical refusals without the full L1-L3 chain or HLA, ledgered
selection rules (tier pool, mutant-vs-wild-type discrimination, overlaps,
capacity), exact back-translation with the pinned standard code, the pinned
codon table's determinism, and the junction audit's flag semantics and
monotone improvement.

Chain fixtures reuse test_present's mini-genome translate products where a
live chain is needed; the selection fixtures hand-write the chain artifacts
in their real column-major shapes and run the genuine run_prioritize on top,
so tier assignment itself stays under the production code.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from neoforge.artifacts import ArtifactStore, StageRun
from neoforge.errors import ArtifactChainError
from neoforge.stages.construct import (
    BINDER_NM,
    JUNCTION_LENGTHS,
    _translate,
    homopolymer_runs,
    improve_order,
    junction_peptides,
    load_codon_usage,
    reverse_translate,
    run_construct,
)
from neoforge.stages.prioritize import run_prioritize

ALLELES = ["HLA-A*02:01", "HLA-B*08:01"]
GENES_3 = ("ENSG1", "ENSG2", "ENSG3")


def write_ingest(store: ArtifactStore, *, alleles: list[str] | None = ALLELES,
                 expression: dict | None = None) -> None:
    payload: dict = {"hla": {"resolved_alleles": alleles or []},
                     "checks": []}
    if expression is not None:
        payload["expression"] = {"tpm_by_gene": expression,
                                 "source_description": "fixture measured TPM"}
    store.write(StageRun(stage="ingest", case="mini"), "ingest", payload)


def write_variants_fixture(store: ArtifactStore, keys: list[str]) -> None:
    store.write(
        StageRun(stage="variants", case="mini"), "variants",
        {"selection": {"ingested_snv_records": len(keys), "selected": len(keys),
                       "excluded": 0, "exclusion_reasons": {}, "indels_seen": 0},
         "selected": {"columns": ["variant_key", "tvaf", "tumor_depth"],
                      "variant_key": keys,
                      "tvaf": [0.5] * len(keys),
                      "tumor_depth": [100] * len(keys)}})


def write_translate(store: ArtifactStore, rows: list[tuple[str, str, str, int]]) -> None:
    store.write(
        StageRun(stage="translate", case="mini"), "translate",
        {"summary": {"windows_derived": len(rows), "transcripts_hit":
                     len({r[1] for r in rows}),
                     "variants_with_peptides": len({r[0] for r in rows}),
                     "variants_processed": len({r[0] for r in rows})},
         "parameters": {"peptide_lengths": [8, 9, 10, 11]},
         "derived": {"columns": ["variant_key", "transcript_id", "gene_id", "mutant_offset"],
                     "variant_key": [r[0] for r in rows],
                     "transcript_id": [r[1] for r in rows],
                     "gene_id": [r[2] for r in rows],
                     "mutant_offset": [r[3] for r in rows]}})


def pred_row(key: str, gene: str, tid: str, allele: str, peptide: str,
             mutant: float, ref: float, offset: int = 0) -> dict:
    return {"variant_key": key, "transcript_id": tid, "gene_id": gene,
            "peptide_length": len(peptide), "window_start": offset,
            "allele": allele, "peptide": peptide,
            "mhcflurry_affinity_nm": mutant,
            "mhcflurry_presentation_score": 0.5,
            "mhcflurry_presentation_percentile": 1.0,
            "reference_peptide": peptide[:-1] + ("E" if peptide[-1] != "E" else "D"),
            "reference_mhcflurry_affinity_nm": ref}


def write_present(store: ArtifactStore, rows: list[dict]) -> None:
    columns: dict[str, list] = {}
    for r in rows:
        for k, v in r.items():
            columns.setdefault(k, []).append(v)
    store.write(
        StageRun(stage="present", case="mini"), "present",
        {"layer": "L3 (predicted)", "evidence_kind": "predicted",
         "predictions": {"columns": sorted(columns), **columns},
         "predictions_made": len(rows), "candidates_scored": len({r["variant_key"] for r in rows})})


def full_store(tmp_path: Path, name: str = "work-construct") -> ArtifactStore:
    """Three expressed mutations on three genes; V1/V3 have discriminating
    windows, V2's reference also binds (no wt discrimination)."""
    store = ArtifactStore(tmp_path / name)
    write_ingest(store, expression={"ENSG1": 5.0, "ENSG2": 5.0, "ENSG3": 5.0})
    keys = ["chrT:100A>G", "chrT:200C>T", "chrT:300G>A"]
    write_variants_fixture(store, keys)
    write_translate(store, [(k, f"ENST{i + 1}", g, 40 * i)
                            for i, (k, g) in enumerate(zip(keys, GENES_3))])
    write_present(store, [
        pred_row(keys[0], GENES_1_3[0], "ENST1", ALLELES[0], "SYFKPEPTI", 120.0, 800.0),
        pred_row(keys[1], GENES_1_3[1], "ENST2", ALLELES[0], "TEPTQVAAA", 80.0, 90.0),
        pred_row(keys[2], GENES_1_3[2], "ENST3", ALLELES[1], "IMMTKEYLL", 200.0, 900.0),
    ])
    run_prioritize(store, "mini")
    return store


GENES_1_3 = GENES_3


class JunctionFakePredictor:
    """Deterministic junction scorer; every junction peptide gets a fixed
    affinity so flag semantics are exact."""

    supported_alleles = ALLELES + ["HLA-C*07:01"]

    def __init__(self, affinity_of=None):
        self._affinity_of = affinity_of or (lambda p: 750.0)

    def predict(self, peptides, alleles, verbose=0):
        if len(alleles) != 1:
            raise ValueError("wrapper contract: one allele per predict() call")
        return pd.DataFrame({
            "peptide": list(peptides),
            "affinity": [self._affinity_of(p) for p in peptides],
            "presentation_score": [0.4] * len(peptides),
            "presentation_percentile": [0.4] * len(peptides),
        })


def fixture_factory(affinity_of=None):
    predictor = JunctionFakePredictor(affinity_of)
    return predictor, "fake-junction-1.0", "FixturePredictor junction test models"


def test_refusal_when_chain_stops_before_present(tmp_path: Path) -> None:
    store = ArtifactStore(tmp_path / "work-refusal")
    write_ingest(store)
    write_variants_fixture(store, ["chrT:100A>G"])
    write_translate(store, [("chrT:100A>G", "ENST1", "ENSG1", 7)])
    with pytest.raises(ArtifactChainError) as err:
        run_construct(store, "mini", junction_predictor_factory=fixture_factory)
    assert "present" in str(err.value) and "refuses" in str(err.value)


def test_refusal_when_hla_unavailable(tmp_path: Path) -> None:
    store = ArtifactStore(tmp_path / "work-nohla")
    write_ingest(store, alleles=[])
    write_variants_fixture(store, ["chrT:100A>G"])
    write_translate(store, [("chrT:100A>G", "ENST1", "ENSG1", 7)])
    write_present(store, [pred_row("chrT:100A>G", "ENSG1", "ENST1", ALLELES[0],
                                   "SYFKPEPTI", 120.0, 800.0)])
    run_prioritize(store, "mini")
    with pytest.raises(ArtifactChainError) as err:
        run_construct(store, "mini", junction_predictor_factory=fixture_factory)
    assert "HLA" in str(err.value)


def test_refusal_when_no_epitope_passes_rules(tmp_path: Path) -> None:
    store = ArtifactStore(tmp_path / "work-noepitope")
    write_ingest(store, expression={"ENSG1": 5.0})
    write_variants_fixture(store, ["chrT:100A>G"])
    write_translate(store, [("chrT:100A>G", "ENST1", "ENSG1", 7)])
    write_present(store, [pred_row("chrT:100A>G", "ENSG1", "ENST1", ALLELES[0],
                                   "TEPTQVAAA", 80.0, 90.0)])
    run_prioritize(store, "mini")
    with pytest.raises(ArtifactChainError) as err:
        run_construct(store, "mini", junction_predictor_factory=fixture_factory)
    assert "no epitope" in str(err.value)


def test_selection_ledger_and_exclusions(tmp_path: Path) -> None:
    store = full_store(tmp_path)
    out = run_construct(store, "mini", junction_predictor_factory=fixture_factory)
    payload = out["payload"]
    included = payload["selection_ledger"]["included"]
    excluded = payload["selection_ledger"]["excluded"]

    assert [e["variant_key"] for e in included] == ["chrT:100A>G", "chrT:300G>A"]
    reasons = {e["variant_key"]: e["reason"] for e in excluded}
    assert reasons["chrT:200C>T"] == "no_mutant_vs_wild_type_discrimination"
    for e in included:
        assert e["mutant_affinity_nm"] <= BINDER_NM <= e["reference_affinity_nm"]
        assert e["selection_reason"].startswith("tier A")
    assert payload["case_summary"]["tier_a_candidates"] == 3
    assert {("tier A" in e["selection_reason"]) for e in included} == {True}
    excluded_keys = {e["variant_key"] for e in excluded}
    assert excluded_keys == {"chrT:200C>T"}
    ledgered = ({e["variant_key"] for e in included} | excluded_keys)
    assert ledgered == {"chrT:100A>G", "chrT:200C>T", "chrT:300G>A"}, \
        "every mutation must appear in exactly one ledger side"


def test_capacity_cap(tmp_path: Path) -> None:
    n = 8
    store = ArtifactStore(tmp_path / "work-cap")
    write_ingest(store, expression={f"ENSG{i}": 5.0 for i in range(n)})
    keys = [f"chrT:{100 * (i + 1)}A>G" for i in range(n)]
    write_variants_fixture(store, keys)
    write_translate(store, [(k, f"ENST{i + 1}", f"ENSG{i}", 40 * i) for i, k in enumerate(keys)])
    write_present(store, [pred_row(k, f"ENSG{i}", f"ENST{i + 1}",
                                   ALLELES[i % 2], "SYFKPEPTI", 120.0, 800.0)
                          for i, k in enumerate(keys)])
    run_prioritize(store, "mini")
    out = run_construct(store, "mini", max_epitopes=5,
                        junction_predictor_factory=fixture_factory)
    payload = out["payload"]
    assert len(payload["selection_ledger"]["included"]) == 5
    assert sum(1 for e in payload["selection_ledger"]["excluded"]
               if e["reason"] == "capacity_reached") == n - 5


def test_construct_payload_is_hashed_deterministic_design(tmp_path: Path) -> None:
    store = full_store(tmp_path, "work-det")
    first = run_construct(store, "mini", junction_predictor_factory=fixture_factory)
    second = run_construct(store, "mini", junction_predictor_factory=fixture_factory)
    assert first["artifact_id"] == second["artifact_id"]
    assert first["payload"]["stamp"].startswith("research design")
    assert set(first["payload"]["sequence_checks"]) == {
        "translation_matches_design", "internal_stop_codons",
        "cds_starts_with_atg", "stop_codon"}
    checks = first["payload"]["sequence_checks"]
    assert checks["translation_matches_design"] is True
    assert checks["internal_stop_codons"] == 0
    assert checks["cds_starts_with_atg"] is True


def test_junction_audit_flags_and_monotone_improvement(tmp_path: Path) -> None:
    store = full_store(tmp_path, "work-junction")
    calls: list[int] = []

    def factory():
        calls.append(1)
        predictor = JunctionFakePredictor(lambda p: 400.0 if "G" in p else 900.0)
        return predictor, "fake-junction-1.0", "FixturePredictor junction test models"

    out = run_construct(store, "mini", junction_predictor_factory=factory)
    audit = out["payload"]["junction_audit"]
    rows = audit["final_rows"]
    assert rows, "junction audit produced no rows"
    for r in rows:
        assert r["flagged"] == (r["mhcflurry_affinity_nm"] <= BINDER_NM)
        assert r["length"] in JUNCTION_LENGTHS
    assert audit["final_flagged_junction_peptides"] == \
        sum(1 for r in rows if r["flagged"])
    flags = [h["flagged_junction_peptides"] for h in audit["history"]
             if h["accepted"]]
    assert flags == sorted(flags, reverse=True), \
        f"accepted orderings must not worsen flags: {flags}"


def _pairwise_cost(bad_pairs: set[tuple[int, int]]):
    def audit_fn(seq: list[int]) -> tuple[int, list[dict]]:
        flags = sum(1 for a, b in zip(seq, seq[1:]) if (a, b) in bad_pairs)
        return flags, []
    return audit_fn


def test_minimiser_reinsertion_reaches_where_swaps_stall() -> None:
    """A junction flag is a property of an epitope epair, so swap-only search
    is genuine hill-climbing on a pair objective: this audit admits optima
    only [0,2,3,1] and [1,3,2,0], and swap-only search from [0,1,2,3] stalls
    at cost 1 while one reinsertion reaches 0."""
    good = {(0, 2), (2, 0), (2, 3), (3, 2), (3, 1), (1, 3)}
    bad_pairs = {(a, b) for a in range(4) for b in range(4)
                 if a != b and (a, b) not in good}
    audit_fn = _pairwise_cost(bad_pairs)

    swap_only, flags_swap, history_swap = improve_order(
        [0, 1, 2, 3], audit_fn, move_family="swap", max_moves=12)
    assert flags_swap == 1 and swap_only == [1, 0, 2, 3]
    assert history_swap[-1]["action"] == "no improving swap or reinsertion"

    both, flags_both, history_both = improve_order([0, 1, 2, 3], audit_fn)
    assert flags_both == 0 and both in ([0, 2, 3, 1], [1, 3, 2, 0])

    accepted = [h["flagged_junction_peptides"] for h in history_both if h["accepted"]]
    assert accepted == sorted(accepted, reverse=True), accepted
    moves = [h["action"] for h in history_both if h["accepted"]][1:]
    assert len(moves) <= 2 * 4, moves
    again_flags = improve_order([0, 1, 2, 3], audit_fn)[1]
    assert again_flags == flags_both


def test_junction_peptides_straddle_and_include_all_lengths() -> None:
    spans = junction_peptides(["AAAAAAAAAAA", "CCCCCCCCCC", "G", "G"])
    j0 = [(j, s) for j, s in spans if j == 0]
    lengths = sorted({len(s) for _, s in j0})
    assert lengths == sorted(JUNCTION_LENGTHS)
    for _, s in j0:
        assert "A" in s and "C" in s, s


def test_reverse_translation_roundtrip_and_codon_preference() -> None:
    by_aa, _ = load_codon_usage()
    protein = "MKTAYIAKQRQISFVKSHFSRQLEERLRLI"
    cds, decisions = reverse_translate(protein, by_aa)
    assert _translate(cds) == protein
    expect_default = {aa: by_aa[aa][0][0] for aa in set(protein)}
    for d in decisions:
        if not d["note"]:
            assert d["codon"] == expect_default[d["aa"]]
    lys = [d for d in decisions if d["aa"] == "K"]
    assert lys[0]["codon"] == expect_default["K"]


def test_reverse_translation_breaks_homopolymers() -> None:
    by_aa, _ = load_codon_usage()
    assert by_aa["P"][0][0] == "CCC", "P's most-frequent codon must be CCC for this test"
    cds, decisions = reverse_translate("PPPPPP", by_aa)
    assert _translate(cds) == "PPPPPP"
    mitigated = [d for d in decisions if d["note"]]
    assert mitigated, "CC-runs must trigger the documented mitigation"
    assert all("homopolymer mitigation" in d["note"] for d in mitigated)
    runs = homopolymer_runs(cds, 6)
    assert runs == []


def test_codon_table_is_pinned_and_complete() -> None:
    by_aa, digest = load_codon_usage()
    assert len(digest) == 64
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        assert aa in by_aa and by_aa[aa], aa
    counts = [c for aa in by_aa for (_c, c) in by_aa[aa]]
    assert all(c > 0 for c in counts)
    assert "*" in by_aa and by_aa["*"][0][0] == "TGA"


def test_nucleotide_stats_wiring(tmp_path: Path) -> None:
    store = full_store(tmp_path, "work-stats")
    out = run_construct(store, "mini", junction_predictor_factory=fixture_factory)
    stats = out["payload"]["nucleotide_stats"]
    assert stats["cds_length_nt"] == stats["protein_length_aa"] * 3 + 3
    assert stats["restriction_sites"] == {} or all(
        set(v["positions_0based"]) for v in stats["restriction_sites"].values())
    fasta = (store.case_dir("mini") / out["fasta_path"]).read_text()
    assert fasta.startswith(">mini_")
    assert "research design" in fasta.splitlines()[0]
    seq = "".join(fasta.splitlines()[1:])
    assert len(seq) == stats["cds_length_nt"]
