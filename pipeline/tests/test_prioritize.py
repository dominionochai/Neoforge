"""Stage 5 tier-rule tests.

The tier rules are the product: strong binder < 500 nM, TPM >= 1 when a source
is provided, VAF >= 0.3 when measured — gates only, no composite score. These
tests pin the gate semantics that protect against the classic mistake of
letting a missing input silently pass as satisfied.

Chain fixtures reuse test_present's mini-genome translate products; a
configurable fake predictor stands in for MHCflurry so tier boundaries are
exercised with known affinities (one run of the real chain showed 329
candidates: 155 strong binders, 0 Tier A because no expression source exists
for this case — asserted separately in the real-case checks)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from neoforge.artifacts import ArtifactStore, StageRun
from neoforge.errors import ArtifactChainError
from neoforge.sequence import InMemoryGenome
from neoforge.stages.prioritize import run_prioritize
from neoforge.stages.present import run_present
from neoforge.stages.translate import run_translate
from tests.test_translate import CDS, GENOME, _CDS_ON_MINUS, alt_for_both, make_gtf, write_variants

SUPPORTED = ["HLA-A*02:01", "HLA-B*08:01"]
VARIANT = ("chrT", 119, CDS[18], alt_for_both(CDS[18:21], 0, _CDS_ON_MINUS[72:75], 2))


def make_factory(affinity_of) -> callable:
    class FakePredictor:
        supported_alleles = SUPPORTED

        def predict(self, peptides, alleles, verbose=0):
            assert len(alleles) == 1
            return pd.DataFrame({
                "peptide": list(peptides),
                "affinity": [affinity_of(p) for p in peptides],
                "presentation_score": [0.4] * len(peptides),
                "presentation_percentile": [0.5] * len(peptides),
            })

    return lambda: (FakePredictor(), "fake-1.0", "FakePredictor models release test")


def as_rows(payload: dict) -> list[dict]:
    candidates = payload["candidates"]
    cols = candidates["columns"]
    return [dict(zip(cols, vals)) for vals in zip(*[candidates[c] for c in cols])]


def chain(tmp_path: Path, factory, name: str = "work-prio", tvaf=None,
          variants=(), peptide_lengths=(8,)) -> ArtifactStore:
    store = ArtifactStore(tmp_path / name)
    store.write(StageRun(stage="ingest", case="mini"),
                "ingest", {"hla": {"resolved_alleles": SUPPORTED}})
    rows = [VARIANT, *variants]
    write_variants(store, rows)
    if tvaf is not None:
        store_payload = store.load("mini", "variants")
        store_payload["selected"]["tvaf"] = [tvaf] * len(rows)
        store.write(StageRun(stage="variants", case="mini"), "variants", store_payload)
    run_translate(store, "mini", make_gtf(tmp_path, ("T1", "T2")),
                  InMemoryGenome({"chrT": GENOME}), peptide_lengths=peptide_lengths)
    run_present(store, "mini", predictor_factory=factory)
    return store


def test_strong_binder_without_expression_lands_in_tier_b(tmp_path: Path) -> None:
    store = chain(tmp_path, make_factory(lambda p: 120.0))
    out = run_prioritize(store, "mini")
    row = as_rows(out["payload"])[0]
    assert out["payload"]["tier_counts"] == {"B": 1}
    assert row["strong_binder"] is True and row["tier"] == "B"
    assert row["expressed"] is None and row["passes_vaf_gate"] is None
    assert "expression_tpm_unavailable" in row["gate_flags"]
    assert "tumour_vaf_unavailable" in row["gate_flags"]
    assert row["improves_on_reference"] is (row["delta_affinity_nm_vs_reference"] < 0)
    assert row["mhcflurry_affinity_nm"] == min(
        r["mhcflurry_affinity_nm"] for r in as_rows(out["payload"]))


def test_weak_binder_lands_in_tier_c(tmp_path: Path) -> None:
    store = chain(tmp_path, make_factory(lambda p: 800.0))
    out = run_prioritize(store, "mini")
    row = as_rows(out["payload"])[0]
    assert out["payload"]["tier_counts"] == {"C": 1}
    assert "not_strong_binder" in row["gate_flags"]
    assert row["strong_binder"] is False


def test_tier_a_requires_expression_and_passing_vafs(tmp_path: Path) -> None:
    factory = make_factory(lambda p: 120.0)

    out = run_prioritize(chain(tmp_path, factory, tvaf=0.55), "mini",
                         expression={"G": 5.0})
    assert out["payload"]["tier_counts"] == {"A": 1}

    out = run_prioritize(chain(tmp_path, factory, name="work-prio-tvaf-low", tvaf=0.25),
                         "mini", expression={"G": 5.0})
    assert out["payload"]["tier_counts"] == {"B": 1}
    assert as_rows(out["payload"])[0]["passes_vaf_gate"] is False

    out = run_prioritize(chain(tmp_path, factory, name="work-prio-tvaf-missing", tvaf=None), "mini",
                         expression={"G": 5.0})

    gene_absent = run_prioritize(chain(tmp_path, factory, tvaf=0.55), "mini",
                                 expression={"OTHER": 9.0})
    rows = as_rows(gene_absent["payload"])[0]
    assert rows["expressed"] is None and rows["tier"] == "B"
    assert "expression_gene_absent_from_measured_source" in rows["gate_flags"]
    assert out["payload"]["tier_counts"] == {"A": 1}
    assert "tumour_vaf_unavailable" in as_rows(out["payload"])[0]["gate_flags"]


def test_ordering_tiers_first_then_affinity(tmp_path: Path) -> None:
    second = ("chrT", 112, CDS[11], alt_for_both(CDS[9:12], 1, _CDS_ON_MINUS[81:84], 1))
    factory = make_factory(lambda p: 100.0 + (sum(map(ord, p)) % 40) * 10)
    store = chain(tmp_path, factory, variants=[second])
    out = run_prioritize(store, "mini", strong_binder_nm=250.0)
    rows = as_rows(out["payload"])
    assert len(rows) == 2
    tier_rank = {"A": 0, "B": 1, "C": 2}
    assert [tier_rank[r["tier"]] for r in rows] == sorted(tier_rank[r["tier"]] for r in rows)
    for a, b in zip(rows, rows[1:]):
        if tier_rank[a["tier"]] == tier_rank[b["tier"]]:
            assert a["mhcflurry_affinity_nm"] <= b["mhcflurry_affinity_nm"]
            assert a["mhcflurry_affinity_nm"] < 250 or b["mhcflurry_affinity_nm"] >= 250


def test_prioritize_is_deterministic(tmp_path: Path) -> None:
    second = ("chrT", 112, CDS[11], alt_for_both(CDS[9:12], 1, _CDS_ON_MINUS[81:84], 1))
    factory = make_factory(lambda p: 100.0 + (sum(map(ord, p)) % 40) * 10)
    out = run_prioritize(chain(tmp_path, factory, variants=[second]), "mini")
    rerun = run_prioritize(chain(tmp_path, factory, name="work-prio-rerun", variants=[second]),
                           "mini")
    assert rerun["sha256"] == out["sha256"]


def test_missing_present_artifact_raises_chain_error(tmp_path: Path) -> None:
    store = ArtifactStore(tmp_path / "work-empty")
    with pytest.raises(ArtifactChainError):
        run_prioritize(store, "mini")
