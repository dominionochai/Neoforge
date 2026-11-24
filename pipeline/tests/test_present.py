"""Stage 4 wrapper-discipline tests.

present must be a thin, honest wrapper around a prediction model: model used
unmodified, one allele per call, outputs in the mhcflurry_ namespace, paired
reference scoring under identical conditions, and refusals when the chain is
incomplete or an allele is unsupported. A FakePredictor with fixed arithmetic
stands in for MHCflurry so these tests exercise the wrapper contract, not the
model.

Predictions are cross-checked against the translate fixture's exact geometry:
one SNV at CDS[18] (protein position 7 on both transcripts), peptide_length 8
gives windows starting 1..7 per transcript, so 2 transcripts x 7 windows x
2 alleles = 28 prediction rows for 1 distinct variant_key."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from neoforge.artifacts import ArtifactStore, StageRun
from neoforge.errors import ArtifactChainError
from neoforge.sequence import InMemoryGenome
from neoforge.stages.present import run_present
from neoforge.stages.translate import run_translate
from tests.test_translate import CDS, GENOME, _CDS_ON_MINUS, alt_for_both, make_gtf, write_variants

SUPPORTED = ["HLA-A*02:01", "HLA-B*08:01"]


class FakePredictor:
    supported_alleles = SUPPORTED

    def predict(self, peptides, alleles, verbose=0):
        if len(alleles) != 1:
            raise ValueError("wrapper contract: one allele per predict() call")
        return pd.DataFrame({
            "peptide": list(peptides),
            "affinity": [10.0 + sum(map(ord, p)) % 50 for p in peptides],
            "presentation_score": [0.25] * len(peptides),
            "presentation_percentile": [0.5] * len(peptides),
        })


def fake_factory():
    return FakePredictor(), "fake-1.0", "FakePredictor models release test"


VARIANT = ("chrT", 119, CDS[18], alt_for_both(CDS[18:21], 0, _CDS_ON_MINUS[72:75], 2))


def as_rows(payload: dict) -> list[dict]:
    predictions = payload["predictions"]
    cols = predictions["columns"]
    return [dict(zip(cols, vals)) for vals in zip(*[predictions[c] for c in cols])]


def prepared_store(tmp_path: Path, name: str = "work-present",
                   alleles: list[str] | None = None) -> ArtifactStore:
    store = ArtifactStore(tmp_path / name)
    store.write(StageRun(stage="ingest", case="mini"),
                "ingest", {"hla": {"resolved_alleles": alleles if alleles is not None else SUPPORTED}})
    write_variants(store, [VARIANT])
    run_translate(store, "mini", make_gtf(tmp_path, ("T1", "T2")),
                  InMemoryGenome({"chrT": GENOME}), peptide_lengths=(8,))
    return store


def run(tmp_path: Path, factory=fake_factory, **kwargs) -> tuple[ArtifactStore, dict]:
    store = prepared_store(tmp_path)
    return store, run_present(store, "mini", predictor_factory=factory, **kwargs)


def test_predictions_paired_namespaced_and_deterministic(tmp_path: Path) -> None:
    _store, out = run(tmp_path)
    payload = out["payload"]

    n_rows = payload["predictions_made"]
    predicted = as_rows(payload)
    assert n_rows == len(predicted) == 28  # 2 transcripts x 7 windows x 2 alleles
    assert payload["candidates_scored"] == 1

    unique_keys = {(r["variant_key"], r["transcript_id"], r["peptide_length"],
                    r["window_start"], r["allele"]) for r in predicted}
    assert len(unique_keys) == n_rows

    for r in predicted:
        assert r["peptide"] != r["reference_peptide"], "mutant and reference windows collided"
        assert r["mhcflurry_affinity_nm"] == 10.0 + sum(map(ord, r["peptide"])) % 50
        assert r["reference_mhcflurry_affinity_nm"] == 10.0 + sum(map(ord, r["reference_peptide"])) % 50
        assert {r["allele"]} <= set(SUPPORTED)

    evidence = {k: payload[k] for k in ("evidence_kind", "layer", "non_claims")}
    assert evidence["evidence_kind"] == "predicted" and evidence["layer"] == "L3"
    assert any("not measurements" in c for c in evidence["non_claims"])
    assert set(payload["parameters"]["alleles"]) == set(SUPPORTED)
    assert "model_provenance" in payload["parameters"]

    key = VARIANT[0] + f":{VARIANT[1]}{VARIANT[2]}>{VARIANT[3]}"
    affinities = [r["mhcflurry_affinity_nm"] for r in predicted if r["variant_key"] == key]
    assert payload["best_affinity_by_candidate"][key] == min(affinities)

    rerun = run_present(prepared_store(tmp_path, "work-present-rerun"), "mini",
                        predictor_factory=fake_factory)
    assert rerun["sha256"] == out["sha256"], "identical inputs and model must reproduce the artifact"


def test_limit_subset_recorded_honestly(tmp_path: Path) -> None:
    _store, out = run(tmp_path, limit=5)
    payload = out["payload"]
    assert payload["predictions_made"] == 10  # 5 windows x 2 alleles
    assert payload["candidates_scored"] == 1
    note = payload["parameters"]["dev_note"]
    assert "first 5" in note and "of 14" in note
    assert sorted(payload["predictions"]["window_start"]) == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5]
    assert all(a == "HLA-A*02:01" for a in payload["predictions"]["allele"][:5])
    assert sorted({r["allele"] for r in as_rows(payload)}) == sorted(SUPPORTED)


def test_unsupported_allele_refuses_to_predict(tmp_path: Path) -> None:
    store = prepared_store(tmp_path, alleles=[SUPPORTED[0], "HLA-A*99:99"])
    with pytest.raises(ArtifactChainError, match="unsupported"):
        run_present(store, "mini", predictor_factory=fake_factory)


def test_missing_ingest_alleles_block_predictions(tmp_path: Path) -> None:
    store = prepared_store(tmp_path, alleles=[])
    with pytest.raises(ArtifactChainError, match="no HLA alleles"):
        run_present(store, "mini", predictor_factory=fake_factory)


def test_chain_inputs_recorded_in_manifest(tmp_path: Path) -> None:
    store, out = run(tmp_path)
    manifest = store.read_manifest("mini")
    record = manifest["stages"]["present"]
    assert record["inputs"] == {
        "translate": manifest["stages"]["translate"]["artifact_id"],
        "ingest": manifest["stages"]["ingest"]["artifact_id"],
    }
    assert record["artifact_id"] == out["artifact_id"]
