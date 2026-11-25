"""Stage 8 structure tests: pinned published parameters, honest labels.

The stage is deliberately thin arithmetic over a pinned published table, so
the tests pin exactly that: the parameter file's digest, the sanity
classification of extreme sequences (poly-alanine helix, poly-valine sheet,
poly-proline coil), the refusal paths, and determinism (hash-identical
rerun) on a construct produced by the production run_construct.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from neoforge.artifacts import ArtifactStore
from neoforge.errors import ArtifactChainError
from neoforge.stages.structure import (
    load_parameters,
    propensity_labels,
    run_structure,
)

PINNED_PARAMS_SHA256 = "a3b3f9dd27203be45fbfe3709ee78dedc22daccafea70c0d157a07b7005a5d1c"


def _params():
    return load_parameters()


def test_parameter_table_is_pinned() -> None:
    params, digest = _params()
    assert digest == PINNED_PARAMS_SHA256
    assert len(params) == 20
    assert params["A"] == (1.42, 0.83)
    assert params["V"] == (1.06, 1.70)
    assert params["E"] == (1.51, 0.37)


def test_extreme_sequences_classify_sanely() -> None:
    params, _ = _params()
    assert propensity_labels("A" * 12, params)[0] == "H" * 12
    assert propensity_labels("V" * 12, params)[0] == "E" * 12
    assert propensity_labels("P" * 12, params)[0] == "C" * 12


def test_unknown_residue_refuses() -> None:
    params, _ = _params()
    with pytest.raises(ArtifactChainError):
        propensity_labels("AZA", params)


def test_structure_refuses_without_construct(tmp_path: Path) -> None:
    with pytest.raises(ArtifactChainError):
        run_structure(ArtifactStore(tmp_path / "work-empty"), "mini")


def test_structure_artifact_is_deterministic_and_honest(tmp_path: Path) -> None:
    from tests.test_export_ui import _construct_chain

    store = _construct_chain(tmp_path)
    first = run_structure(store, "mini")
    second = run_structure(store, "mini")
    assert first["artifact_id"] == second["artifact_id"]
    assert first["sha256"] == second["sha256"]

    payload = first["payload"]
    protein = payload["protein"]
    labels_str, _table = propensity_labels(protein, load_parameters()[0])
    assert payload["protein_length_aa"] == len(protein)
    assert payload["propensity_table"]["columns"] == ["position", "aa", "pa_mean", "pb_mean", "label"]
    assert payload["propensity_table"]["label"] == list(labels_str)
    assert payload["validation_status"] == "PENDING"
    assert payload["parameters"]["chou_fasman"]["sha256"] == PINNED_PARAMS_SHA256
    counts = payload["summary"]["label_counts"]
    assert sum(counts.values()) == len(protein)
    assert set(counts) <= {"H", "E", "C"}
