"""UI-projection tests on the synthetic mini-genome chain.

The export module is the contract between the artifact store and the frontend:
it must be a pure projection (no new computation), keep the upstream artifact
digests on every block, and scope heavy tables to the shortlist. The tesla
projection is tested against a hand-built store-side benchmark payload — the
harness itself is covered by test_tesla_eval."""
from __future__ import annotations

import pytest

from neoforge.artifacts import ArtifactStore, StageRun
from neoforge.errors import ArtifactChainError
from neoforge.export_ui import export_case, export_tesla
from neoforge.stages.report import run_report
from tests.test_report import full_chain


def _augment_ingest_with_snv_table(store: ArtifactStore) -> None:
    ing = store.load("mini", "ingest")
    ing["variant_tables"] = {
        "somatic_snv_vcf": {
            "chrom": ["chrT"], "pos": [119], "ref": ["G"], "alt": ["T"],
            "tvaf": [0.4], "nvaf": [0.0], "tumor_depth": [50],
        },
    }
    store.write(StageRun(stage="ingest", case="mini"), "ingest", ing)


@pytest.fixture()
def mini_export(tmp_path):
    store = full_chain(tmp_path, "work-export")
    _augment_ingest_with_snv_table(store)
    run_report(store, "mini")
    cases_dir = tmp_path / "cases" / "mini"
    cases_dir.mkdir(parents=True)
    (cases_dir / "case.json").write_text(
        '{"case_id": "mini", "display_name": "Mini", "description": "fixture"}')
    return store, export_case(store, "mini", tmp_path / "cases")


def test_export_carries_upstream_artifact_hashes(mini_export) -> None:
    store, payload = mini_export
    for stage, rec in payload["pipeline_stage_records"].items():
        manifest_rec = store.read_manifest("mini")["stages"][stage]
        assert rec["artifact_id"] == manifest_rec["artifact_id"]
        assert rec["sha256"] == manifest_rec["sha256"]
    assert set(payload["pipeline_stage_records"]) == {
        "ingest", "variants", "translate", "present", "prioritize", "report"}


def test_export_ingest_projection_is_curated(mini_export) -> None:
    _store, payload = mini_export
    ing = payload["ingest"]
    case_id = ing["case"] if isinstance(ing["case"], str) else ing["case"]["case_id"]
    assert case_id == "mini"
    assert ing["variant_counts"]["somatic_snv_positions"] == 1
    assert ing["snv_table"]["columns"] == ["chrom", "pos", "ref", "alt", "tvaf", "nvaf", "tumor_depth"]
    assert ing["snv_table"]["tvaf"] == [0.4]
    assert ing["vaf_source"]["with_tvaf"] == 1
    assert "tpm_by_gene" not in (ing["expression"] or {})


def test_export_scopes_heavy_tables_with_declared_scope(mini_export) -> None:
    _store, payload = mini_export
    preds = payload["present"]["predictions"]
    assert preds["rows_exported"] == preds["rows_total"] > 0
    assert preds["export_scope"] == "prediction rows of the report shortlist variants"
    tra = payload["translate"]
    assert tra["derived_rows_exported"] == tra["derived_rows_total"] > 0
    assert tra["export_scope"] == "windows of shortlisted/prioritised variants only"
    assert set(payload["report"]["shortlist"][0]) >= {"variant_key", "tier", "best_peptide"}


def test_export_refuses_a_case_without_artifacts(tmp_path) -> None:
    store = ArtifactStore(tmp_path / "work-empty")
    (tmp_path / "cases" / "ghost").mkdir(parents=True)
    (tmp_path / "cases" / "ghost" / "case.json").write_text("{}")
    with pytest.raises(ArtifactChainError, match="no ingest artifact"):
        export_case(store, "ghost", tmp_path / "cases")


def test_export_tesla_projects_validated_rows(tmp_path) -> None:
    store = ArtifactStore(tmp_path / "work-tesla")
    payload = {
        "parameters": {"falsification_clause": "ttif_top20pct <= 0.2"},
        "evidence_kind": "measured",
        "layer": "L4",
        "non_claims": ["cohort-level evidence only"],
        "rank_rule_falsified": False,
        "metrics": {"ttif_top20pct": 0.5, "auprc": 0.2, "n_candidates": 2,
                    "ttif_random_expectation": 0.2, "auprc_random_expectation": 0.2},
        "predictions": {
            "columns": ["pmhc_id", "peptide", "allele", "rank", "validated", "unexported_col",
                        "mhcflurry_affinity_nm", "measured_binding_affinity_nm"],
            "pmhc_id": ["HLA-A*02:01_AAALENLKT", "HLA-A*02:01_TTTLENIKK", "HLA-A*02:01_GGGLENLKT"],
            "peptide": ["AAALENLKT", "TTTLENIKK", "GGGLENLKT"],
            "allele": ["HLA-A*02:01"] * 3,
            "rank": [3, 1, 2],
            "validated": [True, False, True],
            "unexported_col": ["x", "y", "z"],
            "mhcflurry_affinity_nm": [10.0, 12.0, 11.0],
            "measured_binding_affinity_nm": [9.0, None, 8.0],
        },
    }
    store.write(StageRun(stage="tesla_eval", case="tesla"), "tesla_eval", payload)

    out = export_tesla(store)
    assert out["layer"] == "L4" and out["rank_rule_falsified"] is False
    pred = out["validated_predictions"]
    assert pred["peptide"] == ["AAALENLKT", "GGGLENLKT"]
    assert pred["measured_binding_affinity_nm"] == [9.0, 8.0]
    assert out["validated_predictions"]["columns_dropped"] == ["unexported_col"]
    assert out["validated_pmhc_ids"] == ["HLA-A*02:01_AAALENLKT", "HLA-A*02:01_GGGLENLKT"]
    assert out["source"].startswith("tesla.tesla_eval.")
