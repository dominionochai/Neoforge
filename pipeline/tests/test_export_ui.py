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
from tests.test_report import SUPPORTED, full_chain


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


# ---------------- construct projection ----------------

def _ingest_for_construct_export(store: ArtifactStore) -> None:
    store.write(
        StageRun(stage="ingest", case="mini"), "ingest",
        {"case": {"case_id": "mini", "display_name": "Mini", "genome_build": "GRCh38",
                  "pairing": {"tumour": "miniT", "normal": "miniN"}},
         "checks": [], "inputs": [{"name": "demo", "path": "demo.vcf.gz",
                                   "sha256": "0" * 64, "role": "demo"}],
         "vcf": {"somatic_snv_vcf": {"counters": {"records": 3, "with_tvaf": 3, "with_depth": 3}},
                 "somatic_indel_vcf": {"counters": {"records": 0}}},
         "variant_tables": {"somatic_snv_vcf": {"chrom": ["chrT"] * 3, "pos": [100, 200, 300],
                                                "ref": ["A"] * 3, "alt": ["G"] * 3,
                                                "tvaf": [0.5] * 3, "nvaf": [0.0] * 3,
                                                "tumor_depth": [100] * 3}},
         "hla": {"resolved_alleles": SUPPORTED, "cell_line": "MINI",
                 "source": {"path": "hla.csv", "sha256": "0" * 64}},
         "expression": {"tpm_by_gene": {"ENSG1": 5.0, "ENSG2": 5.0, "ENSG3": 5.0},
                        "source_description": "fixture measured TPM"}})


def _construct_chain(tmp_path) -> ArtifactStore:
    from tests.test_construct import (ALLELES, GENES_1_3, fixture_factory, write_present,
                                      write_translate, write_variants_fixture)
    from neoforge.stages.construct import run_construct
    from neoforge.stages.prioritize import run_prioritize

    store = ArtifactStore(tmp_path / "work-ui-construct")
    _ingest_for_construct_export(store)
    keys = ["chrT:100A>G", "chrT:200C>T", "chrT:300G>A"]
    write_variants_fixture(store, keys)
    write_translate(store, [(k, f"ENST{i + 1}", g, 40 * i)
                            for i, (k, g) in enumerate(zip(keys, GENES_1_3))])
    write_present(store, [
        {"variant_key": keys[0], "transcript_id": "ENST1", "gene_id": GENES_1_3[0],
         "peptide_length": 9, "window_start": 0, "allele": ALLELES[0],
         "peptide": "SYFKPEPTI", "mhcflurry_affinity_nm": 120.0,
         "mhcflurry_presentation_score": 0.5, "mhcflurry_presentation_percentile": 1.0,
         "reference_peptide": "SYFKPEPTE", "reference_mhcflurry_affinity_nm": 800.0},
        {"variant_key": keys[1], "transcript_id": "ENST2", "gene_id": GENES_1_3[1],
         "peptide_length": 9, "window_start": 0, "allele": ALLELES[0],
         "peptide": "TEPTQVAAA", "mhcflurry_affinity_nm": 80.0,
         "mhcflurry_presentation_score": 0.5, "mhcflurry_presentation_percentile": 1.0,
         "reference_peptide": "TEPTQVAAE", "reference_mhcflurry_affinity_nm": 90.0},
        {"variant_key": keys[2], "transcript_id": "ENST3", "gene_id": GENES_1_3[2],
         "peptide_length": 9, "window_start": 0, "allele": ALLELES[1],
         "peptide": "IMMTKEYLL", "mhcflurry_affinity_nm": 200.0,
         "mhcflurry_presentation_score": 0.5, "mhcflurry_presentation_percentile": 1.0,
         "reference_peptide": "IMMTKEYLE", "reference_mhcflurry_affinity_nm": 900.0},
    ])
    pre = store.load("mini", "present")
    pre["parameters"] = {"binder_nm": 500.0, "predictor": "mhcflurry Class1PresentationPredictor",
                         "mhcflurry_version": "fixture-1.0", "model_provenance": "fixture",
                         "alleles": ALLELES}
    pre["non_claims"] = ["fixture non-claim"]
    store.write(StageRun(stage="present", case="mini"), "present", pre)
    run_prioritize(store, "mini")
    run_construct(store, "mini", junction_predictor_factory=fixture_factory)
    return store


def test_export_projects_construct_block(tmp_path) -> None:
    store = _construct_chain(tmp_path)
    cases_dir = tmp_path / "cases" / "mini"
    cases_dir.mkdir(parents=True)
    (cases_dir / "case.json").write_text('{"case_id": "mini", "display_name": "Mini"}')

    out = export_case(store, "mini", tmp_path / "cases")
    con = out["construct"]
    rec = out["pipeline_stage_records"]["construct"]
    assert set(out["pipeline_stage_records"]) >= {"ingest", "variants", "translate",
                                                  "present", "prioritize", "construct"}
    assert rec["artifact_id"] == store.read_manifest("mini")["stages"]["construct"]["artifact_id"]
    assert con["source"] == f"construct.{rec['artifact_id']}.{rec['sha256']}"
    assert con["evidence_kind"] == "in_silico_design"
    assert con["stamp"].startswith("research design")
    assert [e["variant_key"] for e in con["selection_ledger_included"]] == \
        ["chrT:100A>G", "chrT:300G>A"]
    assert con["exclusion_reason_counts"] == {"no_mutant_vs_wild_type_discrimination": 1}
    assert con["junction_audit"]["final_rows"]
    assert con["sequence_checks"]["translation_matches_design"] is True
    assert con["nucleotide_map"][0]["kind"] == "placeholder"
    assert con["fasta_file"] == f"construct_cds.{rec['artifact_id']}.fasta"
    assert all("reason" in e for e in con["selection_ledger_excluded"])


def test_export_construct_is_none_until_construct_runs(mini_export) -> None:
    _store, payload = mini_export
    assert payload["construct"] is None
    assert "construct" not in payload["pipeline_stage_records"]


def test_export_fasta_projection_matches_store_copy(tmp_path) -> None:
    import hashlib

    store = _construct_chain(tmp_path)
    cases_dir = tmp_path / "cases" / "mini"
    cases_dir.mkdir(parents=True)
    (cases_dir / "case.json").write_text('{"case_id": "mini"}')
    out_dir = tmp_path / "ui"
    out_dir.mkdir()

    import neoforge.export_ui as ex
    payload = ex.export_case(store, "mini", tmp_path / "cases")
    copied = ex._copy_construct_fastas(store, "mini", out_dir, payload)
    assert copied == [payload["construct"]["fasta_file"]]
    fasta = out_dir / copied[0]
    assert hashlib.sha256(fasta.read_bytes()).hexdigest() == payload["construct"]["fasta_sha256"]
    body = fasta.read_text().splitlines()
    assert body[0].startswith(">mini_") and "research design" in body[0]
    assert len("".join(body[1:])) == payload["construct"]["nucleotide_stats"]["cds_length_nt"]
