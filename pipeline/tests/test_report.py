"""Stage 6 report tests on the synthetic mini-genome chain.

The report must be an honest projection of the artifact store: no new numbers,
no wall-clock time, warnings and missing gates lifted verbatim, and an
exportable Markdown rendering derived from the same payload. Tests build the
full five-stage chain with the fake predictor from test_prioritize."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from neoforge.artifacts import ArtifactStore, StageRun
from neoforge.errors import ArtifactChainError
from neoforge.sequence import InMemoryGenome
from neoforge.stages.prioritize import STRONG_BINDER_NM, run_prioritize
from neoforge.stages.present import run_present
from neoforge.stages.report import run_report
from neoforge.stages.translate import run_translate
from tests.test_translate import CDS, GENOME, _CDS_ON_MINUS, alt_for_both, make_gtf, write_variants

SUPPORTED = ["HLA-A*02:01", "HLA-B*08:01"]
VARIANT = ("chrT", 119, CDS[18], alt_for_both(CDS[18:21], 0, _CDS_ON_MINUS[72:75], 2))


def make_factory(affinity_of, supported=SUPPORTED):
    class FakePredictor:
        supported_alleles = supported

        def predict(self, peptides, alleles, verbose=0):
            assert len(alleles) == 1
            return pd.DataFrame({
                "peptide": list(peptides),
                "affinity": [affinity_of(p) for p in peptides],
                "presentation_score": [0.4] * len(peptides),
                "presentation_percentile": [0.5] * len(peptides),
            })

    return lambda: (FakePredictor(), "fake-1.0", "FakePredictor models release test")


def full_chain(tmp_path: Path, name: str, affinity=120.0, ingest_checks=()) -> ArtifactStore:
    store = ArtifactStore(tmp_path / name)
    store.write(StageRun(stage="ingest", case="mini"),
                "ingest", {"case": {"case_id": "mini", "display_name": "Mini fixture",
                                    "genome_build": "GRCh38",
                                    "pairing": {"tumour": "miniT", "normal": "miniN"}},
                           "checks": list(ingest_checks),
                           "inputs": [{"name": "demo", "path": "demo.vcf.gz",
                                       "sha256": "0" * 64, "role": "demo"}],
                           "vcf": {"somatic_snv_vcf": {"counters": {"records": len([VARIANT]),
                                                                    "with_tvaf": 1,
                                                                    "with_depth": 1}},
                                   "somatic_indel_vcf": {"counters": {"records": 0}}},
                           "hla": {"resolved_alleles": SUPPORTED, "cell_line": "MINI",
                                   "source": {"path": "hla.csv", "sha256": "0" * 64}}})
    write_variants(store, [VARIANT])
    run_translate(store, "mini", make_gtf(tmp_path, ("T1", "T2")),
                  InMemoryGenome({"chrT": GENOME}), peptide_lengths=(8,))
    run_present(store, "mini", predictor_factory=make_factory(lambda p: affinity))
    run_prioritize(store, "mini")
    return store


def test_report_projects_the_chain_honestly(tmp_path: Path) -> None:
    store = full_chain(tmp_path, "work-report", ingest_checks=[
        {"id": "demo.warn_check", "status": "warn", "detail": "the honest detail"}])
    out = run_report(store, "mini")
    p = out["payload"]

    assert len(p["shortlist"]) == 1
    row = p["shortlist"][0]
    assert row["tier"] == "B" and row["mhcflurry_affinity_nm"] < STRONG_BINDER_NM
    assert row["gate_flags"] == ["expression_tpm_unavailable", "tumour_vaf_unavailable"]

    ledger = p["evidence_ledger"]
    assert ledger["L1_measured"]["kind"] == "measured"
    assert ledger["L2_derived"]["kind"] == "derived"
    assert ledger["L3_predicted"]["kind"] == "predicted"
    assert ledger["L4_experimental"] == {
        "kind": "missing", "status": "PENDING",
        "note": "no TESLA or other experimental immunogenicity labels are evaluated "
                "in this run; no experimental immunogenicity claim exists for any candidate"}

    sources = [f["source"] for f in p["uncertainty_register"]]
    assert "ingest.demo.warn_check" in sources
    assert "prioritize.expression_tpm_unavailable" in sources
    assert "L4_experimental" in sources
    warn = next(f for f in p["uncertainty_register"]
                if f["source"] == "ingest.demo.warn_check")
    assert warn["detail"] == "the honest detail"

    assert p["molecular_profile"]["peptide_windows"] > 1
    assert all(s in p["provenance_chain"]
               for s in ("ingest", "variants", "translate", "present", "prioritize"))

    md = p["research_report_md"]
    assert "research summary, not a treatment recommendation" in md
    assert "Mini fixture" in md
    assert row["best_peptide"] in md
    assert "Pred. affinity (nM)" in md


def test_report_deterministic_and_artifact_chain_pinned(tmp_path: Path) -> None:
    out = run_report(full_chain(tmp_path, "work-report"), "mini")
    rerun = run_report(full_chain(tmp_path, "work-report-rerun"), "mini")
    assert rerun["sha256"] == out["sha256"]
    a, b = out["payload"]["provenance_chain"], rerun["payload"]["provenance_chain"]
    assert a == b


def test_report_refuses_incomplete_chain(tmp_path: Path) -> None:
    store = ArtifactStore(tmp_path / "work-partial")
    store.write(StageRun(stage="ingest", case="mini"), "ingest", {"hla": {"resolved_alleles": SUPPORTED}})
    with pytest.raises(ArtifactChainError, match="requires completed stages"):
        run_report(store, "mini")
