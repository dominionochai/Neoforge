"""L4 TESLA harness tests on a synthetic label workbook.

The metric math (TTIF, validated_in_top20, AUPRC) is pinned against
hand-computable orderings; the real benchmark run is separate. The fake
predictor gives affinity = 10 + ord(peptide[0]) - ord('A'), so ordering is
engineerable from the peptide column alone: 'A' peptides score 10.0, 'Z'
peptides 35.0. presentation_score deliberately ranks in the OPPOSITE order
(A=8.0, Z=33.0) so the rule-ablation orders differ by construction."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from neoforge.artifacts import ArtifactStore
from neoforge.evaluation.tesla import load_labels, run_tesla_evaluation


class FakePredictor:
    supported_alleles = ["HLA-A*01:01", "HLA-B*07:02"]

    def predict(self, peptides, alleles, verbose=0):
        assert len(alleles) == 1
        return pd.DataFrame({
            "peptide": list(peptides),
            "affinity": [10.0 + ord(p[0]) - ord("A") for p in peptides],
            "presentation_score": [8.0 + (ord(p[0]) - ord("A")) for p in peptides],
            "presentation_percentile": [0.4] * len(peptides),
        })


def factory():
    return FakePredictor(), "fake-1.0", "FakePredictor models release test"


def write_xlsx(path: Path, rows: list[tuple[str, str, bool]]) -> None:
    pd.DataFrame({
        "PMHC": [f"{a}_{p}" for p, a, _v in rows],
        "PATIENT_ID": [1] * len(rows),
        "TISSUE_TYPE": ["PBMC"] * len(rows),
        "MHC": [a for _p, a, _v in rows],
        "ALT_EPI_SEQ": [p for p, _a, _v in rows],
        "PEP_LEN": [len(p) for p, _a, _v in rows],
        "MEASURED_BINDING_AFFINITY": [50.0] * len(rows),
        "VALIDATED": [v for _p, _a, v in rows],
    }).to_excel(path, sheet_name="master-bindings-selected", index=False)


def test_labels_load_normalises_alleles(tmp_path: Path) -> None:
    xlsx = tmp_path / "labels.xlsx"
    write_xlsx(xlsx, [("APEPTIDE", "A*01:01", True), ("BPEPTIDE", "B*07:02", False)])
    labels = load_labels(xlsx)
    assert labels["alleles"] == ["HLA-A*01:01", "HLA-B*07:02"] and labels["n_patients"] == 1
    assert labels["validated"] == [True, False]
    assert labels["measured_binding_affinity_nm"] == [50.0, 50.0]


def test_metrics_match_hand_computed_ordering(tmp_path: Path) -> None:
    xlsx = tmp_path / "labels.xlsx"
    rows = [("APEPTIDE1", "A*01:01", True), ("APEPTIDE2", "A*01:01", True)]
    rows += [(f"ZPEPTIDE{i}", "A*01:01", False) for i in range(8)]
    write_xlsx(xlsx, rows)
    out = run_tesla_evaluation(ArtifactStore(tmp_path / "work"), xlsx, factory)
    m = out["payload"]["metrics"]
    assert m["n_candidates"] == 10 and m["n_validated"] == 2
    assert m["ttif_top20pct_cut_rank"] == 2
    assert m["ttif_top20pct"] == 1.0 and m["validated_in_top20"] == 2
    assert m["auprc"] == 1.0 and m["auprc_random_expectation"] == pytest.approx(0.2)
    assert m["ttif_oracle"] == 1.0 and m["auprc_oracle"] == 1.0
    assert out["payload"]["rank_rule_falsified"] is False


def test_negative_result_records_falsification(tmp_path: Path) -> None:
    xlsx = tmp_path / "labels.xlsx"
    rows = [("APEPTIDE1", "A*01:01", False), ("APEPTIDE2", "A*01:01", False)]
    rows += [(f"ZPEPTIDE{i}", "B*07:02" if i % 2 else "A*01:01", True) for i in range(10)]
    write_xlsx(xlsx, rows)
    out = run_tesla_evaluation(ArtifactStore(tmp_path / "work"), xlsx, factory)
    m = out["payload"]["metrics"]
    assert m["ttif_top20pct"] < 0.2
    assert out["payload"]["rank_rule_falsified"] is True
    assert m["auprc"] < m["auprc_oracle"]
    worst = m["median_rank_of_validated"]
    assert worst > 6


def test_unsupported_allele_falsification_guard(tmp_path: Path) -> None:
    xlsx = tmp_path / "labels.xlsx"
    write_xlsx(xlsx, [("APEPTIDE", "C*07:02", True)])
    with pytest.raises(ValueError, match="unsupported"):
        run_tesla_evaluation(ArtifactStore(tmp_path / "work"), xlsx, factory)


def test_evaluation_is_deterministic(tmp_path: Path) -> None:
    xlsx = tmp_path / "labels.xlsx"
    rows = [("APEPTIDE", "A*01:01", True)]
    rows += [(f"ZPEPTIDE{i}", "A*01:01", False) for i in range(9)]
    write_xlsx(xlsx, rows)
    one = run_tesla_evaluation(ArtifactStore(tmp_path / "w1"), xlsx, factory)
    two = run_tesla_evaluation(ArtifactStore(tmp_path / "w2"), xlsx, factory)
    assert one["sha256"] == two["sha256"]

def write_reverse_xlsx(path: Path) -> None:
    """Two validated 'A' peptides and eight 'Z' decoys on one allele. The fake
    predictor's presentation_score deliberately ranks Z>peptides highest, so
    the presentation-score ordering is the mirror of the affinity ordering."""
    rows = [("APEPTIDE1", "A*01:01", True), ("APEPTIDE2", "A*01:01", True)]
    rows += [(f"ZPEPTIDE{i}", "A*01:01", False) for i in range(8)]
    write_xlsx(path, rows)


def test_rule_ablation_discriminates_criteria(tmp_path: Path) -> None:
    xlsx = tmp_path / "labels.xlsx"
    write_reverse_xlsx(xlsx)
    out = run_tesla_evaluation(ArtifactStore(tmp_path / "work"), xlsx, factory)
    entries = {e["criterion"].split(" (")[0]: e for e in out["payload"]["rule_ablation"]["entries"]}
    aff = entries["ascending mhcflurry_affinity_nm"]
    pres = entries["descending mhcflurry_presentation_score"]
    assert aff["ttif_top20pct"] == 1.0 and pres["ttif_top20pct"] == 0.0
    assert aff["auprc"] == 1.0 and pres["auprc"] == pytest.approx((1 / 2) * (1 / 9) + (1 / 2) * (2 / 10))
    assert aff["median_rank_of_validated"] < pres["median_rank_of_validated"]
    assert out["payload"]["metrics"]["ttif_top20pct"] == 1.0


def test_rank_curve_matches_hand_computed_points(tmp_path: Path) -> None:
    xlsx = tmp_path / "labels.xlsx"
    write_reverse_xlsx(xlsx)
    out = run_tesla_evaluation(ArtifactStore(tmp_path / "work"), xlsx, factory)
    curve = out["payload"]["metrics"]["rank_at_k_curve"]
    assert curve["k_pct"] == [1, 2, 5, 10, 20, 30, 50]
    assert curve["cut_rank"] == [0, 0, 0, 1, 2, 3, 5]
    assert curve["found_fraction"][0] == 0.0
    assert curve["found_fraction"][4] == 1.0 and curve["found_fraction"][6] == 1.0
    assert curve["random_expectation"] == [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5]


def test_ablation_recorded_unavailable_without_presentation_scores(tmp_path: Path) -> None:
    class AffinityOnlyPredictor:
        supported_alleles = ["HLA-A*01:01"]

        def predict(self, peptides, alleles, verbose=0):
            return pd.DataFrame({
                "peptide": list(peptides),
                "affinity": [10.0] * len(peptides),
            })

    xlsx = tmp_path / "labels.xlsx"
    write_reverse_xlsx(xlsx)
    out = run_tesla_evaluation(ArtifactStore(tmp_path / "work"), xlsx,
                               lambda: (AffinityOnlyPredictor(), "aff-only-1.0", "AffinityOnly"))
    ablation = out["payload"]["rule_ablation"]
    assert ablation["entries"] == []
    assert "did not provide presentation_score" in ablation["note"]
    assert out["payload"]["predictions"]["mhcflurry_presentation_score"][0] is None
