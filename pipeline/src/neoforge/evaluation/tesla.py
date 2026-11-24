"""L4 — external benchmark evaluation against TESLA (Wells et al., Cell 2020).

Labels: `mmc4.xlsx` sheet `master-bindings-selected` — 608 candidate pMHC from
6 melanoma patients, 13 class-I alleles, of which 37 carry an experimental
VALIDATED=True label (tetramer / T-cell assays). File and SHA-256 pinned in
DATA_PROVENANCE.md. This module never touches project cases: it is a
cohort-level external audit of the RANKING RULE, not of any patient.

Method (documented exactly in the artifact it writes):
- every (peptide, allele) pair is scored with the same unmodified per-allele
  MHCflurry wrapper the present stage uses — mhcflurry_affinity_nm only;
- ranking = ascending predicted affinity (single citable criterion, ties
  broken deterministically by peptide then allele);
- metrics, each with its analytic random baseline (no invented simulations):
  * ttif_top20pct: fraction of the validated inside the top 20% of the ranked
    list (TESLA's top-immunogenic-fraction convention); random ranking
    expectation = 0.2 exactly by construction;
  * validated_in_top20: count of validated with rank <= 20; random
    expectation = n_validated * 20 / n;
  * auprc: area under the precision-recall curve, candidates consumed in
    ascending-affinity order, label = VALIDATED (step rule, the sklearn
    average_precision_score convention); random expectation = prevalence
    = n_validated / n;
  * oracle parallels (validated consumed first, independent of workbook row
    order) give the headroom of the label set itself.
Falsification clause (EVALUATION_PLAN): ttif_top20pct <= 0.2 falsifies the
prioritisation design; the artifact records the verdict either way.
"""
from __future__ import annotations

_LABEL_SHEET = "master-bindings-selected"


def load_labels(xlsx_path) -> dict:
    import pandas as pd  # noqa: PLC0415 — heavy import, xlsx only

    df = pd.read_excel(xlsx_path, sheet_name=_LABEL_SHEET)
    required = {"PMHC", "PATIENT_ID", "MHC", "ALT_EPI_SEQ", "VALIDATED"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{xlsx_path}: label sheet missing columns {sorted(missing)}")
    return {
        "columns": ["pmhc_id", "patient_id", "allele", "peptide", "validated",
                    "measured_binding_affinity_nm"],
        "pmhc_id": [str(x) for x in df["PMHC"]],
        "patient_id": [int(x) for x in df["PATIENT_ID"]],
        "allele": ["HLA-" + str(a) for a in df["MHC"]],
        "peptide": [str(x) for x in df["ALT_EPI_SEQ"]],
        "validated": [bool(x) for x in df["VALIDATED"]],
        "measured_binding_affinity_nm": [None if v != v else float(v)
                                         for v in df["MEASURED_BINDING_AFFINITY"]],
        "alleles": sorted({("HLA-" + str(a)) for a in df["MHC"]}),
        "n_patients": int(df["PATIENT_ID"].nunique()),
    }


def _auprc(labels: list[bool], consumption_order: list[int]) -> float:
    positives = sum(labels)
    seen, tp, area, prev_recall = 0, 0, 0.0, 0.0
    for idx in consumption_order:
        seen += 1
        if labels[idx]:
            tp += 1
            recall = tp / positives
            area += (recall - prev_recall) * (tp / seen)
            prev_recall = recall
    return area


def _predict_affinity_by_pair(predictor, pairs: list[tuple[str, str]],
                              alleles: list[str]) -> dict[tuple[str, str], float]:
    by_pair: dict[tuple[str, str], float] = {}
    for allele in alleles:
        peptides = [p for p, a in pairs if a == allele]
        if not peptides:
            continue
        df = predictor.predict(peptides, alleles=[allele], verbose=0)
        for peptide, affinity in zip(peptides, df["affinity"]):
            by_pair[(peptide, allele)] = float(affinity)
    return by_pair


def run_tesla_evaluation(store, xlsx_path, predictor_factory=None) -> dict:
    from neoforge.artifacts import StageRun  # noqa: PLC0415
    from neoforge.stages.present import default_predictor_factory  # noqa: PLC0415

    labels = load_labels(xlsx_path)
    predictor, version, provenance = (predictor_factory or default_predictor_factory)()
    supported = set(predictor.supported_alleles)
    unsupported = [a for a in labels["alleles"] if a not in supported]
    if unsupported:
        raise ValueError(
            f"TESLA alleles {unsupported} unsupported by the loaded model ({provenance[:40]}...); "
            "refusing to evaluate with imputed predictions")

    n = len(labels["peptide"])
    pairs = list(zip(labels["peptide"], labels["allele"]))
    affinity_by_pair = _predict_affinity_by_pair(predictor, pairs, labels["alleles"])
    affinity = [affinity_by_pair[pair] for pair in pairs]

    order = sorted(range(n), key=lambda i: (affinity[i], labels["peptide"][i], labels["allele"][i]))
    rank = [0] * n
    for position, i in enumerate(order, start=1):
        rank[i] = position

    validated = labels["validated"]
    n_validated = sum(validated)
    top20pct_cut = (20 * n) // 100
    metrics = {
        "n_candidates": n,
        "n_validated": n_validated,
        "n_patients": labels["n_patients"],
        "n_alleles": len(labels["alleles"]),
        "ttif_top20pct": sum(1 for i in range(n) if validated[i] and rank[i] <= top20pct_cut) / n_validated,
        "ttif_top20pct_cut_rank": top20pct_cut,
        "ttif_random_expectation": 0.2,
        "ttif_oracle": min(1.0, top20pct_cut / n_validated),
        "validated_in_top20": sum(1 for i in range(n) if validated[i] and rank[i] <= 20),
        "validated_in_top20_random_expectation": n_validated * 20 / n,
        "auprc": _auprc(validated, order),
        "auprc_random_expectation": n_validated / n,
        "auprc_oracle": _auprc(validated, sorted(range(n), key=lambda i: not validated[i])),
        "median_rank_of_validated": sorted(rank[i] for i in range(n) if validated[i])[n_validated // 2],
    }
    falsified = not metrics["ttif_top20pct"] > 0.2
    payload = {
        "parameters": {
            "labels": f"{xlsx_path}#{_LABEL_SHEET}",
            "predictor": "mhcflurry Class1PresentationPredictor",
            "mhcflurry_version": version,
            "model_provenance": provenance,
            "ranking": "ascending mhcflurry_affinity_nm; ties by peptide then allele",
            "allele_normalisation": "TESLA 'A*01:01' -> 'HLA-A*01:01' (documented prefix)",
            "falsification_clause": "ttif_top20pct <= 0.2 falsifies the ranking design",
        },
        "evidence_kind": "measured",
        "layer": "L4",
        "non_claims": [
            "benchmark evidence about the ranking rule at cohort level — not about any project case",
            "mhcflurry values are predictions; VALIDATED labels are external experimental results",
        ],
        "metrics": metrics,
        "rank_rule_falsified": falsified,
        "labels": labels,
        "predictions": {
            "columns": ["pmhc_id", "peptide", "allele", "mhcflurry_affinity_nm", "rank",
                        "validated", "measured_binding_affinity_nm"],
            "pmhc_id": labels["pmhc_id"],
            "peptide": labels["peptide"],
            "allele": labels["allele"],
            "mhcflurry_affinity_nm": affinity,
            "rank": rank,
            "validated": validated,
            "measured_binding_affinity_nm": labels["measured_binding_affinity_nm"],
        },
    }
    run = StageRun(stage="tesla_eval", case="tesla",
                   inputs={"labels": str(xlsx_path)},
                   parameters=payload["parameters"])
    artifact_id, digest = store.write(run, "tesla_eval", payload)
    return {"artifact_id": artifact_id, "sha256": digest, "payload": payload}
