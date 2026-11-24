"""Stage 5 — prioritize: explainable tiering of the L3 prediction table.

Tier rules are fixed and published here (SCIENTIFIC_DESIGN §4); there is no
composite or weighted score anywhere in this stage. Gates use only documented
thresholds from the field (pVACtools, TESLA): strong binder < 500 nM best
predicted affinity; expression TPM >= 1 when a source is provided; tumour VAF
>= 0.3 when measured. A gate without its input is recorded as missing with an
explicit reason — it never silently passes as satisfied.

Ordering inside a tier is by ascending best predicted affinity, a single
citable criterion. Everything this stage reads is L1 (measured VAF/depth),
L3 (MHCflurry predictions) and their join; no immunogenicity is implied.
"""
from __future__ import annotations

from neoforge.errors import ArtifactChainError

STRONG_BINDER_NM = 500.0
EXPRESSED_TPM = 1.0
MIN_TUMOUR_VAF = 0.3


def run_prioritize(store, case: str, *,
                   expression: dict[str, float] | None = None,
                   strong_binder_nm: float = STRONG_BINDER_NM,
                   expressed_tpm: float = EXPRESSED_TPM,
                   min_tumour_vaf: float = MIN_TUMOUR_VAF) -> dict:
    from neoforge.artifacts import StageRun  # noqa: PLC0415 — keeps stage dep surface small

    try:
        present = store.load(case, "present")
        selected = store.load(case, "variants")["selected"]
    except FileNotFoundError as exc:
        raise ArtifactChainError(str(exc)) from exc
    if present.get("evidence_kind") != "predicted":
        raise ValueError("prioritize requires the present stage's predicted artifact")
    predictions = present["predictions"]
    cols = predictions["columns"]
    rows = [dict(zip(cols, vals)) for vals in zip(*[predictions[c] for c in cols])]

    tvaf_by_key, depth_by_key = _measured_variant_features(selected)
    expression = expression or {}

    by_candidate: dict[str, list[dict]] = {}
    for r in rows:
        by_candidate.setdefault(r["variant_key"], []).append(r)

    out = {
        "variant_key": [], "tier": [], "prediction_rows": [],
        "best_peptide": [], "best_allele": [],
        "mhcflurry_affinity_nm": [], "mhcflurry_presentation_score": [],
        "mhcflurry_presentation_percentile": [],
        "reference_mhcflurry_affinity_nm": [],
        "delta_affinity_nm_vs_reference": [], "improves_on_reference": [],
        "tvaf": [], "tumor_depth": [], "expression_tpm": [],
        "strong_binder": [], "expressed": [], "passes_vaf_gate": [],
        "gate_flags": [],
    }
    for key in sorted(by_candidate):
        cand_rows = by_candidate[key]
        best = min(cand_rows, key=lambda r: (r["mhcflurry_affinity_nm"], r["peptide"], r["allele"]))
        tvaf = tvaf_by_key.get(key)
        depth = depth_by_key.get(key)
        tpm = expression.get(key)
        flags: list[str] = []
        expressed = _gate(tpm, expressed_tpm, "expression_tpm_unavailable", flags)
        passes_vaf = _gate(tvaf, min_tumour_vaf, "tumour_vaf_unavailable", flags)
        strong = best["mhcflurry_affinity_nm"] < strong_binder_nm

        tier = "A" if strong and expressed is True and passes_vaf is not False else (
            "B" if strong else "C")
        if not strong:
            flags.append("not_strong_binder")

        out["variant_key"].append(key)
        out["tier"].append(tier)
        out["prediction_rows"].append(len(cand_rows))
        out["best_peptide"].append(best["peptide"])
        out["best_allele"].append(best["allele"])
        out["mhcflurry_affinity_nm"].append(best["mhcflurry_affinity_nm"])
        out["mhcflurry_presentation_score"].append(best["mhcflurry_presentation_score"])
        out["mhcflurry_presentation_percentile"].append(best["mhcflurry_presentation_percentile"])
        out["reference_mhcflurry_affinity_nm"].append(best["reference_mhcflurry_affinity_nm"])
        delta = best["mhcflurry_affinity_nm"] - best["reference_mhcflurry_affinity_nm"]
        out["delta_affinity_nm_vs_reference"].append(delta)
        out["improves_on_reference"].append(delta < 0)
        out["tvaf"].append(tvaf)
        out["tumor_depth"].append(depth)
        out["expression_tpm"].append(tpm)
        out["strong_binder"].append(strong)
        out["expressed"].append(expressed)
        out["passes_vaf_gate"].append(passes_vaf)
        out["gate_flags"].append(sorted(flags))

    _order(out)

    tier_counts = {t: out["tier"].count(t) for t in sorted(out["tier"])}
    parameters = {
        "strong_binder_nm": strong_binder_nm,
        "expressed_tpm": expressed_tpm,
        "min_tumour_vaf": min_tumour_vaf,
        "threshold_sources": ["pVACtools strong-binder convention", "TESLA recommended features"],
        "ordering": "tier ascending, then best mhcflurry_affinity_nm ascending",
        "expression_source": "provided" if expression else "none provided",
    }
    payload = {
        "parameters": parameters,
        "evidence_kind": "derived",
        "layer": "L1+L3",
        "non_claims": [
            "tiers rank research priority, not immunogenicity",
            "gates with missing inputs are flagged, never silently passed",
        ],
        "candidates_prioritised": len(out["variant_key"]),
        "tier_counts": tier_counts,
        "candidates": {"columns": list(out), **out},
    }
    run = StageRun(
        stage="prioritize",
        case=case,
        inputs={"present": store.read_manifest(case)["stages"]["present"]["artifact_id"],
                "variants": store.read_manifest(case)["stages"]["variants"]["artifact_id"]},
        parameters=parameters,
    )
    artifact_id, digest = store.write(run, "prioritize", payload)
    return {"artifact_id": artifact_id, "sha256": digest, "payload": payload}


def _gate(value, threshold: float, missing_flag: str, flags: list[str]) -> bool | None:
    if value is None:
        flags.append(missing_flag)
        return None
    return value >= threshold


def _measured_variant_features(measured: dict) -> tuple[dict[str, float], dict[str, float]]:
    tvaf, depth = {}, {}
    for i, key in enumerate(measured["variant_key"]):
        tvaf[key] = measured["tvaf"][i]
        depth[key] = measured["tumor_depth"][i]
    return tvaf, depth


def _order(out: dict) -> None:
    tier_rank = {"A": 0, "B": 1, "C": 2}
    order = sorted(range(len(out["variant_key"])),
                   key=lambda i: (tier_rank[out["tier"][i]], out["mhcflurry_affinity_nm"][i],
                                  out["variant_key"][i]))
    for col in out:
        out[col] = [out[col][i] for i in order]
