"""Stage 4 — present: MHCflurry prediction wrapper (L3 only).

Every number this stage produces is a *prediction*: stored with the
mhcflurry_ namespace, the allele, the peptide, and the predictor's model
provenance. Nothing here is measured or validated immunogenicity.

Per-allele predictions come from single-allele calls of the published
Class1PresentationPredictor (used unmodified: we wrap, never re-implement),
so each (peptide, allele) pair carries affinity (nM), presentation score and
percentile. Each mutant window's paired reference window is scored under
identical conditions so prioritisation can apply the documented
change-vs-reference rule instead of an invented composite score."""
from __future__ import annotations

import time
from typing import Callable

from neoforge.artifacts import ArtifactStore, StageRun
from neoforge.errors import ArtifactChainError

BATCH = 20000


def default_predictor_factory():
    import mhcflurry  # noqa: PLC0415 — heavy import only needed in this stage

    predictor = mhcflurry.Class1PresentationPredictor.load()
    return predictor, mhcflurry.__version__, str(predictor.provenance_string)


def run_present(store: ArtifactStore, case: str,
                predictor_factory: Callable | None = None,
                limit: int | None = None) -> dict:
    predictor_factory = predictor_factory or default_predictor_factory
    derived = store.load(case, "translate")["derived"]
    cols = derived["columns"]
    hla = store.load(case, "ingest")["hla"]
    alleles = hla["resolved_alleles"]
    if not alleles:
        raise ArtifactChainError("ingest resolved no HLA alleles; predictions must not be invented")

    rows = [dict(zip(cols, vals)) for vals in zip(*[derived[c] for c in cols])]
    note = None
    if limit is not None:
        total = len(rows)
        rows = rows[:limit]
        note = f"dev subset: first {len(rows)} windows of {total}"

    predictor, version, provenance = predictor_factory()
    unsupported = [a for a in alleles if a not in set(predictor.supported_alleles)]
    if unsupported:
        raise ArtifactChainError(
            f"alleles {unsupported} unsupported by MHCflurry models ({provenance[:40]}...); "
            "refusing to predict"
        )

    out = {name: [] for name in [
        "variant_key", "transcript_id", "gene_id", "peptide_length", "window_start",
        "allele", "peptide", "mhcflurry_affinity_nm",
        "mhcflurry_presentation_score", "mhcflurry_presentation_percentile",
        "reference_peptide", "reference_mhcflurry_affinity_nm",
    ]}

    t0 = time.time()
    for base in range(0, len(rows), BATCH):
        chunk = rows[base:base + BATCH]
        specs = [(r, allele) for r in chunk for allele in alleles]
        for allele_idx, allele in enumerate(alleles):
            sel = list(range(allele_idx, len(specs), len(alleles)))
            mutant_df = predictor.predict([specs[i][0]["window_seq"] for i in sel],
                                          alleles=[allele], verbose=0)
            ref_df = predictor.predict([specs[i][0]["ref_window_seq"] for i in sel],
                                       alleles=[allele], verbose=0)
            for local_i, global_i in enumerate(sel):
                r, _ = specs[global_i]
                mutant_row, ref_row = mutant_df.iloc[local_i], ref_df.iloc[local_i]
                out["variant_key"].append(r["variant_key"])
                out["transcript_id"].append(r["transcript_id"])
                out["gene_id"].append(r["gene_id"])
                out["peptide_length"].append(r["peptide_length"])
                out["window_start"].append(r["window_start"])
                out["allele"].append(allele)
                out["peptide"].append(r["window_seq"])
                out["mhcflurry_affinity_nm"].append(float(mutant_row["affinity"]))
                out["mhcflurry_presentation_score"].append(float(mutant_row["presentation_score"]))
                out["mhcflurry_presentation_percentile"].append(float(mutant_row["presentation_percentile"]))
                out["reference_peptide"].append(r["ref_window_seq"])
                out["reference_mhcflurry_affinity_nm"].append(float(ref_row["affinity"]))

    best: dict[str, float] = {}
    for i, key in enumerate(out["variant_key"]):
        best[key] = min(best.get(key, float("inf")), out["mhcflurry_affinity_nm"][i])
    parameters = {
        "predictor": "mhcflurry Class1PresentationPredictor",
        "mhcflurry_version": version,
        "model_provenance": provenance,
        "alleles": alleles,
        "batch_size": BATCH,
        "limit": limit,
        "dev_note": note,
    }
    run = StageRun(
        stage="present",
        case=case,
        inputs={
            "translate": store.read_manifest(case)["stages"]["translate"]["artifact_id"],
            "ingest": store.read_manifest(case)["stages"]["ingest"]["artifact_id"],
        },
        parameters=parameters,
    )
    payload = {
        "parameters": parameters,
        "evidence_kind": "predicted",
        "layer": "L3",
        "non_claims": [
            "MHC binding and presentation scores are predictions, not measurements",
            "no immunogenicity of any candidate is claimed",
        ],
        "candidates_scored": len(best),
        "predictions_made": len(out["variant_key"]),
        "predictions": {"columns": list(out), **out, "uniqueness_key":
                        ["variant_key", "transcript_id", "gene_id", "peptide_length", "window_start", "allele"]},
        "best_affinity_by_candidate": best,
    }
    artifact_id, digest = store.write(run, "present", payload)
    return {"artifact_id": artifact_id, "sha256": digest, "payload": payload}
