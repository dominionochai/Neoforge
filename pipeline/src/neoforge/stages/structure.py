"""Stage 8 - structure: sequence-derived secondary-structure propensity.

Consumes the construct artifact (Stage 7) and annotates the designed
polyepitope protein, residue by residue, with a secondary-structure
propensity label computed from the published Chou-Fasman (1978)
conformational parameters, pinned as package data with their SHA-256
recorded in the artifact.

What this stage is NOT:
- not a 3-D model: no coordinates, no contact map, no docking, no loops;
- not a measurement and not a machine-learning prediction: each label is a
  fixed, documented function (windowed means + argmax) of the published
  parameters, so it is "derived" in the pipeline's evidence-layer sense;
- not validated against any experimental structure: that comparison is
  validation-pending and is declared as such.

Mechanical refusal (never a degraded annotation): missing construct artifact
(the chain must reach the design layer) or a protein with a residue absent
from the pinned parameter table.
"""
from __future__ import annotations

import csv
import hashlib
import time
from pathlib import Path

from neoforge.artifacts import ArtifactStore, StageRun
from neoforge.errors import ArtifactChainError

PARAMS_NAME = "chou_fasman_1978.csv"
CITATION = ("Chou & Fasman (1978), Advances in Enzymology 47:45-148; "
            "conformational parameters Pa (alpha-helix) and Pb (beta-sheet) "
            "(values reproduced in the ×100 prowl.rockefeller.edu/aainfo/chou.htm table)")
WINDOW = 6
COIL_P = 1.0
LABELS = ("H", "E", "C")
STAMP = ("sequence-derived propensity computed with pinned published parameters; "
         "not a measured or model-predicted 3-D structure")
NON_CLAIMS = [
    "not a 3-D model: no coordinates, no contact map, no loop modelling, no docking",
    "propensity != structure: a labelled window may fold differently in context",
    ("parameters were derived from 1970s globular-protein statistics, not "
     "from polyepitope vaccine strings; the string is out of the table's domain"),
    "coil is the fixed reference propensity 1.0, the same constant for every residue",
    "no comparison to an experimental structure has been performed (validation pending)",
]
WINDOW_DETAIL = f"truncated centred {WINDOW}-residue window (positions r-{WINDOW // 2}..r+{WINDOW // 2 - 1} inclusive)"


def load_parameters(path: str | Path | None = None) -> tuple[dict[str, tuple[float, float]], str]:
    file = Path(path) if path else Path(__file__).resolve().parents[1] / "data" / PARAMS_NAME
    digest = hashlib.sha256(file.read_bytes()).hexdigest()
    with file.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    params = {r["aa"].strip(): (float(r["pa"]), float(r["pb"])) for r in rows}
    if len(params) != 20 or any(k not in params for k in "ACDEFGHIKLMNPQRSTVWY"):
        raise ArtifactChainError(f"{file}: parameter table is not the 20-amino-acid set")
    return params, digest


def propensity_labels(protein: str, params: dict[str, tuple[float, float]],
                      window: int = WINDOW) -> tuple[str, list[dict]]:
    unknown = sorted(set(protein) - set(params))
    if unknown:
        raise ArtifactChainError(
            "structure refuses: residue(s) "
            f"{', '.join(unknown)} absent from the pinned parameter table")
    n = len(protein)
    half = window // 2
    labels: list[str] = []
    table: list[dict] = []
    for idx, aa in enumerate(protein):
        lo, hi = max(0, idx - (window - 1 - half)), min(n, idx + half + 1)
        seg = protein[lo:hi]
        pa = sum(params[a][0] for a in seg) / len(seg)
        pb = sum(params[a][1] for a in seg) / len(seg)
        best, label = COIL_P, "C"
        if pa > best:
            best, label = pa, "H"
        if pb > best:
            best, label = pb, "E"
        labels.append(label)
        table.append({"position": idx + 1, "aa": aa,
                      "pa_mean": round(pa, 4), "pb_mean": round(pb, 4),
                      "label": label})
    return "".join(labels), table


def _runs(labels: str, label: str, min_run: int = 3) -> list[dict]:
    out = []
    start = None
    for i, lab in enumerate(labels + "-"):
        if lab == label and start is None:
            start = i
        elif lab != label and start is not None:
            if i - start >= min_run:
                out.append({"start": start + 1, "end": i, "length": i - start})
            start = None
    return out


def run_structure(store: ArtifactStore, case: str, *, params_path: str | Path | None = None) -> dict:
    started = time.time()
    manifest = store.read_manifest(case)
    if "construct" not in manifest["stages"]:
        raise ArtifactChainError(
            f"structure refuses: no construct artifact for case '{case}' "
            "(the chain must reach the construct design layer)")

    con = store.load(case, "construct")
    protein = con.get("protein")
    if not protein:
        raise ArtifactChainError("structure refuses: construct artifact carries no designed protein")

    params, params_sha = load_parameters(params_path)
    labels, table = propensity_labels(protein, params)
    counts = {lab: labels.count(lab) for lab in LABELS}

    payload = {
        "stage": "structure",
        "case": case,
        "layer": "derived (fixed documented function of published parameters)",
        "evidence_kind": "computed_propensity",
        "stamp": STAMP,
        "non_claims": list(NON_CLAIMS),
        "inputs": {"construct": {"artifact_id": manifest["stages"]["construct"]["artifact_id"],
                                 "sha256": manifest["stages"]["construct"]["sha256"]}},
        "parameters": {
            "method": "per-residue argmax of windowed mean propensity over {H, E, C}, "
                      "tie-broken in the fixed order H, E, C; coil propensity is the "
                      f"constant {COIL_P}",
            "window_detail": WINDOW_DETAIL,
            "window": WINDOW,
            "chou_fasman": {"citation": CITATION, "file": PARAMS_NAME, "sha256": params_sha},
        },
        "protein": protein,
        "protein_length_aa": len(protein),
        "validation_status": "PENDING",
        "validation_note": ("labels have not been compared to any experimental structure; "
                            "the comparison is planned, not performed"),
        "propensity_table": {"columns": ["position", "aa", "pa_mean", "pb_mean", "label"],
                             **{k: [t[k] for t in table] for k in
                                ("position", "aa", "pa_mean", "pb_mean", "label")}},
        "summary": {"label_counts": counts,
                    "helix_runs_ge3": _runs(labels, "H"),
                    "sheet_runs_ge3": _runs(labels, "E")},
        "uncertainty_register": [
            "propensity labels are sequence statistics, not structure measurements",
            "the designed protein exists only in silico; no protein has been expressed or imaged",
            ("linker regions inherit the same propensity treatment, though their real "
             "conformation depends on proteolytic processing, which is not simulated"),
        ],
    }

    run = StageRun(stage="structure", case=case,
                   inputs={"construct": manifest["stages"]["construct"]["artifact_id"]},
                   parameters={"window": WINDOW, "params_sha256": params_sha})
    artifact_id, digest = store.write(run, "structure", payload)
    return {"artifact_id": artifact_id, "sha256": digest, "payload": payload,
            "wall": round(time.time() - started, 3)}
