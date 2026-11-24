"""UI projection — the bridge between the artifact store and the frontend.

The frontend renders ONLY what this module exports; it never computes science
of its own. Every value here is copied or trivially re-shaped from verified
artifacts; per-block `source` strings name the artifact (id + SHA-256) each
section was read from, so the interface can show the exact upstream digest.

Not a pipeline stage: it writes nothing into the artifact store and adds no
scientific claims. Run: `python -m neoforge.export_ui --work-dir ../work \
--cases-dir ../cases --out ../frontend/data`.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from neoforge.artifacts import ArtifactStore
from neoforge.errors import ArtifactChainError


def _load(store: ArtifactStore, case: str, stage: str) -> tuple[dict, dict]:
    """Return (payload, manifest record) for one stage, or (None, None)."""
    record = store.read_manifest(case)["stages"].get(stage)
    if record is None:
        return None, None
    payload = store.load(case, stage)
    record = {**record, "stage": stage, "case": case}
    return payload, record


def _project_ingest(ing: dict, record: dict) -> dict:
    tables = ing["variant_tables"]
    snv = tables.get("somatic_snv_vcf") or tables.get("somatic_maf") or {"pos": []}
    indel_table = tables.get("somatic_indel_vcf") or tables.get("somatic_indel_maf")
    indel_kinds = [k for k in (indel_table or {}).get("filters", [])]
    return {
        "stage": "ingest",
        "source": f"{record['stage']}.{record['artifact_id']}.{record['sha256']}",
        "case": ing["case"],
        "checks": ing["checks"],
        "vaf_source": (ing["vcf"].get("somatic_snv_vcf") or {}).get("counters") if ing.get("vcf") else None,
        "maf": ing.get("maf"),
        "inputs": ing["inputs"],
        "variant_counts": {
            "somatic_snv_positions": len(snv["pos"]),
            "somatic_indel_positions": len((indel_table or {}).get("pos", []) or []),
            "indel_kind_counts": {k: indel_kinds.count(k) for k in sorted(set(indel_kinds))},
        },
        "snv_table": {"columns": ["chrom", "pos", "ref", "alt", "tvaf", "nvaf", "tumor_depth"],
                      **{c: snv[c] for c in ("chrom", "pos", "ref", "alt", "tvaf", "nvaf",
                                             "tumor_depth") if c in snv}},
        "hla": ing["hla"],
        "expression": ({k: v for k, v in ing["expression"].items()
                        if k != "tpm_by_gene"} | {"n_genes": len(ing["expression"].get("tpm_by_gene", {}))})
                      if ing.get("expression") else None,
    }


def _shortlist_keys(report: dict | None, prioritize: dict | None) -> list[str]:
    if report and report.get("shortlist"):
        return [r["variant_key"] for r in report["shortlist"]]
    if prioritize:
        return prioritize["candidates"]["variant_key"]
    return []


def _project_translate(tra: dict | None, record: dict, keep_keys: set[str] | None) -> dict | None:
    if tra is None:
        return None
    derived = tra.get("derived", {})
    keys = derived.get("variant_key", [])
    idx = [i for i, k in enumerate(keys) if keep_keys is None or k in keep_keys]
    return {
        "stage": "translate",
        "source": f"{record['stage']}.{record['artifact_id']}.{record['sha256']}"
                  if record else "translate artifact",
        "summary": tra["summary"],
        "derived": {"columns": derived["columns"],
                    **{c: [derived[c][i] for i in idx] for c in derived["columns"]}},
        "derived_rows_total": len(keys),
        "derived_rows_exported": len(idx),
        "export_scope": "all windows" if keep_keys is None else "windows of shortlisted/prioritised variants only",
    }


def export_case(store: ArtifactStore, case: str, cases_dir: Path) -> dict:
    ing, ing_rec = _load(store, case, "ingest")
    if ing is None:
        raise ArtifactChainError(f"no ingest artifact for case '{case}'; run the pipeline first")
    var, var_rec = _load(store, case, "variants")
    tra, tra_rec = _load(store, case, "translate")
    pre, pre_rec = _load(store, case, "present")
    pri, pri_rec = _load(store, case, "prioritize")
    rep, rep_rec = _load(store, case, "report")

    bundle = json.loads((cases_dir / case / "case.json").read_text())
    shortlist_keys = set(_shortlist_keys(rep, pri))

    # predictions: only what shortlisted/prioritised variants need on screen —
    # the full 200k-row prediction table stays in the artifact store (inspectable
    # hashes listed in the workbench view); the frontend workload is bounded.
    predictions_export = None
    if pre is not None:
        preds = pre["predictions"]
        idx = [i for i, k in enumerate(preds["variant_key"]) if not shortlist_keys or k in shortlist_keys]
        predictions_export = {
            "columns": preds["columns"],
            **{c: [preds[c][i] for i in idx] for c in preds["columns"]},
            "rows_total": len(preds["variant_key"]),
            "rows_exported": len(idx),
            "export_scope": ("all prediction rows" if not shortlist_keys
                             else "prediction rows of the report shortlist variants"),
        }

    payload = {
        "exported_by": "neoforge.export_ui (no scientific computation; projections of hashed artifacts)",
        "pipeline_stage_records": {s: r for s, r in (("ingest", ing_rec), ("variants", var_rec),
                                                     ("translate", tra_rec), ("present", pre_rec),
                                                     ("prioritize", pri_rec), ("report", rep_rec)) if r},
        "case_bundle": bundle,
        "report": rep,
        "prioritize": ({"stage": "prioritize",
                        "source": f"{pri_rec['stage']}.{pri_rec['artifact_id']}.{pri_rec['sha256']}",
                        "tier_counts": pri["tier_counts"],
                        "parameters": pri["parameters"],
                        "candidates": pri["candidates"],
                        "candidates_prioritised": pri.get("candidates_prioritised"),
                        "non_claims": pri["non_claims"]} if pri else None),
        "present": ({"stage": "present",
                     "source": f"{pre_rec['stage']}.{pre_rec['artifact_id']}.{pre_rec['sha256']}",
                     "parameters": pre["parameters"],
                     "candidates_scored": pre["candidates_scored"],
                     "predictions_made": pre["predictions_made"],
                     "best_affinity_by_candidate": pre.get("best_affinity_by_candidate"),
                     "predictions": predictions_export,
                     "non_claims": pre["non_claims"]} if pre else None),
        "variants": ({"stage": "variants",
                      "source": f"{var_rec['stage']}.{var_rec['artifact_id']}.{var_rec['sha256']}",
                      "selection": var["selection"],
                      "selected": var["selected"]} if var else None),
        "translate": None,
        "ingest": None,
    }
    payload["translate"] = _project_translate(tra, tra_rec, shortlist_keys or None)
    payload["ingest"] = _project_ingest(ing, ing_rec)
    return payload


def export_tesla(store: ArtifactStore) -> dict | None:
    record = store.read_manifest("tesla")["stages"].get("tesla_eval")
    if record is None:
        return None
    payload = store.load("tesla", "tesla_eval")
    rows = payload["predictions"]
    n = len(rows["peptide"])
    validated_rows = [rows["pmhc_id"][i] for i in range(n) if rows["validated"][i]]
    return {
        "stage": "tesla_eval",
        "source": f"tesla.{record['stage']}.{record['artifact_id']}.{record['sha256']}",
        "parameters": payload["parameters"],
        "evidence_kind": payload["evidence_kind"],
        "layer": payload["layer"],
        "non_claims": payload["non_claims"],
        "metrics": payload["metrics"],
        "rank_rule_falsified": payload["rank_rule_falsified"],
        "validated_predictions": {
            "columns": rows["columns"],
            "columns_dropped": [c for c in rows["columns"] if c not in
                                ("pmhc_id", "peptide", "allele", "mhcflurry_affinity_nm",
                                 "rank", "validated", "measured_binding_affinity_nm")],
            **{c: [rows[c][i] for i in range(n) if rows["validated"][i]] for c in rows["columns"]},
        },
        "validated_pmhc_ids": validated_rows,
    }


def main() -> int:
    ap = argparse.ArgumentParser(prog="neoforge.export_ui")
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument("--cases-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--case", action="append", default=[],
                    help="case to export (repeatable); default: every case with an ingest artifact")
    args = ap.parse_args()

    store = ArtifactStore(args.work_dir)
    args.out.mkdir(parents=True, exist_ok=True)
    cases = args.case or sorted(p.name for p in args.cases_dir.iterdir() if p.is_dir())
    for case in cases:
        try:
            payload = export_case(store, case, args.cases_dir)
        except ArtifactChainError as exc:
            print(f"export_ui: skipped case '{case}' ({exc})")
            continue
        path = args.out / f"case_{case}.json"
        path.write_text(json.dumps(payload, sort_keys=True, separators=(",", ":")))
        print(f"export_ui: {path} ({path.stat().st_size} bytes)")
    tesla = export_tesla(store)
    if tesla:
        path = args.out / "benchmark_tesla.json"
        path.write_text(json.dumps(tesla, sort_keys=True,
                                   separators=(",", ":")))
        print(f"export_ui: {path} ({path.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
