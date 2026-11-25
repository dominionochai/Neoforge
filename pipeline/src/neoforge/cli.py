"""NEOFORGE pipeline command line: one stage per invocation, artifacts in a work dir."""
from __future__ import annotations

import argparse
from pathlib import Path

from neoforge import __version__
from neoforge.artifacts import ArtifactStore

STAGES = ("ingest", "variants", "translate", "present", "prioritize", "report", "tesla_eval")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="neoforge")
    ap.add_argument("--version", action="version", version=__version__)
    ap.add_argument("stage", choices=STAGES)
    ap.add_argument("--case", help="case name; all stages except tesla_eval")
    ap.add_argument("--tesla-xlsx", type=Path, help="tesla_eval: TESLA labels workbook (mmc4.xlsx)")
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument("--cases-dir", type=Path, default=Path("cases"))
    ap.add_argument("--refs-dir", type=Path, required=True)
    ap.add_argument("--limit", type=int, help="dev: cap windows for the present stage")
    args = ap.parse_args(argv)
    if args.stage == "tesla_eval":
        if not args.tesla_xlsx:
            ap.error("tesla_eval requires --tesla-xlsx (TESLA mmc4.xlsx labels workbook)")
    elif not args.case:
        ap.error(f"--case is required for stage '{args.stage}'")

    from neoforge.stages.ingest import run_ingest  # noqa: PLC0415 — stages import lazily

    store = ArtifactStore(args.work_dir)
    if args.stage == "tesla_eval":
        from neoforge.evaluation.tesla import run_tesla_evaluation

        return _emit_tesla(run_tesla_evaluation(store, args.tesla_xlsx))
    bundle = args.cases_dir / args.case / "case.json"
    if args.stage == "ingest":
        return _emit_ingest(run_ingest(bundle, args.refs_dir, store))
    if args.stage == "variants":
        from neoforge.stages.variants import run_variants

        return _emit_selection(run_variants(store, args.case))
    if args.stage == "translate":
        from neoforge.sequence import IndexedFasta
        from neoforge.stages.translate import run_translate

        gtf = args.refs_dir / "gencode/gencode.v46.annotation.gtf.gz"
        fasta_gz = args.refs_dir / "gencode/GRCh38.primary_assembly.genome.fa.gz"
        genome = IndexedFasta(IndexedFasta.decompress_if_needed(fasta_gz, args.refs_dir / "gencode"))
        result = run_translate(store, args.case, gtf, genome)
        print(f"translate: artifact {result['artifact_id']}")
        for key, value in result["payload"]["summary"].items():
            print(f"  {key}: {value}")
        return 0
    if args.stage == "present":
        from neoforge.stages.present import run_present

        result = run_present(store, args.case, limit=args.limit)
        record = store.read_manifest(args.case)["stages"]["present"]
        print(f"present: wall_seconds={record['wall_seconds']}")
        return _emit_present(result)
    if args.stage == "prioritize":
        from neoforge.stages.prioritize import run_prioritize

        return _emit_prioritize(run_prioritize(store, args.case))
    if args.stage == "report":
        from neoforge.stages.report import run_report

        return _emit_report(run_report(store, args.case))
    raise NotImplementedError(f"stage '{args.stage}' connects when implemented")


def _emit_tesla(result: dict) -> int:
    p = result["payload"]
    m = p["metrics"]
    verdict = "FALSIFIED" if p["rank_rule_falsified"] else "not falsified"
    print(f"tesla_eval: artifact {result['artifact_id']} [{p['layer']}, evidence={p['evidence_kind']}]")
    print(f"  candidates: {m['n_candidates']} ({m['n_validated']} validated) across "
          f"{m['n_patients']} patients / {m['n_alleles']} alleles")
    print(f"  ttif_top20pct: {m['ttif_top20pct']:.3f} (random {m['ttif_random_expectation']:.3f}, "
          f"oracle {m['ttif_oracle']:.3f})")
    print(f"  auprc: {m['auprc']:.3f} (random {m['auprc_random_expectation']:.3f}, "
          f"oracle {m['auprc_oracle']:.3f})")
    print(f"  median rank of validated: {m['median_rank_of_validated']}")
    for e in (p.get("rule_ablation") or {}).get("entries", []):
        print(f"  ablation [{e['criterion']}]: ttif {e['ttif_top20pct']:.3f}, auprc {e['auprc']:.3f}, "
              f"median rank {e['median_rank_of_validated']}")
    print(f"  falsification verdict ({p['parameters']['falsification_clause']}): {verdict}")
    return 0


def _emit_report(result: dict) -> int:
    p = result["payload"]
    print(f"report: artifact {result['artifact_id']} [{p['layer']}, evidence={p['evidence_kind']}]")
    print(f"  shortlist: {len(p['shortlist'])} of {sum(p['tier_counts'].values())} candidates")
    return 0


def _emit_prioritize(result: dict) -> int:
    p = result["payload"]
    print(f"prioritize: artifact {result['artifact_id']} [{p['layer']}, evidence={p['evidence_kind']}]")
    print(f"  candidates: {p['candidates_prioritised']}; tiers: {p['tier_counts']}")
    return 0


def _emit_present(result: dict) -> int:
    p = result["payload"]
    print(f"present: artifact {result['artifact_id']} [{p['layer']}, evidence={p['evidence_kind']}]")
    print(f"  candidates scored: {p['candidates_scored']}; predictions made: {p['predictions_made']}")
    return 0


def _emit(result: dict, lines: list[str]) -> int:
    print(f"{result['payload'].get('stage', '')} artifact {result['artifact_id']}")
    for line in lines:
        print(" ", line)
    return 0


def _emit_ingest(result: dict) -> int:
    checks = {c["id"]: c["status"] for c in result["payload"]["checks"]}
    failed = [k for k, v in checks.items() if v == "fail"]
    print(f"ingest: artifact {result['artifact_id']} sha256={result['sha256'][:16]}…")
    for check_id, status in checks.items():
        print(f"  [{status}] {check_id}")
    return 1 if failed else 0


def _emit_selection(result: dict) -> int:
    sel = result["payload"]["selection"]
    print(f"variants: artifact {result['artifact_id']}")
    print(f"  selected {sel['selected']} of {sel['ingested_snv_records']} SNVs; "
          f"exclusions: {sel['exclusion_reasons'] or 'none'}; indels seen: {sel['indels_seen']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
