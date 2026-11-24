"""NEOFORGE pipeline command line: one stage per invocation, artifacts in a work dir."""
from __future__ import annotations

import argparse
from pathlib import Path

from neoforge import __version__
from neoforge.artifacts import ArtifactStore

STAGES = ("ingest", "variants", "translate", "present", "prioritize", "report")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="neoforge")
    ap.add_argument("--version", action="version", version=__version__)
    ap.add_argument("stage", choices=STAGES)
    ap.add_argument("--case", required=True)
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument("--cases-dir", type=Path, default=Path("cases"))
    ap.add_argument("--refs-dir", type=Path, required=True)
    ap.add_argument("--limit", type=int, help="dev: cap windows for the present stage")
    args = ap.parse_args(argv)

    from neoforge.stages.ingest import run_ingest  # noqa: PLC0415 — stages import lazily

    bundle = args.cases_dir / args.case / "case.json"
    store = ArtifactStore(args.work_dir)
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

        return _emit_present(run_present(store, args.case, limit=args.limit))
    raise NotImplementedError(f"stage '{args.stage}' connects when implemented")


def _emit_present(result: dict) -> int:
    p = result["payload"]
    print(f"present: artifact {result['artifact_id']} [{p['layer']}, evidence={p['evidence_kind']}]")
    print(f"  candidates scored: {p['candidates_scored']}; predictions made: {p['predictions_made']} "
          f"in {p['wall_seconds']}s")
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
