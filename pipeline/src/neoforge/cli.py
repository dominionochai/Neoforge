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
    args = ap.parse_args(argv)

    from neoforge.stages.ingest import run_ingest  # noqa: PLC0415 — stages import lazily

    bundle = args.cases_dir / args.case / "case.json"
    store = ArtifactStore(args.work_dir)
    result = run_ingest(bundle, args.refs_dir, store)
    checks = {c["id"]: c["status"] for c in result["payload"]["checks"]}
    failed = [k for k, v in checks.items() if v != "pass"]
    print(f"ingest {args.case}: artifact {result['artifact_id']} sha256={result['sha256'][:16]}…")
    for check_id, status in checks.items():
        print(f"  [{status}] {check_id}")
    if failed:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
