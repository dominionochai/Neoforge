#!/usr/bin/env python3
"""NEOFORGE one-command demo: pipeline stages -> exports -> identity
verification -> local server. Cross-platform: runs identically under
PowerShell (`python scripts/run_demo.py`), bash and any other shell with
Python 3.11+ on PATH.

Every number shown in the interface is recomputed from the reference data in
REFS (deterministic, artefact-hashed); the export step then verifies that the
projections are identical to the tracked frontend/data/ files (modulo per-run
wall-clock timings, which legitimately vary), so what you browse is provably
what the pipeline produced.

Usage:
  python scripts/run_demo.py                 # run everything, serve on 127.0.0.1:8000
  python scripts/run_demo.py --no-serve      # run everything, skip the server
  python scripts/run_demo.py --fetch-refs    # fetch missing reference data first

Environment:
  NEOFORGE_REFS   reference-data directory (no portable default: derives from
                  NEOFORGE_REFS, else ./refs if it exists, else /workspace/refs;
                  build it with fetch/fetch_manifests.py or --fetch-refs)
  PORT            server port (default 8000)
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CASES = ("hcc1395", "skcm_d9a3z1")
HCC_STAGES = ("ingest", "variants", "translate", "present", "prioritize", "construct", "report")
SKCM_STAGES = ("ingest", "variants", "translate")  # honestly stops at L2: no HLA genotype
FRONTEND_DATA = ROOT / "frontend" / "data"


def die(msg: str) -> None:
    print(f"run_demo: {msg}", file=sys.stderr)
    raise SystemExit(1)


def resolve_refs(explicit: str | None) -> Path:
    if explicit:
        return Path(explicit)
    env = os.environ.get("NEOFORGE_REFS")
    if env:
        return Path(env)
    if (ROOT / "refs").is_dir():
        return ROOT / "refs"
    return Path("/workspace/refs")


def preflight_install() -> None:
    probe = subprocess.run([sys.executable, "-c", "import neoforge"],
                           capture_output=True)
    if probe.returncode != 0:
        print("run_demo: installing the pipeline (pip install -e ./pipeline)...")
        subprocess.run([sys.executable, "-m", "pip", "install", "-e", str(ROOT / "pipeline")],
                       check=True)


def preflight_mhcflurry() -> None:
    """present/construct need mhcflurry AND its downloaded models; fail with
    the exact remediation instead of an unexplained traceback mid-pipeline."""
    probe = (
        "from mhcflurry import Class1PresentationPredictor\n"
        "Class1PresentationPredictor.load().predict(['SYFPEITHI'], ['HLA-A*02:01'])\n"
    )
    got = subprocess.run([sys.executable, "-c", probe], capture_output=True, text=True)
    if got.returncode != 0:
        tail = (got.stderr or got.stdout).strip().splitlines()[-4:]
        print("\n".join(f"    {line}" for line in tail), file=sys.stderr)
        die("mhcflurry is unusable; install it and fetch its models:\n"
            f"    {sys.executable} -m pip install mhcflurry\n"
            "    mhcflurry-downloads fetch models_class1_presentation")


def run_stage(refs: Path, case: str, stage: str) -> None:
    print(f"==> [{case}] {stage}")
    cmd = [sys.executable, "-m", "neoforge.cli", stage, "--case", case,
           "--work-dir", str(ROOT / "work"), "--cases-dir", str(ROOT / "cases"),
           "--refs-dir", str(refs)]
    got = subprocess.run(cmd)
    if got.returncode != 0:
        die(f"stage '{stage}' failed for case '{case}' (exit {got.returncode})")


def load_projection(path: Path):
    if path.name.endswith(".fasta"):
        return path.read_bytes()
    def strip_volatile(obj):
        if isinstance(obj, dict):
            return {k: strip_volatile(v) for k, v in obj.items() if k != "wall_seconds"}
        if isinstance(obj, list):
            return [strip_volatile(v) for v in obj]
        return obj
    raw = gzip.open(path).read() if path.name.endswith(".json.gz") else path.read_bytes()
    return strip_volatile(json.loads(raw))


def export_and_verify() -> list[str]:
    print("==> export projections")
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "data"
        got = subprocess.run([sys.executable, "-m", "neoforge.export_ui",
                              "--work-dir", str(ROOT / "work"),
                              "--cases-dir", str(ROOT / "cases"), "--out", str(out)])
        if got.returncode != 0:
            die("export_ui failed")
        produced = sorted(out.iterdir())
        if not produced:
            die("export produced no files at all")
        print("==> verifying exports match frontend/data, modulo per-run wall-clock")
        for src in produced:
            tracked = FRONTEND_DATA / src.name
            if not tracked.exists():
                die(f"export produced {src.name} but frontend/data/{src.name} is not tracked; "
                    "regenerate with neoforge.export_ui before serving")
            if load_projection(src) != load_projection(tracked):
                die(f"{src.name} differs from the tracked export (beyond wall-clock timings) — "
                    "regenerate with neoforge.export_ui before serving")
            print(f"    ok: {src.name}")
    return [p.name for p in produced]


def serve(port: int) -> None:
    print(f"==> serving {ROOT / 'frontend'} at http://127.0.0.1:{port} (Ctrl-C to stop)")
    import functools
    from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler

    handler = functools.partial(SimpleHTTPRequestHandler, directory=str(ROOT / "frontend"))
    ThreadingHTTPServer(("127.0.0.1", port), handler).serve_forever()


def main() -> int:
    ap = argparse.ArgumentParser(prog="run_demo")
    ap.add_argument("--no-serve", action="store_true", help="skip the local server")
    ap.add_argument("--port", type=int, default=int(os.environ.get("PORT", "8000")))
    ap.add_argument("--refs-dir", default=None, help="overrides NEOFORGE_REFS")
    ap.add_argument("--fetch-refs", action="store_true",
                    help="fetch/refresh reference data (pinned SHA-256) before running")
    args = ap.parse_args()

    refs = resolve_refs(args.refs_dir)
    wanted = ["gencode/gencode.v46.annotation.gtf.gz",
              "gencode/GRCh38.primary_assembly.genome.fa.gz",
              "hla/cell-line-class1-mhc-alleles.csv",
              "mhci2018_mmc4.xlsx", "tesla/mmc4.xlsx"]
    missing = [w for w in wanted if not (refs / w).exists()]
    if missing and not args.fetch_refs:
        die("reference data not found in "
            f"'{refs}' (missing: {', '.join(missing)}).\n"
            "Set NEOFORGE_REFS, pass --refs-dir, or rebuild it:\n"
            f"    python fetch/fetch_manifests.py --out-dir <dir>\n"
            "    python scripts/run_demo.py --refs-dir <dir> --no-serve\n"
            "Note: the genome and GTF downloads total roughly 900 MB.")
    if missing:
        print(f"==> fetching reference data into {refs} ({len(missing)} file(s) missing)")
        got = subprocess.run([sys.executable, str(ROOT / "fetch" / "fetch_manifests.py"),
                              "--out-dir", str(refs)])
        if got.returncode != 0:
            die("reference fetch failed (checksums are pinned in fetch/manifests.csv)")

    preflight_install()
    preflight_mhcflurry()

    # construct runs before report so the research summary can embed the design;
    # skcm stops after translate: without HLA the later stages refuse mechanically.
    for stage in HCC_STAGES:
        run_stage(refs, "hcc1395", stage)
    for stage in SKCM_STAGES:
        run_stage(refs, "skcm_d9a3z1", stage)
    print("==> [tesla] tesla_eval")
    got = subprocess.run([sys.executable, "-m", "neoforge.cli", "tesla_eval",
                          "--tesla-xlsx", str(refs / "tesla" / "mmc4.xlsx"),
                          "--work-dir", str(ROOT / "work"),
                          "--cases-dir", str(ROOT / "cases"), "--refs-dir", str(refs)])
    if got.returncode != 0:
        die("tesla_eval failed")

    export_and_verify()
    if args.no_serve:
        print("run_demo: done (--no-serve)")
    else:
        serve(args.port)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
