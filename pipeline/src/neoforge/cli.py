from __future__ import annotations

import argparse
from pathlib import Path

from neoforge import __version__


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="neoforge")
    ap.add_argument("--version", action="version", version=__version__)
    ap.add_argument("stage", choices=["ingest", "variants", "translate", "present", "prioritize", "report"])
    ap.add_argument("--case", required=True)
    ap.add_argument("--work-dir", type=Path, required=True)
    args = ap.parse_args(argv)
    print(f"stage '{args.stage}' for case '{args.case}' (placeholder wiring lands with pipeline run)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
