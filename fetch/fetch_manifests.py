"""Download reference artifacts and verify pinned SHA-256 checksums."""
from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path

import requests


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(rows: list[dict], out_dir: Path, allow_unpinned: bool) -> int:
    failures = 0
    for row in rows:
        name, url, sha = row["name"], row["url"].split()[0], row["sha256"]
        if url.startswith("(") or not url.startswith("http"):
            print(f"derived {name} (built locally by a fetch/ script; not downloadable)")
            continue
        target = out_dir / name
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and sha != "pinned_at_fetch" and sha256_file(target) == sha:
            print(f"ok    {name}")
            continue
        print(f"fetch {name} ...", flush=True)
        with requests.get(url, stream=True, timeout=600) as resp:
            resp.raise_for_status()
            with target.open("wb") as fh:
                for chunk in resp.iter_content(1 << 21):
                    fh.write(chunk)
        digest = sha256_file(target)
        if sha == "pinned_at_fetch":
            print(f"UNPINNED {name} sha256={digest}")
            if not allow_unpinned:
                failures += 1
        elif digest != sha:
            print(f"FAIL  {name}: {digest} != {sha}")
            failures += 1
        else:
            print(f"ok    {name}")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, default=Path(__file__).parent / "manifests.csv")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--allow-unpinned", action="store_true")
    args = ap.parse_args()
    with args.manifest.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    return 1 if fetch(rows, args.out_dir, args.allow_unpinned) else 0


if __name__ == "__main__":
    sys.exit(main())
