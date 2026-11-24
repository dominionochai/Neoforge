"""Case bundle: a self-describing JSON document for one research case.

A bundle pins every input file path (relative to a refs root), its SHA-256 and
its scientific role, plus the tumour/normal pairing claim that every later
artifact inherits. Loading it validates structure and hash well-formedness;
byte-level integrity is verified at ingest time, not here."""
from __future__ import annotations

import json
import re
from pathlib import Path

from neoforge.errors import CaseBundleError

_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")

REQUIRED_TOP = ("case_id", "display_name", "genome_build", "input_files", "hla", "pairing")
REQUIRED_INPUT_ENTRY = ("path", "sha256", "role")


def load_case_bundle(path: str | Path) -> dict:
    path = Path(path)
    if not path.exists():
        raise CaseBundleError(f"case bundle not found: {path}")
    bundle = json.loads(path.read_text())
    if not bundle.get("case_id"):
        raise CaseBundleError("case_id missing")
    for key in REQUIRED_TOP:
        if key not in bundle:
            raise CaseBundleError(f"{path.name}: required key '{key}' missing")
    if not isinstance(bundle["input_files"], dict) or not bundle["input_files"]:
        raise CaseBundleError("input_files must be a non-empty object")
    for name, entry in bundle["input_files"].items():
        if not isinstance(entry, dict):
            raise CaseBundleError(f"input_files.{name}: must be an object")
        for key in REQUIRED_INPUT_ENTRY:
            if key not in entry:
                raise CaseBundleError(f"input_files.{name}: required key '{key}' missing")
        sha = entry["sha256"]
        if not _SHA256_RE.match(sha):
            raise CaseBundleError(f"input_files.{name}: sha256 must be 64 lowercase hex chars")
        if "path" in entry and entry["path"].startswith("/"):
            raise CaseBundleError("input_files paths must be relative to the refs root")
    if bundle["genome_build"] != "GRCh38":
        raise CaseBundleError(f"unsupported genome_build {bundle['genome_build']!r} (this pipeline targets GRCh38)")
    pairing = bundle["pairing"]
    if not all(("cell_line" in pairing.get("tumour", {}), "cell_line" in pairing.get("normal", {}))):
        raise CaseBundleError("pairing must declare tumour and normal cell lines")
    return bundle


def resolve_input(refs_root: Path, entry: dict) -> Path:
    resolved = (refs_root / entry["path"]).resolve()
    if refs_root.resolve() not in resolved.parents and resolved.parent != refs_root.resolve():
        raise CaseBundleError(f"input path {entry['path']!r} escapes the refs root")
    return resolved
