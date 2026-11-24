"""Immutable artifact store with SHA-256 content hashing and provenance.
Both the pipeline stages and the web app bind to these artifacts."""
from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _canonical_json(obj: Any) -> bytes:
    return json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()


@dataclass
class StageRun:
    stage: str
    case: str
    inputs: dict[str, str] = field(default_factory=dict)
    parameters: dict[str, Any] = field(default_factory=dict)
    started: float = field(default_factory=time.time)
    software: dict[str, str] = field(default_factory=dict)

    def record_software(self, texts: dict[str, str]) -> None:
        self.software.update(texts)

    def environment(self) -> dict[str, str]:
        return {
            "python_version": sys.version.split()[0],
            "platform": platform.platform(),
        }


class ArtifactStore:
    def __init__(self, work_dir: Path) -> None:
        self.work_dir = work_dir
        self.work_dir.mkdir(parents=True, exist_ok=True)

    def case_dir(self, case: str) -> Path:
        d = self.work_dir / case
        d.mkdir(parents=True, exist_ok=True)
        return d

    def read_manifest(self, case: str) -> dict:
        p = self.case_dir(case) / "manifest.json"
        if not p.exists():
            return {"stages": {}}
        return json.loads(p.read_text())

    def write(self, run: StageRun, name: str, payload: dict) -> tuple[str, str]:
        body = _canonical_json(payload)
        artifact_id = sha256_bytes(body)[:16]
        case = self.case_dir(run.case)
        path = case / f"{name}.{artifact_id}.json"
        path.write_bytes(body)
        digest = sha256_bytes(body)
        wall = round(time.time() - run.started, 3)
        record = {
            "stage": run.stage,
            "artifact_id": artifact_id,
            "path": path.name,
            "sha256": digest,
            "inputs": run.inputs,
            "parameters": run.parameters,
            "software": {**run.software, "python": sys.version.split()[0], "platform": platform.platform()},
            "wall_seconds": wall,
        }
        manifest = self.read_manifest(run.case)
        manifest.setdefault("stages", {})[name] = record
        (case / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
        return artifact_id, digest

    def load(self, case: str, name: str) -> dict:
        case_dir = self.case_dir(case)
        record = self.read_manifest(case)["stages"].get(name)
        if record is None:
            raise FileNotFoundError(f"artifact '{name}' not produced for case '{case}'")
        path = case_dir / record["path"]
        return json.loads(path.read_text())
