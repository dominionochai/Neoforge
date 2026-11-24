"""End-to-end six-stage CLI determinism check on the real case.

Runs the CLI for every stage (with a dev `--limit` on present so the suite
stays minutes, not hours) and then reruns the downstream stages: every
artifact id must be unchanged, which pins the SCIENTIFIC_DESIGN §5 claim that
a rerun with identical inputs reproduces hash-identical outputs."""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from neoforge.artifacts import ArtifactStore
from neoforge.cli import main

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(
        "NEOFORGE_REFS" not in os.environ
        and not Path(os.environ.get("NEOFORGE_REFS", "/workspace/refs"), "seqc2").exists(),
        reason="real reference data not present",
    ),
]

REFS = Path(os.environ.get("NEOFORGE_REFS", "/workspace/refs"))
REPO = Path(os.environ.get("NEOFORGE_REPO", Path(__file__).resolve().parents[3]))
BUNDLE = REPO / "cases/hcc1395/case.json"
CASES_DIR = REPO / "cases"
LIMIT = 2000


def run_all(work: Path) -> None:
    for stage in ("ingest", "variants", "translate", "present", "prioritize", "report"):
        argv = [stage, "--case", "hcc1395", "--work-dir", str(work),
                "--cases-dir", str(CASES_DIR), "--refs-dir", str(REFS)]
        if stage == "present":
            argv += ["--limit", str(LIMIT)]
        assert main(argv) == 0, f"CLI {stage} failed"


def test_six_stage_end_to_end_is_deterministic(tmp_path: Path) -> None:
    first, second = tmp_path / "run1", tmp_path / "run2"
    run_all(first)
    run_all(second)

    first_store, second_store = ArtifactStore(first), ArtifactStore(second)
    m1 = first_store.read_manifest("hcc1395")["stages"]
    m2 = second_store.read_manifest("hcc1395")["stages"]
    for stage in ("ingest", "variants", "translate", "present", "prioritize", "report"):
        assert m1[stage]["artifact_id"] == m2[stage]["artifact_id"], \
            f"{stage} did not reproduce a hash-identical artifact"

    report = second_store.load("hcc1395", "report")
    assert report["molecular_profile"]["peptide_windows"] >= LIMIT
    assert report["evidence_ledger"]["L4_experimental"]["status"] == "PENDING"
    chain = report["provenance_chain"]
    assert all(chain[s]["artifact_id"] == m1[s]["artifact_id"]
               for s in ("ingest", "variants", "translate", "present", "prioritize"))
