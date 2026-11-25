"""End-to-end CLI determinism check on the real case.

Runs the CLI for every stage (with a dev `--limit` on present so the suite
stays minutes, not hours) and then reruns the downstream stages: every
artifact id must be unchanged, which pins the SCIENTIFIC_DESIGN §5 claim that
a rerun with identical inputs reproduces hash-identical outputs. The construct
stage is exercised last-gated: with the dev limit some windows may fail the
mutant-vs-wild-type discrimination check, so refusal is a legal outcome and is
asserted honestly either way."""
from __future__ import annotations

import os
from pathlib import Path

import pytest
from neoforge.artifacts import ArtifactStore
from neoforge.cli import main
from neoforge.errors import ArtifactChainError

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(
        "NEOFORGE_REFS" not in os.environ
        and not Path(os.environ.get("NEOFORGE_REFS", "/workspace/refs"), "seqc2").exists(),
        reason="real reference data not present",
    ),
]

REFS = Path(os.environ.get("NEOFORGE_REFS", "/workspace/refs"))
REPO = Path(os.environ.get("NEOFORGE_REPO", Path(__file__).resolve().parents[2]))
BUNDLE = REPO / "cases/hcc1395/case.json"
CASES_DIR = REPO / "cases"
LIMIT = 2000


def run_all(work: Path) -> bool:
    construct_ran = True
    for stage in ("ingest", "variants", "translate", "present", "prioritize",
                  "construct", "report"):
        argv = [stage, "--case", "hcc1395", "--work-dir", str(work),
                "--cases-dir", str(CASES_DIR), "--refs-dir", str(REFS)]
        if stage == "present":
            argv += ["--limit", str(LIMIT)]
        try:
            rc = main(argv)
        except ArtifactChainError:
            if stage == "construct":
                construct_ran = False
                continue
            raise
        assert rc == 0, f"CLI {stage} failed"
    return construct_ran


def test_pipeline_end_to_end_is_deterministic(tmp_path: Path) -> None:
    first, second = tmp_path / "run1", tmp_path / "run2"
    ran1, ran2 = run_all(first), run_all(second)
    assert ran1 == ran2, "construct refusal is nondeterministic between runs"

    first_store, second_store = ArtifactStore(first), ArtifactStore(second)
    m1 = first_store.read_manifest("hcc1395")["stages"]
    m2 = second_store.read_manifest("hcc1395")["stages"]
    stages = ("ingest", "variants", "translate", "present", "prioritize", "report")
    if ran1:
        stages = stages[:-1] + ("construct", "report")
    elif "construct" in m1 or "construct" in m2:
        pytest.fail("construct run flag says refused but a construct artifact exists")
    for stage in stages:
        assert m1[stage]["artifact_id"] == m2[stage]["artifact_id"], \
            f"{stage} did not reproduce a hash-identical artifact"

    report = second_store.load("hcc1395", "report")
    assert report["molecular_profile"]["peptide_windows"] >= LIMIT
    assert report["evidence_ledger"]["L4_experimental"]["status"] == "PENDING"
    chain = report["provenance_chain"]
    chain_stages = ("ingest", "variants", "translate", "present", "prioritize")
    if ran1:
        chain_stages += ("construct",)
        assert ran1 == ("construct" in chain), \
            "provenance chain disagrees with the construct outcome"
    else:
        assert "construct" not in chain, "refused construct must not appear in the chain"
    assert all(chain[s]["artifact_id"] == m1[s]["artifact_id"]
               for s in chain_stages)
