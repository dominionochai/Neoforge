"""Real-case (HCC1395/SEQC2) validation exercises.

These are the EVALUATION_PLAN's L1/L2 checks on genuine public data. They are
marked slow, run from the repo work dir / reference tree, and skipped when
those are absent (paths overridable via NEOFORGE_REFS / NEOFORGE_WORK):

- the SEQC2 SNV truth set must re-verify 100% position containment in the v1.2
  high-confidence regions (its own release claim);
- every derived peptide window must equal the slice of the GENCODE v46 official
  protein translation — an independent check of the whole translation path.
"""
from __future__ import annotations

import gzip
import json
import os
from pathlib import Path

import pytest

from neoforge.artifacts import ArtifactStore
from neoforge.stages.ingest import run_ingest

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
WORK = Path(os.environ.get("NEOFORGE_WORK", REPO / "work"))


def load_work(store: ArtifactStore, stage: str) -> dict:
    record = store.read_manifest("hcc1395")["stages"][stage]
    return json.loads((store.case_dir("hcc1395") / record["path"]).read_text())


def test_real_ingest_snv_containment(tmp_path: Path) -> None:
    store = ArtifactStore(tmp_path / "work")
    run = run_ingest(BUNDLE, REFS, store)
    checks = {c["id"]: c for c in run["payload"]["checks"]}
    assert checks["somatic_snv_vcf.hc_region_containment"]["status"] == "pass"
    counters = run["payload"]["vcf"]["somatic_snv_vcf"]["counters"]
    assert counters["records"] == 39447
    assert counters["outside_hc_region_positions"] == 0


def test_real_translation_matches_gencode_proteins() -> None:
    store = ArtifactStore(WORK)
    if "translate" not in store.read_manifest("hcc1395")["stages"]:
        pytest.skip("translate artifact not built in work dir")
    derived = load_work(store, "translate")["derived"]
    cols = derived["columns"]

    proteins: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    fasta = REFS / "gencode/gencode.v46.pc_translations.fa.gz"
    with gzip.open(fasta, "rt") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    proteins[name] = "".join(chunks)
                name, chunks = line.split("|")[1], []
            else:
                chunks.append(line)
        if name is not None:
            proteins[name] = "".join(chunks)

    checked = 0
    unknown: set[str] = set()
    for row in zip(*[derived[c] for c in cols]):
        r = dict(zip(cols, row))
        protein = proteins.get(r["transcript_id"])
        if protein is None:
            unknown.add(r["transcript_id"])
            continue
        start, k = r["window_start"], r["peptide_length"]
        expected = protein[start - 1:start - 1 + k]
        mutant = expected[:r["mutant_offset"] - 1] + r["alt_aa"] + expected[r["mutant_offset"]:]
        assert r["ref_window_seq"] == expected, f"ref window mismatch {r['variant_key']} {r['transcript_id']}"
        assert r["window_seq"] == mutant, f"mutant window mismatch {r['variant_key']} {r['transcript_id']}"
        checked += 1
    assert not unknown
    assert checked > 0
