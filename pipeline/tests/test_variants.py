"""Ledger tests for the variants stage: conservation of records, honest reasons."""
from __future__ import annotations

from pathlib import Path

import pytest

from neoforge.artifacts import ArtifactStore
from neoforge.errors import ArtifactChainError
from neoforge.stages.ingest import run_ingest
from neoforge.stages.variants import run_variants
from tests.test_ingest import make_case


def seeded_store(tmp_path: Path, vcf_body: str) -> tuple[ArtifactStore, str]:
    bundle = make_case(tmp_path, vcf_body, "chr1\t0\t1000\n", "TUM-A,HLA-A*29:02;HLA-B*08:01")
    store = ArtifactStore(tmp_path / "work")
    run_ingest(bundle, tmp_path / "refs", store)
    return store, "mini"


VCF_MIXED = "\n".join([
    "chr1\t100\t.\tA\tG\t.\tPASS\tSOMATIC;TVAF=0.55;NVAF=0.01;bwaDP=40,60",
    "chr1\t200\t.\tC\tT,N\t.\tPASS\tSOMATIC",
    "chr1\t300\t.\tG\tC\t.\tPASS\tSOMATIC",
])


def test_selection_conservation(tmp_path: Path) -> None:
    store, case = seeded_store(tmp_path, VCF_MIXED)
    run = run_variants(store, case)
    ledger = run["payload"]["selection"]
    assert ledger["ingested_snv_records"] == 3
    assert ledger["selected"] == 2
    assert ledger["excluded"] == 1
    assert ledger["exclusion_reasons"]["multi_allelic_unresolvable_v0.1"] == 1
    assert run["payload"]["selected"]["variant_key"] == ["chr1:100A>G", "chr1:300G>C"]
    kept = run["payload"]["selected"]
    assert kept["tvaf"] == [0.55, None]
    assert store.load("mini", "variants")["selection"] == ledger


def test_variants_requires_ingest(tmp_path: Path) -> None:
    store = ArtifactStore(tmp_path / "work")
    with pytest.raises(ArtifactChainError, match="ingest"):
        run_variants(store, "mini")
