"""Tests for the L1 ingest stage: parsers, integrity gates, evidence honesty."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from neoforge.artifacts import ArtifactStore
from neoforge.bed import IntervalIndex
from neoforge.case_config import load_case_bundle
from neoforge.errors import CaseBundleError, InputIntegrityError
from neoforge.hla import normalise_alleles
from neoforge.stages.ingest import run_ingest
from neoforge.vcf import VcfFormatError, VcfReader

VCF_HEADER = "\n".join([
    "##fileformat=VCFv4.1",
    "##reference=GRCh38.d1.vd1",
    "##INFO=<ID=TVAF,Number=1,Type=Float,Description=\"tumor VAF\">",
    "##INFO=<ID=NVAF,Number=1,Type=Float,Description=\"normal VAF\">",
    "##INFO=<ID=bwaDP,Number=2,Type=Integer,Description=\"vaf depth,total\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
])


def make_case(tmp_path: Path, vcf_body: str, bed_body: str, alleles: str) -> Path:
    refs = tmp_path / "refs"
    case = tmp_path / "cases" / "mini"
    (refs / "seqc2").mkdir(parents=True)
    (refs / "hla").mkdir(parents=True)
    case.mkdir(parents=True)
    (refs / "seqc2" / "mini.vcf").write_text(VCF_HEADER + "\n" + vcf_body)
    (refs / "seqc2" / "mini.bed").write_text(bed_body)
    (refs / "hla" / "alleles.csv").write_text(alleles)
    entry = {
        "somatic_snv_vcf": {
            "path": "seqc2/mini.vcf",
            "sha256": hashlib.sha256((refs / "seqc2" / "mini.vcf").read_bytes()).hexdigest(),
            "role": "test SNV set",
        },
        "hc_regions_bed": {
            "path": "seqc2/mini.bed",
            "sha256": hashlib.sha256((refs / "seqc2" / "mini.bed").read_bytes()).hexdigest(),
            "role": "HC regions",
        },
        "hla_csv": {
            "path": "hla/alleles.csv",
            "sha256": hashlib.sha256((refs / "hla" / "alleles.csv").read_bytes()).hexdigest(),
            "role": "HLA source",
        },
    }
    (case / "case.json").write_text(json.dumps({
        "case_id": "mini",
        "display_name": "Synthetic mini case",
        "genome_build": "GRCh38",
        "pairing": {
            "tumour": {"cell_line": "TUM-A", "type": "test tumour"},
            "normal": {"cell_line": "NRM-A", "type": "test normal"},
        },
        "input_files": entry,
        "hla": {"cell_line": "TUM-A"},
    }))
    return case / "case.json"


BED = "chr1\t0\t1000\n"
VCF_FULL = "\n".join([
    "chr1\t100\t.\tA\tG\t.\tPASS;HighConf\tSOMATIC;TVAF=0.55;NVAF=0.01;bwaDP=40,60",
    "chr1\t200\t.\tC\tT\t.\tPASS\tSOMATIC",
    "chr1\t300\t.\tG\tC,G\t.\tLowConf\tSOMATIC",
])
VCF_OUTSIDE = "chr2\t100\t.\tA\tT\t.\tPASS\tSOMATIC"
ALLELES = "TUM-A,HLA-A*29:02;HLA-B*08:01"


def read_table(run: dict, vcf_name: str) -> dict:
    return run["payload"]["variant_tables"][vcf_name]


def test_ingest_parses_and_reports_vaf_availability(tmp_path: Path) -> None:
    bundle = make_case(tmp_path, VCF_FULL, BED, ALLELES)
    store = ArtifactStore(tmp_path / "work")
    run = run_ingest(bundle, tmp_path / "refs", store)
    table = read_table(run, "somatic_snv_vcf")
    assert table["tvaf"] == [0.55, None, None]
    assert table["nvaf"] == [0.01, None, None]
    assert table["tumor_depth"] == [40, None, None]
    assert run["payload"]["vcf"]["somatic_snv_vcf"]["counters"]["records"] == 3
    assert run["payload"]["vcf"]["somatic_snv_vcf"]["counters"]["multi_allelic"] == 1
    assert run["payload"]["vcf"]["somatic_snv_vcf"]["counters"]["with_tvaf"] == 1
    assert run["payload"]["hla"]["resolved_alleles"] == ["HLA-A*29:02", "HLA-B*08:01"]
    assert store.load("mini", "ingest") == run["payload"]


def test_ingest_fails_hc_containment_check(tmp_path: Path) -> None:
    bundle = make_case(tmp_path, VCF_OUTSIDE, BED, ALLELES)
    run = run_ingest(bundle, tmp_path / "refs", ArtifactStore(tmp_path / "work"))
    check = next(c for c in run["payload"]["checks"] if "containment" in c["id"])
    assert check["status"] == "fail" and "chr2:100" in check["detail"]


def test_ingest_rejects_checksum_mismatch(tmp_path: Path) -> None:
    bundle = make_case(tmp_path, VCF_FULL, BED, ALLELES)
    tampered = json.loads(bundle.read_text())
    tampered["input_files"]["somatic_snv_vcf"]["sha256"] = "0" * 64
    bundle.write_text(json.dumps(tampered))
    with pytest.raises(InputIntegrityError, match="pinned checksum"):
        run_ingest(bundle, tmp_path / "refs", ArtifactStore(tmp_path / "work"))


def test_ingest_is_deterministic(tmp_path: Path) -> None:
    bundle = make_case(tmp_path, VCF_FULL, BED, ALLELES)
    r1 = run_ingest(bundle, tmp_path / "refs", ArtifactStore(tmp_path / "work-a"))
    r2 = run_ingest(bundle, tmp_path / "refs", ArtifactStore(tmp_path / "work-b"))
    assert r1["sha256"] == r2["sha256"]


def test_vcf_reader_rejects_malformed_row(tmp_path: Path) -> None:
    body = "chr1\tXX\t.\tA\tG\t.\tPASS\tSOMATIC"
    bundle = make_case(tmp_path, body, BED, ALLELES)
    refs = tmp_path / "refs"
    (refs / "seqc2" / "mini.vcf").write_text(VCF_HEADER + "\n" + body)
    import hashlib as _h
    b = json.loads(bundle.read_text())
    b["input_files"]["somatic_snv_vcf"]["sha256"] = _h.sha256((refs / "seqc2" / "mini.vcf").read_bytes()).hexdigest()
    bundle.write_text(json.dumps(b))
    with pytest.raises(VcfFormatError, match="non-integer POS"):
        list(VcfReader(refs / "seqc2" / "mini.vcf").records())


def test_bed_index_containment(tmp_path: Path) -> None:
    bed = tmp_path / "x.bed"
    bed.write_text("chr1\t100\t200\nchr1\t300\t400\n")
    idx = IntervalIndex.from_bed(bed)
    assert idx.n_intervals == 2
    assert idx.contains("chr1", 101)          # inside first (1-based)
    assert idx.contains("chr1", 200)          # end base (BED half-open end 200)
    assert not idx.contains("chr1", 201)
    assert not idx.contains("chr1", 201, 399)  # start outside rejects the range
    assert idx.contains("chr1", 301, 400)
    assert not idx.contains("chrY", 150)
    wide = tmp_path / "wide.bed"
    wide.write_text("chr1\t0\t300\n")
    wide_idx = IntervalIndex.from_bed(wide)
    assert wide_idx.contains("chr1", 298, 300)   # span ends on interval edge
    assert not wide_idx.contains("chr1", 298, 301)  # span crossing the edge


def test_case_bundle_validation(tmp_path: Path) -> None:
    bundle = make_case(tmp_path, VCF_FULL, BED, ALLELES)
    good = json.loads(bundle.read_text())
    _ = load_case_bundle(bundle)
    bad_build = dict(good, genome_build="hg19")
    bad_path = bundle.with_name("bad_build.json")
    bad_path.write_text(json.dumps(bad_build))
    with pytest.raises(CaseBundleError, match="genome_build"):
        load_case_bundle(bad_path)
    bad_pair = dict(good, pairing={"tumour": {}})
    bad_path2 = bundle.with_name("bad_pair.json")
    bad_path2.write_text(json.dumps(bad_pair))
    with pytest.raises(CaseBundleError, match="pairing"):
        load_case_bundle(bad_path2)


def test_allele_normalisation_honesty() -> None:
    statuses = normalise_alleles(["HLA-A*29:02", "a*02:01", "HLA-B*08:01", "HLA-B*08:01"])
    assert [s.normalised for s in statuses] == ["HLA-A*29:02", "HLA-A*02:01", "HLA-B*08:01", None]
    dup_note = statuses[-1].note
    assert dup_note and "duplicate" in dup_note
