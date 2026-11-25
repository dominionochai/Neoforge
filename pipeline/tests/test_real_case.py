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

import csv
import gzip
import json
import re
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
REPO = Path(os.environ.get("NEOFORGE_REPO", Path(__file__).resolve().parents[2]))
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


SKCM_BUNDLE = REPO / "cases/skcm_d9a3z1/case.json"
SKCM_MAF = REFS / "tcga/d9a3z1.wxs.aliquot_ensemble_masked.maf.gz"

_THREE1 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
           "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
           "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*"}
_HGVP = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")


def test_real_skcm_ingest_same_patient_and_boundary(tmp_path: Path) -> None:
    store = ArtifactStore(tmp_path / "work")
    run = run_ingest(SKCM_BUNDLE, REFS, store)
    checks = {c["id"]: c["status"] for c in run["payload"]["checks"]}
    assert checks["somatic_maf.same_patient"] == "pass"
    assert checks["hla.unavailable"] == "warn"
    assert run["payload"]["expression"]["tpm_by_gene"]["ENSG00000141510"] > 0
    tables = run["payload"]["variant_tables"]
    assert len(tables["somatic_maf"]["pos"]) == 588
    assert len(tables.get("somatic_indel_maf", {}).get("pos", [])) == 3


def test_real_skcm_translation_matches_maf_hgvsp_annotations(tmp_path: Path) -> None:

    from neoforge.sequence import IndexedFasta  # noqa: PLC0415
    from neoforge.stages.translate import run_translate  # noqa: PLC0415

    store = ArtifactStore(tmp_path / "work")
    run_ingest(SKCM_BUNDLE, REFS, store)
    from neoforge.stages.variants import run_variants
    run_variants(store, "skcm_d9a3z1")
    fasta_gz = REFS / "gencode/GRCh38.primary_assembly.genome.fa.gz"
    genome = IndexedFasta(IndexedFasta.decompress_if_needed(fasta_gz, REFS / "gencode"))
    run_translate(store, "skcm_d9a3z1",
                  REFS / "gencode/gencode.v46.annotation.gtf.gz", genome)
    from neoforge.stages.present import run_present
    from neoforge.errors import ArtifactChainError
    with pytest.raises(ArtifactChainError, match="no HLA alleles"):
        run_present(store, "skcm_d9a3z1")

    derived = store.load("skcm_d9a3z1", "translate")["derived"]

    observed: dict[str, set] = {}
    for i in range(len(derived["variant_key"])):
        observed.setdefault(derived["variant_key"][i], set()).add(
            (derived["transcript_id"][i].split(".", 1)[0],
             derived["protein_pos"][i], derived["ref_aa"][i], derived["alt_aa"][i]))

    with gzip.open(SKCM_MAF, "rt") as f:
        lines = [line for line in f if not line.startswith("#")]
    reader = csv.reader(lines, delimiter="\t")
    header = next(reader)
    rows = [dict(zip(header, r)) for r in reader]

    expected_missense = [r for r in rows
                         if r["Variant_Type"] == "SNP" and r["Variant_Classification"] == "Missense_Mutation"]
    compared = same_transcript = 0
    for r in expected_missense:
        m = _HGVP.match(r["HGVSp"] or "")
        if not m:
            continue
        key = f"{r['Chromosome']}:{r['Start_Position']}{r['Reference_Allele']}>{r['Tumor_Seq_Allele2']}"
        windows = observed.get(key)
        if not windows:
            continue
        compared += 1
        want = (int(m.group(2)), _THREE1[m.group(1)], _THREE1[m.group(3)])
        if want in {(pos, ref, alt) for tid, pos, ref, alt in windows
                    if tid == r["Transcript_ID"]}:
            same_transcript += 1
        else:
            verified = {(ref, alt) for _, _, ref, alt in windows}
            assert (want[1], want[2]) in verified, (
                f"{key}: MAF HGVSp {r['HGVSp']} contradicts every derived window")
    assert compared >= 350
    assert same_transcript / compared > 0.90
