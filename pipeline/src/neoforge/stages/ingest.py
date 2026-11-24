"""Stage 1 — ingest.

Loads one case bundle, byte-verifies every declared input against its pinned
SHA-256, parses the somatic VCFs, and independently re-checks the structural
claims the source makes about itself (SEQC2 v1.2.1: every call lies inside the
v1.2 high-confidence regions). Produces the L1 `ingest` artifact: a
column-major variant table plus per-check statuses — the only thing downstream
stages are allowed to consume. VAF, when the source carries it, is kept with
its source string; when it does not, the field is missing with that reason."""
from __future__ import annotations

import hashlib
from pathlib import Path

from neoforge import hla
from neoforge.artifacts import ArtifactStore, StageRun
from neoforge.bed import IntervalIndex
from neoforge.case_config import load_case_bundle, resolve_input
from neoforge.errors import InputIntegrityError
from neoforge.vcf import VcfReader

_HC_CONTAINMENT_MIN = 1.0


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 21), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _verify_inputs(refs_root: Path, bundle: dict) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for name, entry in bundle["input_files"].items():
        path = resolve_input(refs_root, entry)
        if not path.exists():
            raise InputIntegrityError(f"input '{name}' missing at {path}")
        actual = sha256_file(path)
        if actual != entry["sha256"]:
            raise InputIntegrityError(
                f"input '{name}' failed pinned checksum: {actual} != {entry['sha256']}"
            )
        hashes[name] = actual
    return hashes


def _parse_vcf_table(reader: VcfReader, hc_index: IntervalIndex) -> tuple[dict, dict]:
    """Counters measure two containment notions separately: SEQC2's own claim
    ('all calls inside HC regions') refers to the call position; multi-base
    REF spans can legitimately cross an interval edge, so span containment is
    reported as an observation, not a gate."""
    chroms: list[str] = []
    positions: list[int] = []
    refs: list[str] = []
    alts: list[str] = []
    filter_tokens: list[str] = []
    tvaf: list[float | None] = []
    nvaf: list[float | None] = []
    depth: list[int | None] = []
    counters = {
        "records": 0,
        "biallelic_snv": 0,
        "non_snv": 0,
        "multi_allelic": 0,
        "non_acgt_alleles": 0,
        "outside_hc_regions": [],
        "ref_span_outside_hc_regions": 0,
        "with_tvaf": 0,
        "with_nvaf": 0,
        "with_depth": 0,
    }
    for rec in reader.records():
        counters["records"] += 1
        if len(rec.alts) > 1:
            counters["multi_allelic"] += 1
        elif rec.is_snv:
            counters["biallelic_snv"] += 1
        else:
            counters["non_snv"] += 1
        if any(set(x) - set("ACGT") for x in (rec.ref, *rec.alts)):
            counters["non_acgt_alleles"] += 1
        if not hc_index.contains(rec.chrom, rec.pos):
            counters["outside_hc_regions"].append(f"{rec.chrom}:{rec.pos}")
        elif len(rec.ref) > 1 and not hc_index.contains(rec.chrom, rec.pos, rec.pos + len(rec.ref) - 1):
            counters["ref_span_outside_hc_regions"] += 1
        tvaf_value = rec.info.get("TVAF")
        nvaf_value = rec.info.get("NVAF")
        if isinstance(tvaf_value, str):
            tvaf.append(float(tvaf_value))
            counters["with_tvaf"] += 1
        else:
            tvaf.append(None)
        if isinstance(nvaf_value, str):
            nvaf.append(float(nvaf_value))
            counters["with_nvaf"] += 1
        else:
            nvaf.append(None)
        dp = rec.info.get("bwaDP")
        if isinstance(dp, str):
            counters["with_depth"] += 1
            first = dp.split(",", 1)[0]
            depth.append(int(float(first)))
        else:
            depth.append(None)
        chroms.append(rec.chrom)
        positions.append(rec.pos)
        refs.append(rec.ref)
        alts.append(rec.alts[0])
        filter_tokens.append(";".join(rec.filters))
    table = {
        "columns": ["chrom", "pos", "ref", "alt", "filters", "tvaf", "nvaf", "tumor_depth"],
        "chrom": chroms,
        "pos": positions,
        "ref": refs,
        "alt": alts,
        "filters": filter_tokens,
        "tvaf": tvaf,
        "nvaf": nvaf,
        "tumor_depth": depth,
        "null_encoding": "null = value unavailable in source, not zero",
    }
    return table, counters


def run_ingest(bundle_path: Path, refs_root: Path, store: ArtifactStore) -> dict:
    bundle = load_case_bundle(bundle_path)
    inputs = bundle["input_files"]
    verified_hashes = _verify_inputs(refs_root, bundle)
    run = StageRun(
        stage="ingest",
        case=bundle["case_id"],
        inputs=verified_hashes,
        parameters={
            "hc_containment_min": _HC_CONTAINMENT_MIN,
            "vaf_source_field": "INFO/TVAF",
        },
    )

    checks: list[dict] = [
        {"id": "pinned_sha256", "status": "pass", "detail": f"{len(verified_hashes)} inputs match pinned digests"}
    ]

    hc_index = IntervalIndex.from_bed(resolve_input(refs_root, inputs["hc_regions_bed"]))
    per_vcf: dict[str, dict] = {}
    tables: dict[str, dict] = {}
    for name in ("somatic_snv_vcf", "somatic_indel_vcf"):
        if name not in inputs:
            continue
        reader = VcfReader(resolve_input(refs_root, inputs[name]))
        meta = reader.read_metadata()
        if meta.reference and "grch38" not in meta.reference.lower():
            raise InputIntegrityError(
                f"{name}: VCF declares reference {meta.reference!r}, bundle declares {bundle['genome_build']}"
            )
        table, counters = _parse_vcf_table(reader, hc_index)
        containment = 1.0 - len(counters["outside_hc_regions"]) / max(counters["records"], 1)
        checks.append({
            "id": f"{name}.hc_region_containment",
            "status": "pass" if containment >= _HC_CONTAINMENT_MIN else "warn",
            "detail": f"{containment:.4f} of {counters['records']} record positions inside HC regions; "
                      f"multi-base REF spans crossing an interval edge: {counters['ref_span_outside_hc_regions']} (observed, not gated)"
                      + (f"; outside-positions: {counters['outside_hc_regions'][:5]}" if counters["outside_hc_regions"] else ""),
        })
        per_vcf[name] = {
            "file_format": meta.file_format,
            "reference_declared": meta.reference,
            "n_header_lines": meta.n_header_lines,
            "counters": {k: v for k, v in counters.items() if k != "outside_hc_regions"},
        }
        tables[name] = table

    allele_source = hla.load_cell_line_alleles(
        resolve_input(refs_root, {"path": inputs["hla_csv"]["path"]}), bundle["hla"]["cell_line"]
    )
    allele_statuses = hla.normalise_alleles(list(allele_source.alleles))
    if not all(s.wellformed for s in allele_statuses):
        bad = [s.note for s in allele_statuses if not s.wellformed]
        raise InputIntegrityError(f"HLA genotype failed normalisation: {bad}")
    loci = sorted({s.locus for s in allele_statuses})
    checks.append({
        "id": "hla.loci",
        "status": "pass",
        "detail": f"loci present: {','.join(loci)}; retained verbatim from source, "
                  "missing loci/alleles are not imputed",
    })

    payload = {
        "case": {
            "case_id": bundle["case_id"],
            "display_name": bundle["display_name"],
            "description": bundle.get("description"),
            "genome_build": bundle["genome_build"],
            "pairing": bundle["pairing"],
            "pairing_provenance": bundle.get("pairing_provenance"),
        },
        "inputs": [
            {"name": name, "path": entry["path"], "sha256": verified_hashes[name], "role": entry["role"]}
            for name, entry in inputs.items()
        ],
        "vcf": per_vcf,
        "variant_tables": tables,
        "hla": {
            "cell_line": allele_source.cell_line,
            "allele_statuses": [s.__dict__ for s in allele_statuses],
            "resolved_alleles": [s.normalised for s in allele_statuses if s.wellformed],
            "source": {"path": inputs["hla_csv"]["path"], "sha256": allele_source.sha256},
        },
        "checks": checks,
    }
    artifact_id, digest = store.write(run, "ingest", payload)
    return {"artifact_id": artifact_id, "sha256": digest, "payload": payload}
