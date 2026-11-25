"""Stage 1 — ingest.

Loads one case bundle, byte-verifies every declared input against its pinned
SHA-256, and parses the somatic variation source: the SEQC2-style somatic VCFs
(independently re-checking the structural claim the source makes about itself —
v1.2.1: every call lies inside the v1.2 high-confidence regions) or the open
per-case GDC masked somatic MAF (TCGA cases); both parse into one column-major
variant table shape. Produces the L1 `ingest` artifact: the variant table plus
per-check statuses — the only thing downstream stages are allowed to consume.

VAF, when the source carries it, is kept with its source string; when the
source implies it (MAF depth/alt counts), it is computed by the single
documented transform recorded in the table; when unavailable, the field stays
missing with that reason. The HLA block is conditional: a bundle either pins
the same-sample HLA genotype (`hla_csv` input) or declares
`hla.status = "unavailable"` with its documented reason — presentation inputs
are never invented, and never borrowed from a different sample or source.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

from neoforge import hla
from neoforge import maf as maf_module
from neoforge.artifacts import ArtifactStore, StageRun
from neoforge.bed import IntervalIndex
from neoforge.case_config import load_case_bundle, resolve_input
from neoforge.errors import InputIntegrityError
from neoforge.vcf import VcfReader

_HC_CONTAINMENT_MIN = 1.0
_VARIANT_COLUMNS = ["chrom", "pos", "ref", "alt", "filters", "tvaf", "nvaf", "tumor_depth"]


def _parse_expression(refs_root: Path, inputs: dict, declare: dict) -> dict:
    from neoforge import expression  # noqa: PLC0415 — keeps other stages' dep surface small

    required = ("expression_tpm_log2p1", "ensembl_entrez_map")
    missing = [n for n in required if n not in inputs]
    if missing:
        raise InputIntegrityError(f"expression inputs missing: {missing}")
    ensg_by_entrez, ambiguous = expression.crosswalk(resolve_input(refs_root, inputs["ensembl_entrez_map"]))
    depmap_id, tpm_by_gene, counts = expression.gene_tpm(
        resolve_input(refs_root, inputs["expression_tpm_log2p1"]), ensg_by_entrez)
    return {
        "depmap_id": depmap_id,
        "cell_line": declare.get("cell_line"),
        "matrix": declare.get("matrix"),
        "source": {"path": inputs["expression_tpm_log2p1"]["path"],
                   "sha256": inputs["expression_tpm_log2p1"]["sha256"],
                   "transform": "tpm = 2^v - 1 (source stores log2(TPM+1))"},
        "source_description": "ingest artifact (measured DepMap TPM)",
        "crosswalk": {"path": inputs["ensembl_entrez_map"]["path"],
                      "sha256": inputs["ensembl_entrez_map"]["sha256"],
                      "ambiguous_entrez_dropped": ambiguous},
        "tpm_by_gene": tpm_by_gene,
        "counts": {**counts, "ambiguous_entrez": ambiguous},
    }


def _parse_gdc_expression(refs_root: Path, inputs: dict, declare: dict) -> dict:
    from neoforge import expression  # noqa: PLC0415

    tpm_by_gene, counts = expression.gene_tpm_gdc_star(
        resolve_input(refs_root, inputs["expression_star_tsv"]))
    return {
        "patient_sample": declare.get("patient_sample"),
        "matrix": declare.get("matrix"),
        "source": {"path": inputs["expression_star_tsv"]["path"],
                   "sha256": inputs["expression_star_tsv"]["sha256"],
                   "transform": "none; tpm_unstranded copied verbatim from the GDC STAR output"},
        "source_description": "ingest artifact (measured GDC STAR TPM, same patient)",
        "tpm_by_gene": tpm_by_gene,
        "counts": counts,
    }


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
        alts.append(",".join(rec.alts))
        filter_tokens.append(";".join(rec.filters))
    table = {
        "columns": _VARIANT_COLUMNS,
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

    checks: list[dict] = [
        {"id": "pinned_sha256", "status": "pass", "detail": f"{len(verified_hashes)} inputs match pinned digests"}
    ]

    hc_index = (IntervalIndex.from_bed(resolve_input(refs_root, inputs["hc_regions_bed"]))
                if "hc_regions_bed" in inputs else None)
    vaf_source = ("INFO/TVAF" if "somatic_snv_vcf" in inputs
                  else "computed: t_alt_count / t_depth" if "somatic_maf" in inputs
                  else "none")
    run = StageRun(
        stage="ingest",
        case=bundle["case_id"],
        inputs=verified_hashes,
        parameters={
            "hc_containment_min": _HC_CONTAINMENT_MIN,
            "vaf_source_field": vaf_source,
        },
    )

    per_vcf: dict[str, dict] = {}
    tables: dict[str, dict] = {}
    parsed_maf: dict | None = None
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
        counters_clean = {k: v for k, v in counters.items() if k != "outside_hc_regions"}
        counters_clean["outside_hc_region_positions"] = len(counters["outside_hc_regions"])
        per_vcf[name] = {
            "file_format": meta.file_format,
            "reference_declared": meta.reference,
            "n_header_lines": meta.n_header_lines,
            "counters": counters_clean,
        }
        tables[name] = table

    if "somatic_maf" in inputs:
        parsed_maf = maf_module.read_gdc_maf(resolve_input(refs_root, inputs["somatic_maf"]))
        full = parsed_maf["table"]
        snv = {c: [] for c in _VARIANT_COLUMNS}
        indel = {c: [] for c in _VARIANT_COLUMNS}
        for i, kind in enumerate(full["filters"]):
            target = snv if kind == "SNP" else indel
            for c in _VARIANT_COLUMNS:
                target[c].append(full[c][i])
        tables["somatic_maf"] = {**snv, "columns": _VARIANT_COLUMNS,
                                 "null_encoding": full["null_encoding"],
                                 "tvaf_source": full["tvaf_source"]}
        if any(indel["chrom"]):
            tables["somatic_indel_maf"] = {**indel, "columns": _VARIANT_COLUMNS,
                                           "null_encoding": full["null_encoding"],
                                           "tvaf_source": full["tvaf_source"]}
        pairing = bundle["pairing"]
        declared_samples = {"tumour": pairing["tumour"].get("sample"),
                            "normal": pairing["normal"].get("sample")}
        if declared_samples["tumour"] and declared_samples["tumour"] != parsed_maf["tumour_sample"]:
            raise InputIntegrityError(
                f"MAF tumour sample {parsed_maf['tumour_sample']!r} != bundle pairing {declared_samples['tumour']!r}")
        if declared_samples["normal"] and declared_samples["normal"] != parsed_maf["normal_sample"]:
            raise InputIntegrityError(
                f"MAF matched normal {parsed_maf['normal_sample']!r} != bundle pairing {declared_samples['normal']!r}")
        checks.append({
            "id": "somatic_maf.same_patient",
            "status": "pass",
            "detail": f"single tumour barcode {parsed_maf['tumour_sample']} + single matched normal "
                      f"{parsed_maf['normal_sample']}; counters {parsed_maf['counters']}",
        })

    hla_payload: dict
    if "hla_csv" in inputs:
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
        hla_payload = {
            "cell_line": allele_source.cell_line,
            "allele_statuses": [s.__dict__ for s in allele_statuses],
            "resolved_alleles": [s.normalised for s in allele_statuses if s.wellformed],
            "source": {"path": inputs["hla_csv"]["path"], "sha256": allele_source.sha256},
        }
    else:
        if bundle["hla"].get("status") != "unavailable" or not bundle["hla"].get("reason"):
            raise InputIntegrityError(
                "case bundle without hla_csv must declare hla.status='unavailable' with a reason")
        checks.append({
            "id": "hla.unavailable",
            "status": "warn",
            "detail": f"no same-sample HLA genotype input: {bundle['hla']['reason']} "
                      "(presentation/prioritisation stages will be refused rather than imputed)",
        })
        hla_payload = {"resolved_alleles": [], "status": "unavailable",
                       "reason": bundle["hla"]["reason"]}

    expression_block = None
    if "expression_tpm_log2p1" in inputs:
        expression_block = _parse_expression(refs_root, inputs, bundle.get("expression", {}))
        counts = expression_block["counts"]
        checks.append({
            "id": "expression.source",
            "status": "pass",
            "detail": (
                f"measured TPM for {expression_block['depmap_id']} "
                f"({expression_block['cell_line']}): {len(expression_block['tpm_by_gene'])} genes; "
                f"columns unmapped to Ensembl: {counts['unmapped_entrez']}; "
                f"ambiguous crosswalk columns dropped: {counts['ambiguous_entrez']}; "
                "single documented transform tpm = 2^(v) - 1"
            ),
        })
    elif "expression_star_tsv" in inputs:
        expression_block = _parse_gdc_expression(refs_root, inputs, bundle.get("expression", {}))
        counts = expression_block["counts"]
        checks.append({
            "id": "expression.source",
            "status": "pass",
            "detail": (
                f"measured GDC STAR TPM ({expression_block['patient_sample']}): "
                f"{len(expression_block['tpm_by_gene'])} genes; value gaps: {counts['value_gaps']}; "
                f"duplicate (PAR) gene ids resolved to their primary record: {counts['duplicate_gene_id']}; "
                "no transform applied"
            ),
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
        "hla": hla_payload,
        "expression": expression_block,
        "maf": ({"counters": parsed_maf["counters"],
                 "tumour_sample": parsed_maf["tumour_sample"],
                 "normal_sample": parsed_maf["normal_sample"]} if parsed_maf else None),
        "checks": checks,
    }
    artifact_id, digest = store.write(run, "ingest", payload)
    return {"artifact_id": artifact_id, "sha256": digest, "payload": payload}
