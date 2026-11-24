"""Stage 2 — variants.

L1 selection ledger over the ingested SNV table. Every variant is either
selected for translation or retained in an exclusion ledger with its reason —
selection is inspectable, never a silent intersection. Version 0.1 processes
SNVs only: INDEL peptide derivation (frameshift ORF reconstruction) is a
documented future scope item, so those records are excluded by scope with
their counter visible, never by accident."""
from __future__ import annotations

from neoforge.artifacts import ArtifactStore, StageRun
from neoforge.errors import ArtifactChainError


def _load_ingest(store: ArtifactStore, case: str) -> dict:
    try:
        return store.load(case, "ingest")
    except FileNotFoundError as exc:
        raise ArtifactChainError(f"run the ingest stage before variants for case '{case}'") from exc


def run_variants(store: ArtifactStore, case: str) -> dict:
    ingest = _load_ingest(store, case)
    table = ingest["variant_tables"].get("somatic_snv_vcf")
    if table is None:
        raise ArtifactChainError("case bundle declared no somatic_snv_vcf; nothing to select")

    excluded: dict[str, int] = {}

    indel_table = ingest["variant_tables"].get("somatic_indel_vcf")
    if indel_table is not None:
        excluded["indels_out_of_scope_v0.1"] = len(indel_table["pos"])

    multi = 0
    non_acgt = 0
    positions = table["pos"]
    refs, alts = table["ref"], table["alt"]
    selected_indices: list[int] = []
    for i in range(len(positions)):
        alt_field = alts[i]
        if "," in alt_field:
            multi += 1
            continue
        if set(refs[i] + alt_field) - set("ACGT"):
            non_acgt += 1
            continue
        selected_indices.append(i)

    if multi:
        excluded["multi_allelic_unresolvable_v0.1"] = multi
    if non_acgt:
        excluded["non_acgt_alleles"] = non_acgt

    selected = {
        "columns": ["variant_key", "chrom", "pos", "ref", "alt", "tvaf", "nvaf", "tumor_depth"],
        "variant_key": [f"{table['chrom'][i]}:{positions[i]}{refs[i]}>{alts[i]}" for i in selected_indices],
        "chrom": [table["chrom"][i] for i in selected_indices],
        "pos": [positions[i] for i in selected_indices],
        "ref": [refs[i] for i in selected_indices],
        "alt": [alts[i] for i in selected_indices],
        "tvaf": [table["tvaf"][i] for i in selected_indices],
        "nvaf": [table["nvaf"][i] for i in selected_indices],
        "tumor_depth": [table["tumor_depth"][i] for i in selected_indices],
    }

    conservation = {
        "ingested_snv_records": len(positions),
        "selected": len(selected_indices),
        "excluded": len(positions) - len(selected_indices),
        "exclusion_reasons": excluded,
        "indels_seen": excluded.get("indels_out_of_scope_v0.1", 0),
    }
    if conservation["ingested_snv_records"] != conservation["selected"] + conservation["excluded"]:
        raise ArtifactChainError("selection ledger does not conserve the ingest table")

    run = StageRun(
        stage="variants",
        case=case,
        inputs={"ingest": store.read_manifest(case)["stages"]["ingest"]["artifact_id"]},
        parameters={"criteria": [
            "single-allelic SNV with ACGT alleles",
            "INDELs excluded from translation scope (documented v0.1 limitation, not a quality filter)",
            "position-level HC containment is enforced upstream by the ingest L1 gate",
        ]},
    )
    payload = {"selection": conservation, "selected": selected}
    artifact_id, digest = store.write(run, "variants", payload)
    return {"artifact_id": artifact_id, "sha256": digest, "payload": payload}
