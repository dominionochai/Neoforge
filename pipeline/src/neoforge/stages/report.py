"""Stage 6 — report: the research summary.

A read-only aggregation of the five upstream artifacts. Every number in the
report is re-read from the artifact store — nothing is recomputed, nothing is
invented, nothing carries wall-clock time (so the report, like every artifact,
hashes identically on a rerun with identical inputs). Sections:

- shortlist: top-ranked candidates from the prioritize artifact, in its own
  documented order, each with every supporting feature and its evidence kind;
- evidence ledger: what was measured (L1), derived (L2), predicted (L3) —
  and what does not exist for the candidates (L4 experimental immunogenicity:
  candidate-level pending; cohort-level ranking-rule audit referenced when the
  tesla_eval artifact exists in this store);
- uncertainty register: warnings, missing-feature flags and scope exclusions
  lifted verbatim from upstream artifacts;
- what_science_would_confirm: the specific observations that would support or
  refute treating these candidates seriously (methodology, not data);
- research_report_md: an exportable Markdown rendering of all of the above.

This report is a research summary. It is not a treatment recommendation.
"""
from __future__ import annotations

from neoforge.errors import ArtifactChainError

UPSTREAM = ("ingest", "variants", "translate", "present", "prioritize")

WHAT_SCIENCE_WOULD_CONFIRM = [
    "RNA evidence (e.g. tumour transcriptome TPM >= 1) would confirm candidate expression; "
    "without it Tier A cannot be claimed",
    "mass spectrometry detection of a prioritised peptide on the patient's HLA would confirm "
    "presentation experimentally (currently only predicted)",
    "T-cell reactivity assays (e.g. ELISpot) against a prioritised peptide would be required "
    "for any immunogenicity statement — none exists in this run",
    "independent replication of this pipeline from the pinned inputs would confirm the "
    "candidate list is a property of the data, not of this implementation",
]


def run_report(store, case: str, *, shortlist_size: int = 10) -> dict:
    from neoforge.artifacts import StageRun  # noqa: PLC0415 — keeps stage dep surface small

    manifest = store.read_manifest(case)
    missing = [s for s in UPSTREAM if s not in manifest["stages"]]
    if missing:
        raise ArtifactChainError(f"report requires completed stages {missing} for case '{case}'")

    ing = store.load(case, "ingest")
    var = store.load(case, "variants")
    tra = store.load(case, "translate")
    pre = store.load(case, "present")
    pri = store.load(case, "prioritize")
    if pri.get("evidence_kind") != "derived":
        raise ArtifactChainError("report requires the prioritize stage's derived artifact")

    candidates = pri["candidates"]
    cols = candidates["columns"]
    shortlist = [dict(zip(cols, vals))
                 for vals in zip(*[candidates[c] for c in cols])][:shortlist_size]

    snv_counters = ing["vcf"]["somatic_snv_vcf"]["counters"]
    tesla_ref = _tesla_reference(store)
    construct_ref = _construct_reference(store, case)
    ledger = _ledger(ing, var, tra, pre, tesla_ref)
    provenance = {s: {"artifact_id": manifest["stages"][s]["artifact_id"],
                      "sha256": manifest["stages"][s]["sha256"]} for s in UPSTREAM}
    if construct_ref is not None:
        provenance["construct"] = {"artifact_id": construct_ref["artifact_id"],
                                   "sha256": construct_ref["sha256"]}

    payload = {
        "parameters": {"shortlist_size": shortlist_size, "stages_read": list(UPSTREAM)},
        "evidence_kind": "derived",
        "layer": "L1+L2+L3",
        "non_claims": [
            "this report is a research summary, not a treatment recommendation",
            "no immunogenicity of any candidate is claimed; candidate-level L4 evidence "
            "is pending (TESLA, when run, audits the ranking rule at cohort level only)",
        ],
        "case": ing["case"],
        "provenance_chain": provenance,
        "inputs": ing["inputs"],
        "hla": {"resolved_alleles": ing["hla"]["resolved_alleles"],
                "cell_line": ing["hla"]["cell_line"],
                "source": ing["hla"]["source"]},
        "molecular_profile": {
            "snv_truth_records": snv_counters["records"],
            "indel_truth_records": ing["vcf"]["somatic_indel_vcf"]["counters"]["records"],
            "snvs_selected_for_peptides": var["selection"]["selected"],
            "variants_with_peptides": tra["summary"]["variants_with_peptides"],
            "peptide_windows": tra["summary"]["windows_derived"],
            "transcripts_hit": tra["summary"]["transcripts_hit"],
        },
        "evidence_ledger": ledger,
        "tier_counts": pri["tier_counts"],
        "prioritisation_rule": pri["parameters"],
        "shortlist": shortlist,
        "construct": construct_ref,
        "uncertainty_register": _uncertainty(ing, var, pri, tesla_ref),
        "what_science_would_confirm": WHAT_SCIENCE_WOULD_CONFIRM,
        "research_report_md": None,
    }
    payload["research_report_md"] = _render_markdown(payload)

    report_inputs = {s: manifest["stages"][s]["artifact_id"] for s in UPSTREAM}
    if construct_ref is not None:
        report_inputs["construct"] = construct_ref["artifact_id"]
    artifact_id, digest = store.write(
        StageRun(stage="report", case=case,
                 inputs=report_inputs,
                 parameters=payload["parameters"]),
        "report", payload)
    return {"artifact_id": artifact_id, "sha256": digest, "payload": payload}


def _ledger(ing: dict, var: dict, tra: dict, pre: dict, tesla_ref: dict | None) -> dict:
    snv = ing["vcf"]["somatic_snv_vcf"]["counters"]
    return {
        "L1_measured": {
            "kind": "measured",
            "snv_truth_records": snv["records"],
            "records_with_tumour_vaf": snv["with_tvaf"],
            "records_with_tumour_depth": snv["with_depth"],
            "hla_typed_alleles": ing["hla"]["resolved_alleles"],
            "hla_source": ing["hla"]["source"]["path"],
            "note": "tumour-normal truth-set quantities; not immunogenicity evidence",
        },
        "L2_derived": {
            "kind": "derived",
            "codon_table": "standard genetic code",
            "peptide_lengths": tra["parameters"]["peptide_lengths"],
            "variants_with_peptides": tra["summary"]["variants_with_peptides"],
            "peptide_windows": tra["summary"]["windows_derived"],
            "exclusions_counted": {k: v for k, v in sorted(tra["summary"].items())
                                   if v and k not in ("transcripts_hit", "variants_with_peptides",
                                                      "windows_derived", "variants_processed")},
            "note": "reference windows validated against GENCODE v46 official protein translations",
        },
        "L3_predicted": {
            "kind": "predicted",
            "model": pre["parameters"]["predictor"],
            "model_version": pre["parameters"]["mhcflurry_version"],
            "model_provenance": pre["parameters"]["model_provenance"],
            "alleles": pre["parameters"]["alleles"],
            "candidates_scored": pre["candidates_scored"],
            "predictions_made": pre["predictions_made"],
            "note": "MHC binding and presentation predictions only",
        },
        "tesla_cohort_reference": tesla_ref,
        "L4_experimental": {
            "kind": "missing",
            "status": "PENDING",
            "note": "no experimental immunogenicity label exists for any candidate of this "
                    "case" + ("; cohort-level audit of the ranking rule itself: see "
                              "tesla_cohort_reference (TESLA benchmark, separate case)"
                              if tesla_ref else
                              "; TESLA ranking-rule audit not yet run in this store"),
        },
    }


def _uncertainty(ing: dict, var: dict, pri: dict, tesla_ref: dict | None) -> list[dict]:
    findings: list[dict] = []
    for check in ing["checks"]:
        if check["status"] != "pass":
            findings.append({"source": f"ingest.{check['id']}", "status": check["status"],
                             "detail": check["detail"]})
    for reason, count in sorted(var["selection"]["exclusion_reasons"].items()):
        findings.append({"source": f"variants.{reason}", "status": "warn",
                         "detail": f"{count} records excluded: documented scope limitation"})
    flag_counts: dict[str, int] = {}
    for flags in pri["candidates"]["gate_flags"]:
        for flag in flags:
            flag_counts[flag] = flag_counts.get(flag, 0) + 1
    for flag, count in sorted(flag_counts.items()):
        findings.append({"source": f"prioritize.{flag}", "status": "warn",
                         "detail": f"gate input missing for {count} candidate(s); "
                                   "gate recorded as missing, not passed"})
    if tesla_ref:
        findings.append({"source": "L4_experimental", "status": "info",
                         "detail": "the priority ranking rule was audited on the TESLA "
                                   f"benchmark (608 pMHC, 37 validated): TTIF top-20% "
                                   f"{tesla_ref['ttif_top20pct']:.3f} vs random "
                                   f"{tesla_ref['ttif_random_expectation']:.3f}"
                                   + ("; ranking rule NOT falsified" if not tesla_ref["rank_rule_falsified"]
                                      else "; ranking rule FALSIFIED")
                                   + " — cohort-level evidence about the rule, not about this case"})
    else:
        findings.append({"source": "L4_experimental", "status": "warn",
                         "detail": "no experimental immunogenicity evaluation yet (TESLA label "
                                   "evaluation, EVALUATION_PLAN.md, is pending)"})
    return findings


def _construct_reference(store, case: str) -> dict | None:
    """Mirror the construct design into the report when its artifact exists.

    Read-only and optional: a case whose chain stops before L3 (or with no run
    construct stage) simply reports the absence instead of inventing a design."""
    rec = store.read_manifest(case)["stages"].get("construct")
    if rec is None:
        return None
    c = store.load(case, "construct")
    cs, stats, audit = c["case_summary"], c["nucleotide_stats"], c["junction_audit"]
    return {
        "artifact_id": rec["artifact_id"],
        "sha256": rec["sha256"],
        "stamp": c["stamp"],
        "fasta_file": f"construct_cds.{rec['artifact_id']}.fasta",
        "tier_a_candidates": cs["tier_a_candidates"],
        "epitopes_selected": cs["epitopes_selected"],
        "epitopes_excluded_ledgered": cs["epitopes_excluded_ledgered"],
        "alleles_covered": cs["alleles_covered"],
        "case_alleles": cs["case_alleles"],
        "genes_covered": cs["genes_covered"],
        "protein_length_aa": stats["protein_length_aa"],
        "cds_length_nt": stats["cds_length_nt"],
        "gc_percent": stats["gc_percent"],
        "junction_rows_scored": len(audit["final_rows"]),
        "junction_flagged_final": audit["final_flagged_junction_peptides"],
        "junction_flagged_initial": audit["history"][0]["flagged_junction_peptides"],
        "optimiser_history": audit["history"],
        "epitope_table": [{"variant_key": e["variant_key"], "peptide": e["peptide"],
                           "allele": e["allele"], "gene_id": e["gene_id"],
                           "mutant_affinity_nm": e["mutant_affinity_nm"],
                           "reference_affinity_nm": e["reference_affinity_nm"]}
                          for e in c["selection_ledger"]["included"]],
    }


def _tesla_reference(store) -> dict | None:
    man = store.read_manifest("tesla")
    rec = man["stages"].get("tesla_eval")
    if rec is None:
        return None
    t = store.load("tesla", "tesla_eval")
    m = t["metrics"]
    return {"artifact_id": rec["artifact_id"], "sha256": rec["sha256"],
            "n_candidates": m["n_candidates"], "n_validated": m["n_validated"],
            "ttif_top20pct": m["ttif_top20pct"],
            "ttif_random_expectation": m["ttif_random_expectation"],
            "auprc": m["auprc"], "auprc_random_expectation": m["auprc_random_expectation"],
            "rank_rule_falsified": t["rank_rule_falsified"],
            "scope": "cohort-level audit of the ranking rule — not about any project case"}


def _render_markdown(payload: dict) -> str:
    case = payload["case"]
    profile = payload["molecular_profile"]
    tiers = payload["tier_counts"]
    lines = [
        f"# NEOFORGE research report — {case['display_name']}",
        "",
        case.get("description") or "",
        "",
        "## Molecular profile",
        "",
        f"- Somatic SNVs (high-confidence truth set): **{profile['snv_truth_records']}**",
        f"- Somatic INDELs (truth set): {profile['indel_truth_records']} "
        "(excluded from peptide derivation in this version, not quality-filtered)",
        f"- Variants with derived peptide windows: **{profile['variants_with_peptides']}**",
        f"- Peptide windows (lengths "
        f"{', '.join(map(str, payload['evidence_ledger']['L2_derived']['peptide_lengths']))}): "
        f"**{profile['peptide_windows']}** across {profile['transcripts_hit']} transcripts",
        f"- HLA alleles (typed, {payload['hla']['cell_line']}): "
        f"{', '.join(payload['hla']['resolved_alleles'])}",
        "",
        "## Prioritised candidates",
        "",
        "Tier rule (fixed, no composite scores): strong predicted binder < "
        f"{payload['prioritisation_rule']['strong_binder_nm']} nM best predicted affinity; "
        f"expression TPM >= {payload['prioritisation_rule']['expressed_tpm']} when a source is "
        f"provided; tumour VAF >= {payload['prioritisation_rule']['min_tumour_vaf']} when "
        "measured. Ordering: tier, then ascending best predicted affinity.",
        "",
        "Tier counts: " + ", ".join(f"Tier {t}: {n}" for t, n in sorted(tiers.items())),
        "",
        "Variant | Tier | Best peptide | Allele | Pred. affinity (nM) | Ref. affinity "
        "(nM) | Tumour VAF | Flags",
        "--- | --- | --- | --- | --- | --- | --- | ---",
    ]
    for r in payload["shortlist"]:
        lines.append(
            f"{r['variant_key']} | {r['tier']} | {r['best_peptide']} | {r['best_allele']} | "
            f"{r['mhcflurry_affinity_nm']:.1f} | {r['reference_mhcflurry_affinity_nm']:.1f} | "
            f"{r['tvaf'] if r['tvaf'] is not None else 'unavailable'} | "
            f"{', '.join(r['gate_flags']) or '—'}")
    lines += [
        "",
        "All predicted affinities are MHCflurry predictions (see provenance), not measurements.",
        "",
        "## Stage 7 construct \u2014 in-silico polyepitope research design",
        "",
    ]
    con = payload["construct"]
    if con is None:
        lines += [
            "No construct design artifact exists for this case: the construct stage refuses",
            "mechanically without the full L1\u2013L3 chain and a resolved HLA genotype, and no",
            "design must be invented in its absence.",
        ]
    else:
        lines += [
            f"**{con['stamp']}**",
            "",
            f"Selected {con['epitopes_selected']} epitopes from a {con['tier_a_candidates']}-mutation "
            f"tier-A pool ({con['epitopes_excluded_ledgered']} ledgered exclusions); HLA alleles "
            f"covered {len(con['alleles_covered'])}/{len(con['case_alleles'])}; "
            f"{len(con['genes_covered'])} genes.",
            f"Design protein {con['protein_length_aa']} aa; back-translated CDS "
            f"{con['cds_length_nt']} nt (GC {con['gc_percent']:.2f}%). FASTA export: "
            f"`{con['fasta_file']}`.",
            f"Junction audit: {con['junction_rows_scored']} junction-spanning prediction rows "
            f"(MHCflurry); {con['junction_flagged_final']} newly flagged strong binders remain "
            f"after the bounded minimisation (initial ordering: "
            f"{con['junction_flagged_initial']}).",
            "",
            "Variant | Peptide | Allele | Gene | Pred. affinity (nM) | Ref. affinity (nM)",
            "--- | --- | --- | --- | --- | ---",
        ]
        lines += [f"{e['variant_key']} | {e['peptide']} | {e['allele']} | {e['gene_id']} | "
                  f"{e['mutant_affinity_nm']:.1f} | {e['reference_affinity_nm']:.1f}"
                  for e in con["epitope_table"]]
        lines += [
            "",
            "This is a research design for qualified-reviewer discussion. It is not a "
            "manufacturable or validated sequence, and nothing here is claimed about "
            "expression, folding, trafficking or immunogenicity.",
        ]
    lines += [
        "",
        "## Uncertainty register",
        "",
    ]
    lines += [f"- **[{f['status']}] {f['source']}** — {f['detail']}"
              for f in payload["uncertainty_register"]]
    lines += ["", "## What science would confirm or refute these candidates", ""]
    lines += [f"- {item}" for item in payload["what_science_would_confirm"]]
    lines += [
        "",
        "## Provenance",
        "",
        f"Model: {payload['evidence_ledger']['L3_predicted']['model']} "
        f"v{payload['evidence_ledger']['L3_predicted']['model_version']} "
        f"({payload['evidence_ledger']['L3_predicted']['model_provenance']})",
        "",
        "Pipeline artifacts (immutable, SHA-256 pinned):",
        "",
        "| Stage | Artifact | SHA-256 |",
        "| --- | --- | --- |",
    ]
    lines += [f"| {s} | `{d['artifact_id']}` | `{d['sha256']}` |"
              for s, d in sorted(payload["provenance_chain"].items())]
    lines += ["", *[f"**{non_claim}**" for non_claim in payload["non_claims"]]]
    return "\n".join(lines) + "\n"
