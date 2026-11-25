"""Stage 7 — construct: research-only in-silico polyepitope design.

Input is the prioritize artifact of a completed L1-L3 chain (Tier A
candidates, HLA-resolved). This is a design layer, not a prediction layer:
epitope selection, linker assembly, a junction audit against the same case
alleles used in stage 4, and a reverse translation from a pinned codon table.
Every artifact carries the stamp "research design, not a manufacturable or
validated sequence".

Mechanical refusals (never a degraded design):
- any missing chain stage (SKCM stops at L2 by construction) -> refusal;
- no resolved HLA alleles -> refusal;
- unsupported case allele in the junction predictor -> refusal.

Selection rules, each ledgered with the reason that applied:
- source pool: prioritize Tier A candidates only;
- one epitope per mutation, window chosen from stage 4 predictions;
- the mutant window must be a predicted strong binder (<= BINDER_NM, the same
  pVACtools/TESLA convention as stage 5) while its paired reference window
  must NOT be (>= BINDER_NM): no wild-type discrimination -> excluded;
- windows within MAX_EPITOPE_LEN residues of a selected window on the same
  transcript are excluded as overlapping;
- greedy fill maximises marginal HLA-allele coverage, then gene diversity,
  then affinity; capacity default 15 (brief cap 10-20).

Assembly junctions are audited by scoring every junction-spanning peptide
(lengths 8-11, including the linker) against every case allele with
MHCflurry; epitope order is then improved by a deterministic bounded
minimisation over adjacent swaps and single-epitope reinsertions of flagged
new binders. An order whose residual flags survive every accepted move is
recorded, not resolved: the minimisation is audit-grade, not an optimizer
claim. Reverse translation uses the pinned
GENCODE-v46-derived human codon table in neoforge/data (raw counts from the
same reference the translate stage reads) with a deterministic homopolymer
mitigation pass. UTRs and poly(A) are labelled placeholders; mRNA chemistry
(N1-methylpseudouridine and friends) is not simulated and is declared.
"""
from __future__ import annotations

import hashlib
import time
from pathlib import Path
from typing import Callable

from neoforge.artifacts import ArtifactStore, StageRun
from neoforge.errors import ArtifactChainError

BINDER_NM = 500.0          # strong-binder convention shared with stage 5
MAX_EPITOPE_LEN = 11       # MHC-I window lengths scored in stage 4
JUNCTION_LENGTHS = (8, 9, 10, 11)
LINKER = "AAY"             # published MHC-I polyepitope linker convention
DEFAULT_SIGNAL_PEPTIDE = "tPA_leader"
DEFAULT_TRAFFICKING_TAG = "none"
MAX_EPITOPES = 15          # brief cap 10-20
MAX_HOMOPOLYMER = 5        # reverse translation avoids longer same-base runs
CODON_TABLE_NAME = "codon_usage_human_GRCh38_GENCODE_v46_CDS.tsv"

DISCLAIMER = "research design, not a manufacturable or validated sequence"

# Verified sequences only; every entry cites its authoritative record. These
# are design assumptions the user may swap, not computed results.
SIGNAL_PEPTIDES: dict[str, dict] = {
    "tPA_leader": {
        "aa": "MDAMKRGLCCVLLLCGAVFVSP",
        "source": "UniProt P00750 (TPA_HUMAN) FT SIGNAL 1..22, entry v283",
        "note": "secretory leader routing the translated product through the ER; a route assumption, not an effect claim",
    },
}
TRAFFICKING_TAGS: dict[str, dict] = {
    "none": {
        "aa": None,
        "source": "parameter default",
        "note": "no tag: subcellular routing after ER entry is not simulated",
    },
    "ubiquitin_C76": {
        "aa": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "source": "UniProt P0CG48 (UBC_HUMAN) CHAIN 1..76, entry v139",
        "note": "proteasome-routing convention from DNA-vaccine literature; appended C-terminally, exact fusion/scission geometry not simulated",
    },
}

RESTRICTION_ENZYMES = {
    "EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT",
    "XhoI": "CTCGAG", "KpnI": "GGTACC", "SacI": "GAGCTC",
    "SalI": "GTCGAC", "NotI": "GCGGCCGC", "BsaI": "GGTCTC",
    "BsmBI": "CGTCTC", "SapI": "GCTCTTC", "AgeI": "ACCGGT",
}

NON_CLAIMS = (
    "construct is a research design, not a manufacturable or validated sequence",
    "predicted MHC binding is not demonstrated immunogenicity",
    "no clinical, efficacy or treatment claim of any kind is made",
    "mRNA chemistry (N1-methylpseudouridine, lipids, dose) is not simulated",
)


def load_codon_usage(path: str | Path | None = None) -> tuple[dict[str, list[tuple[str, int]]], str]:
    path = Path(path) if path else Path(__file__).resolve().parents[1] / "data" / CODON_TABLE_NAME
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    by_aa: dict[str, list[tuple[str, int]]] = {}
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        codon, aa, count = line.split("\t")
        by_aa.setdefault(aa, []).append((codon, int(count)))
    for aa in by_aa:
        by_aa[aa].sort(key=lambda t: (-t[1], t[0]))
    return by_aa, digest


def _stage_payload(store: ArtifactStore, case: str, name: str) -> dict:
    try:
        return store.load(case, name)
    except FileNotFoundError as exc:
        raise ArtifactChainError(
            f"stage '{name}' has no artifact for case '{case}'; construction "
            f"requires the full L1-L3 chain ({exc})") from exc


def _refuse(case: str, ingest: dict, manifest: dict) -> None:
    stages = manifest.get("stages", {})
    missing = [s for s in ("ingest", "variants", "translate", "present", "prioritize")
               if s not in stages]
    if missing:
        raise ArtifactChainError(
            f"construct refuses case '{case}': chain stopped before {missing[0]}; "
            "no degraded construction without the full L1-L3 chain")
    if not ingest.get("hla", {}).get("resolved_alleles"):
        raise ArtifactChainError(
            f"construct refuses case '{case}': ingest resolved no HLA alleles; "
            "HLA-restricted design must not be invented")


def _col_rows(columns: dict) -> list[dict]:
    cols = columns["columns"]
    return [dict(zip(cols, vals)) for vals in zip(*[columns[c] for c in cols])]


def select_epitopes(candidates: list[dict], predictions_by_variant: dict[str, list[dict]],
                    offsets_by_key: dict[tuple[str, str], int],
                    *, max_epitopes: int, binder_nm: float,
                    case_alleles: list[str]) -> tuple[list[dict], list[dict]]:
    """Greedy allele-coverage then gene-diversity then affinity fill.

    Returns (included, excluded) ledger rows; every mutation gets an entry."""
    included: list[dict] = []
    excluded: list[dict] = []
    used_alleles: set[str] = set()
    used_genes: set[str] = set()
    occupied: list[tuple[str, int]] = []

    tier_a = [c for c in candidates if c["tier"] == "A"]
    for c in candidates:
        if c["tier"] != "A":
            excluded.append({"variant_key": c["variant_key"], "scope": "mutation",
                             "reason": "not_tier_A", "detail": f"tier={c['tier']}"})
    tier_a.sort(key=lambda c: (c["mhcflurry_affinity_nm"], c["variant_key"]))

    for cand in tier_a:
        if len(included) >= max_epitopes:
            excluded.append({"variant_key": cand["variant_key"], "scope": "mutation",
                             "reason": "capacity_reached", "detail": f"cap {max_epitopes}"})
            continue
        key = cand["variant_key"]
        windows = predictions_by_variant.get(key, [])
        valid = [w for w in windows
                 if w["mhcflurry_affinity_nm"] <= binder_nm
                 and w["reference_mhcflurry_affinity_nm"] >= binder_nm]
        if not valid:
            excluded.append({"variant_key": key, "scope": "mutation",
                             "reason": "no_mutant_vs_wild_type_discrimination",
                             "detail": "no window with mutant <= %.0f nM and reference >= %.0f nM"
                                       % (binder_nm, binder_nm)})
            continue
        scored = []
        for w in valid:
            offset = offsets_by_key.get((key, w["transcript_id"]))
            if offset is None:
                continue
            overlaps = [(tid, pos) for tid, pos in occupied
                        if tid == w["transcript_id"] and abs(pos - offset) < MAX_EPITOPE_LEN]
            if overlaps:
                continue
            gain_alleles = 0 if w["allele"] in used_alleles else 1
            gain_gene = 0 if w["gene_id"] in used_genes else 1
            scored.append((-gain_alleles, -gain_gene, w["mhcflurry_affinity_nm"],
                           w["peptide"], w["allele"], w, offset, gain_alleles, gain_gene))
        if not scored:
            reason = "overlaps_selected_window" if _near_occupied(valid, offsets_by_key, key, occupied) \
                else "transcript_offset_unresolved"
            excluded.append({"variant_key": key, "scope": "mutation", "reason": reason,
                             "detail": "no admissible window remained after ledgered exclusions"})
            continue
        scored.sort(key=lambda t: t[:5])
        _, _, _, _, _, w, offset, gain_alleles, gain_gene = scored[0]
        occupied.append((w["transcript_id"], offset))
        used_alleles.add(w["allele"])
        used_genes.add(w["gene_id"])
        included.append({
            "variant_key": key, "gene_id": w["gene_id"], "transcript_id": w["transcript_id"],
            "peptide": w["peptide"], "peptide_length": w["peptide_length"], "protein_offset": offset,
            "allele": w["allele"], "mutant_affinity_nm": w["mhcflurry_affinity_nm"],
            "reference_affinity_nm": w["reference_mhcflurry_affinity_nm"],
            "reference_peptide": w["reference_peptide"],
            "presentation_score": w["mhcflurry_presentation_score"],
            "selected_allele_new": bool(gain_alleles), "selected_gene_new": bool(gain_gene),
            "selection_reason": (f"tier A; mutant {w['mhcflurry_affinity_nm']:.1f} nM <= "
                                 f"{binder_nm:.0f}, reference {w['reference_mhcflurry_affinity_nm']:.1f} nM"
                                 + ("; new allele coverage" if gain_alleles else "")
                                 + ("; new gene" if gain_gene else "")),
        })
    return included, excluded


def _near_occupied(valid, offsets_by_key, key, occupied) -> bool:
    for w in valid:
        offset = offsets_by_key.get((key, w["transcript_id"]))
        if offset is None:
            continue
        if any(tid == w["transcript_id"] and abs(pos - offset) < MAX_EPITOPE_LEN
               for tid, pos in occupied):
            return True
    return False


def _selection_order(epitopes: list[dict]) -> list[int]:
    """Round-robin across genes (group order by best affinity), within a group
    by affinity: similar sequences are less likely to sit adjacent."""
    by_gene: dict[str, list[tuple[float, int]]] = {}
    for i, e in enumerate(epitopes):
        by_gene.setdefault(e["gene_id"], []).append((e["mutant_affinity_nm"], i))
    groups = sorted(by_gene.values(), key=lambda g: min(x[0] for x in g))
    for g in groups:
        g.sort(key=lambda t: (t[0], epitopes[t[1]]["variant_key"]))
    order: list[int] = []
    while any(groups):
        for g in groups:
            if g:
                order.append(g.pop(0)[1])
    return order


def junction_peptides(parts: list[str]) -> list[tuple[int, str]]:
    """Peptides of JUNCTION_LENGTHS spanning each junction, linker included.

    junction i is between parts[i] and parts[i+1]. A window must straddle the
    boundary and contain only standard residue letters (they do, by
    construction, for A/T/B gates enforced upstream).
    """
    out: list[tuple[int, str]] = []
    std = set("ACDEFGHIKLMNPQRSTVWY")
    for j in range(len(parts) - 1):
        context = parts[j][-MAX_EPITOPE_LEN:] + parts[j + 1][:MAX_EPITOPE_LEN]
        boundary = len(parts[j][-MAX_EPITOPE_LEN:])
        for ln in JUNCTION_LENGTHS:
            for start in range(max(0, boundary - ln + 1), min(boundary, len(context) - ln + 1)):
                seq = context[start:start + ln]
                if start < boundary < start + ln and set(seq) <= std:
                    out.append((j, seq))
    return out


def improve_order(order: list[int], audit_fn: Callable[[list[int]], tuple[int, list[dict]]],
                  *, move_family: str = "swap+reinsert", max_moves: int | None = None
                  ) -> tuple[list[int], int, list[dict]]:
    """Deterministic bounded minimisation of audit_fn's flagged-row count.

    move_family "swap+reinsert" offers two move kinds per round: an adjacent
    swap and a single-epitope reinsertion (remove at i, insert at j). Each
    accepted move strictly lowers the flagged count (monotone accepted
    history), the best improvement wins with first-enumerated tie-break, and
    the search stops on a no-improvement round or after max_moves accepted
    moves (default 2 x population). The result is therefore audit-grade: a
    documented, reproducible reduction of flagged junctions, with anything
    residual reported as residual.
    """
    if move_family not in ("swap", "swap+reinsert"):
        raise ValueError(f"unknown move_family '{move_family}'")

    def apply(cur: list[int], kind: str, i: int, j: int) -> list[int]:
        if kind == "swap":
            return cur[:i] + [cur[i + 1], cur[i]] + cur[i + 2:]
        trial = cur[:i] + cur[i + 1:]
        return trial[:j] + [cur[i]] + trial[j:]

    n = len(order)
    moves: list[tuple[str, int, int]] = [("swap", i, i + 1) for i in range(n - 1)]
    if move_family == "swap+reinsert":
        moves += [("move", i, j) for i in range(n) for j in range(n)
                  if j not in (i - 1, i, i + 1)]

    current = list(order)
    flagged, _rows = audit_fn(current)
    history: list[dict] = [{"action": "initial allele-coverage ordering",
                            "flagged_junction_peptides": flagged, "accepted": True}]
    budget = 2 * n if max_moves is None else max_moves
    while budget > 0:
        best_f, best_seq, best_action = None, None, None
        for kind, i, j in moves:
            trial = apply(current, kind, i, j)
            f, _rows = audit_fn(trial)
            if f < flagged and (best_f is None or f < best_f):
                best_f, best_seq, best_action = f, trial, (kind, i, j)
        if best_seq is None:
            history.append({"action": "no improving swap or reinsertion",
                            "flagged_junction_peptides": flagged, "accepted": False})
            break
        flagged, current = best_f, best_seq
        kind, i, j = best_action
        label = (f"swapped positions {i} and {i + 1}" if kind == "swap"
                 else f"moved epitope from position {i} to {j}")
        history.append({"action": label,
                        "flagged_junction_peptides": flagged, "accepted": True})
        budget -= 1
    return current, flagged, history


def reverse_translate(protein: str, by_aa: dict[str, list[tuple[str, int]]]) -> tuple[str, list[dict]]:
    """Most-frequent human codon per residue with a deterministic homopolymer fix:

    the default codon is replaced by the next-most-frequent encoding whenever
    it would extend a same-nucleotide run beyond MAX_HOMOPOLYMER. Returns
    (CDS without stop, decision ledger)."""
    cds_chars: list[str] = []
    decisions: list[dict] = []
    for idx, aa in enumerate(protein):
        options = by_aa.get(aa)
        if not options:
            raise ArtifactChainError(f"codon table has no codon for amino acid '{aa}'")
        chosen, note = options[0][0], ""
        for codon, _count in options:
            tail_override = max_run("".join(cds_chars[-4:]) + codon)
            if tail_override <= MAX_HOMOPOLYMER:
                chosen = codon
                if codon != options[0][0]:
                    note = (f"homopolymer mitigation: default codon {options[0][0]} "
                            f"would exceed a {MAX_HOMOPOLYMER}-base run")
                break
        cds_chars.append(chosen)
        decisions.append({"residue_index": idx, "aa": aa, "codon": chosen,
                          "default": chosen == options[0][0], "note": note})
    return "".join(cds_chars), decisions


def max_run(seq: str) -> int:
    best = run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        best = max(best, run)
    return best if seq else 0


def homopolymer_runs(seq: str, threshold: int) -> list[dict]:
    out = []
    i = 0
    while i < len(seq):
        j = i
        while j + 1 < len(seq) and seq[j + 1] == seq[i]:
            j += 1
        if j - i + 1 >= threshold:
            out.append({"base": seq[i], "start": i, "end": j, "length": j - i + 1})
        i = j + 1
    return out


def restriction_sites(seq: str) -> dict[str, dict]:
    out = {}
    for enzyme, site in sorted(RESTRICTION_ENZYMES.items()):
        positions = [i for i in range(len(seq) - len(site) + 1) if seq[i:i + len(site)] == site]
        if positions:
            out[enzyme] = {"site": site, "positions_0based": positions}
    return out


def _default_junction_factory():
    import mhcflurry  # noqa: PLC0415 - heavy import only on the real path

    predictor = mhcflurry.Class1PresentationPredictor.load()
    return predictor, mhcflurry.__version__, str(predictor.provenance_string)


def _translate(nt: str) -> str:
    from neoforge.stages.translate import CODON_TABLE  # noqa: PLC0415 - reuse standard code

    return "".join(CODON_TABLE[nt[i:i + 3]] for i in range(0, len(nt) - 2, 3))


def run_construct(store: ArtifactStore, case: str, *,
                  max_epitopes: int = MAX_EPITOPES,
                  linker: str = LINKER,
                  signal_peptide: str = DEFAULT_SIGNAL_PEPTIDE,
                  trafficking_tag: str = DEFAULT_TRAFFICKING_TAG,
                  binder_nm: float = BINDER_NM,
                  codon_table_path: str | Path | None = None,
                  junction_predictor_factory: Callable | None = None) -> dict:
    junction_predictor_factory = junction_predictor_factory or _default_junction_factory
    started = time.time()

    manifest = store.read_manifest(case)
    ingest = _stage_payload(store, case, "ingest")
    _refuse(case, ingest, manifest)
    prior = _stage_payload(store, case, "prioritize")
    present = _stage_payload(store, case, "present")
    translate = _stage_payload(store, case, "translate")

    if signal_peptide not in SIGNAL_PEPTIDES:
        raise ArtifactChainError(f"unknown signal peptide option '{signal_peptide}'")
    if trafficking_tag not in TRAFFICKING_TAGS:
        raise ArtifactChainError(f"unknown trafficking tag option '{trafficking_tag}'")

    candidates = _col_rows(prior["candidates"])
    predictions = _col_rows(present["predictions"])
    predictions_by_variant: dict[str, list[dict]] = {}
    for r in predictions:
        predictions_by_variant.setdefault(r["variant_key"], []).append(r)
    offsets_by_key = {(r["variant_key"], r["transcript_id"]): r["mutant_offset"]
                      for r in _col_rows(translate["derived"])}
    case_alleles = sorted({r["allele"] for r in predictions})

    included, excluded = select_epitopes(candidates, predictions_by_variant, offsets_by_key,
                                         case_alleles=case_alleles,
                                         max_epitopes=max_epitopes, binder_nm=binder_nm)
    if not included:
        raise ArtifactChainError("construct refuses: no epitope passed the documented rules")

    order = _selection_order(included)
    predictor, j_version, j_provenance = junction_predictor_factory()
    unsupported = [a for a in case_alleles if a not in set(predictor.supported_alleles)]
    if unsupported:
        raise ArtifactChainError(
            f"junction audit: alleles {unsupported} unsupported by the junction "
            "predictor; refusing to design")

    cache: dict[tuple[str, str], tuple[float, int]] = {}

    def parts_of(seq: list[int]) -> list[str]:
        epi = [included[i]["peptide"] for i in seq]
        return [p for pair in zip(epi, [linker] * len(epi)) for p in pair][:-1]

    def audit(seq: list[int]) -> tuple[int, list[dict]]:
        peptides: dict[int, list[str]] = {}
        for j, s in junction_peptides(parts_of(seq)):
            peptides.setdefault(j, []).append(s)
        uniq = sorted({s for lst in peptides.values() for s in lst})
        if not uniq:
            return 0, []
        rows: list[dict] = []
        for allele in case_alleles:
            missing = [s for s in uniq if (s, allele) not in cache]
            if missing:
                df = predictor.predict(missing, alleles=[allele], verbose=0)
                for peptide_value, affinity_value in zip(df["peptide"], df["affinity"]):
                    cache[(peptide_value, allele)] = (float(affinity_value), len(peptide_value))
        for j in sorted(peptides):
            for s in sorted(set(peptides[j])):
                for allele in case_alleles:
                    aff, ln = cache[(s, allele)]
                    rows.append({"junction": j, "peptide": s, "length": ln, "allele": allele,
                                 "mhcflurry_affinity_nm": aff, "flagged": aff <= binder_nm})
        return sum(1 for r in rows if r["flagged"]), rows

    current, flagged, history = improve_order(list(order), audit)
    flagged_final, junction_rows = audit(current)
    assert flagged_final == flagged, \
        "minimisation history disagrees with the final junction audit"

    final_epitopes = [included[i] for i in current]
    sig = SIGNAL_PEPTIDES[signal_peptide]
    tag = TRAFFICKING_TAGS[trafficking_tag]
    protein_parts = [sig["aa"]] + [epi["peptide"] for epi in final_epitopes]
    protein = protein_parts[0]
    for part in protein_parts[1:]:
        protein += linker + part
    if tag["aa"]:
        protein += linker + tag["aa"]

    by_aa, table_sha = load_codon_usage(codon_table_path)
    stop_codon = by_aa.get("*", [("TGA", 0)])[0][0]
    cds_nostop, codon_decisions = reverse_translate(protein, by_aa)
    cds = cds_nostop + stop_codon

    protein_check = _translate(cds_nostop)
    checks = {
        "translation_matches_design": protein_check == protein,
        "internal_stop_codons": sum(1 for aa in protein_check if aa == "*"),
        "cds_starts_with_atg": cds.startswith("ATG"),
        "stop_codon": stop_codon,
    }
    if not checks["translation_matches_design"]:
        raise ArtifactChainError("reverse translation does not reproduce the designed protein")

    mitigations = sum(1 for d in codon_decisions if d["note"])
    stats = {
        "cds_length_nt": len(cds),
        "protein_length_aa": len(protein),
        "gc_percent": round(100.0 * sum(cds.count(b) for b in "GC") / len(cds), 2),
        "homopolymer_runs_ge5": homopolymer_runs(cds, 5),
        "restriction_sites": restriction_sites(cds),
        "codon_decisions_mitigated": mitigations,
        "linker": linker, "linker_count": len(final_epitopes) - 1,
    }

    segments: list[dict] = [{
        "name": "5'UTR", "kind": "placeholder", "sequence": None,
        "note": "labelled placeholder: UTR is a vector/mRNA-platform decision, not simulated here",
    }, {
        "name": "signal peptide", "kind": "design_assumption", "aa": sig["aa"],
        "source": sig["source"], "note": sig["note"],
    }]
    for i, epi in enumerate(final_epitopes):
        segments.append({"name": f"epitope {i + 1}: {epi['peptide']}", "kind": "selected_epitope",
                         "variant_key": epi["variant_key"], "allele": epi["allele"],
                         "gene_id": epi["gene_id"], "aa": epi["peptide"],
                         "mutant_affinity_nm": epi["mutant_affinity_nm"]})
        if i < len(final_epitopes) - 1:
            segments.append({"name": "linker AAY", "kind": "linker", "aa": linker,
                             "source": "published MHC-I polyepitope linker convention"})
    if tag["aa"]:
        segments.append({"name": "trafficking tag", "kind": "design_assumption", "aa": tag["aa"],
                         "source": tag["source"], "note": tag["note"]})
    segments += [{
        "name": "3'UTR", "kind": "placeholder", "sequence": None,
        "note": "labelled placeholder: UTR is a vector/mRNA-platform decision, not simulated here",
    }, {
        "name": "poly(A)", "kind": "placeholder", "sequence": None, "length_nominal": 100,
        "note": "labelled placeholder: length and position are platform decisions, not simulated here",
    }]

    payload = {
        "stage": "construct",
        "case": case,
        "layer": "design (from L1-L3 chain)",
        "evidence_kind": "in_silico_design",
        "non_claims": list(NON_CLAIMS),
        "stamp": DISCLAIMER,
        "inputs": {name: {"artifact_id": manifest["stages"][name]["artifact_id"],
                          "sha256": manifest["stages"][name]["sha256"]}
                   for name in ("ingest", "variants", "translate", "present", "prioritize")},
        "parameters": {
            "max_epitopes": max_epitopes,
            "binder_nm": binder_nm,
            "threshold_source": "pVACtools strong-binder convention (shared with prioritize)",
            "linker": linker,
            "linker_source": "published MHC-I polyepitope linker convention (e.g. AAY)",
            "signal_peptide": signal_peptide, "signal_peptide_source": sig["source"],
            "signal_peptide_note": sig["note"],
            "trafficking_tag": trafficking_tag, "trafficking_tag_source": tag["source"],
            "trafficking_tag_note": tag["note"],
            "signal_peptide_options": {k: v["source"] for k, v in SIGNAL_PEPTIDES.items()},
            "trafficking_tag_options": {k: v["source"] for k, v in TRAFFICKING_TAGS.items()},
            "junction_lengths": list(JUNCTION_LENGTHS),
            "junction_flag_threshold_nm": binder_nm,
            "max_homopolymer_run": MAX_HOMOPOLYMER,
            "codon_table": {"file": CODON_TABLE_NAME, "sha256": table_sha,
                            "derivation": "raw CDS codon counts from pinned GENCODE v46 GRCh38 (scripts/build_codon_usage.py)"},
        },
        "design_assumptions": [
            "signal peptide routes the translated product through the ER; other routes are not modelled",
            "AAY linkers are assumed proteolytically permissive; processing of linker-embedded epitopes is not measured",
            "trafficking tag default 'none': routing after ER entry is an untested assumption listed as a swap option",
            "UTRs and poly(A) are labelled placeholders; no vector sequences are implied",
            "mRNA chemistry (N1-methylpseudouridine and lipid formulation) is not simulated",
        ],
        "case_summary": {"tier_a_candidates": sum(1 for c in candidates if c["tier"] == "A"),
                         "epitopes_selected": len(final_epitopes),
                         "epitopes_excluded_ledgered": len(excluded),
                         "case_alleles": case_alleles,
                         "alleles_covered": sorted({e["allele"] for e in final_epitopes}),
                         "genes_covered": sorted({e["gene_id"] for e in final_epitopes})},
        "protein": protein,
        "selection_ledger": {"included": final_epitopes, "excluded": excluded},
        "junction_audit": {
            "predictor": {"software": j_version, "provenance": j_provenance[:80]},
            "final_flagged_junction_peptides": flagged_final,
            "final_rows": junction_rows,
            "history": history,
        },
        "reverse_translation": {
            "aa_length": len(protein),
            "codon_choices_summary": (f"{mitigations} codon choices mitigated against the "
                                      f"{MAX_HOMOPOLYMER}-base homopolymer threshold"
                                      if mitigations else "no homopolymer mitigations needed"),
        },
        "nucleotide_stats": stats,
        "sequence_checks": checks,
        "uncertainty_register": [
            "every junction score is a prediction; none is measured",
            "processing of linker-embedded epitopes by the proteasome is not simulated",
            "codon choices maximise table frequency only; no codon-pair optimisation was performed",
            "subcellular routing after signal-peptide entry is an assumption, not a measurement",
            "no DLA (canine) handling exists in this human-allele pipeline",
        ],
        "nucleotide_map": segments,
    }

    run = StageRun(stage="construct", case=case,
                   inputs={"prioritize": manifest["stages"]["prioritize"]["artifact_id"],
                           "present": manifest["stages"]["present"]["artifact_id"],
                           "translate": manifest["stages"]["translate"]["artifact_id"]},
                   parameters={"max_epitopes": max_epitopes, "linker": linker,
                               "binder_nm": binder_nm, "signal_peptide": signal_peptide,
                               "trafficking_tag": trafficking_tag})
    run.record_software({"mhcflurry": j_version})
    artifact_id, digest = store.write(run, "construct", payload)

    fasta_path = store.case_dir(case) / f"construct_cds.{artifact_id}.fasta"
    fasta_path.write_text(
        f">{case}_{artifact_id} neoforge {DISCLAIMER}\n"
        + "\n".join(cds[i:i + 60] for i in range(0, len(cds), 60)) + "\n")
    return {"artifact_id": artifact_id, "sha256": digest, "payload": payload,
            "fasta_path": fasta_path.name, "wall": round(time.time() - started, 3)}
