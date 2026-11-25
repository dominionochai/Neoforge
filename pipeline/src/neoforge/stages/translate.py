"""Stage 3 — translate: codon-exact mutant peptide derivation.

Given the selected SNV table, finds every protein-coding CDS overlap, assembles
the real codon from the reference genome (handling exon junctions and the minus
strand), and derives reference/mutant peptide windows of configurable lengths —
all windows of the mutant protein that contain the mutated residue, following
the published neoantigen-window convention. Genomic integrity gates:

- REF allele must equal the FASTA base (strand-complemented on minus strand);
- REF codon must assemble without N;
- 5'-incomplete CDS (nonzero first frame) is excluded entirely.

Every excluded variant is counted with its reason; nothing is silently dropped.
"""
from __future__ import annotations

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

_AAS = ("TTT F  TTC F  TTA L  TTG L  CTT L  CTC L  CTA L  CTG L  ATT I  ATC I  "
        "ATA I  ATG M  GTT V  GTC V  GTA V  GTG V  TCT S  TCC S  TCA S  TCG S  "
        "CCT P  CCC P  CCA P  CCG P  ACT T  ACC T  ACA T  ACG T  GCT A  GCC A  "
        "GCA A  GCG A  TAT Y  TAC Y  TAA *  TAG *  CAT H  CAC H  CAA Q  CAG Q  "
        "AAT N  AAC N  AAA K  AAG K  GAT D  GAC D  GAA E  GAG E  TGT C  TGC C  "
        "TGA *  TGG W  CGT R  CGC R  CGA R  CGG R  AGT S  AGC S  AGA R  AGG R  "
        "GGT G  GGC G  GGA G  GGG G ")
CODON_TABLE: dict[str, str] = {
    _AAS[i:i + 3]: _AAS[i + 4] for i in range(0, len(_AAS), 7)
}


def complement(base: str) -> str:
    return base.translate(_COMPLEMENT)


def run_translate(store, case: str, gtf_path, genome, *, peptide_lengths=(8, 9, 10, 11)) -> dict:
    from neoforge.gtf import GtfTranscriptIndex  # noqa: PLC0415 — keeps GTF load out of other stages

    variants = store.load(case, "variants")["selected"]
    index = GtfTranscriptIndex.from_gtf(gtf_path)

    stats: dict[str, int] = {}
    transcripts_hit: set[str] = set()

    cols = {
        "variant_key": [], "transcript_id": [], "gene_id": [], "strand": [],
        "protein_pos": [], "ref_aa": [], "alt_aa": [],
        "peptide_length": [], "window_start": [], "mutant_offset": [],
        "window_seq": [], "ref_window_seq": [],
    }
    protein_cache: dict[str, str] = {}

    def count(reason: str, n: int = 1) -> None:
        stats[reason] = stats.get(reason, 0) + n

    def protein_of(t, tid: str) -> str:
        seq = protein_cache.get(tid)
        if seq is None:
            parts = []
            for start, end in t.exons:
                part = genome.fetch(t.chrom, start - 1, end)
                if t.strand == "-":
                    part = part.translate(_COMPLEMENT)[::-1]
                parts.append(part)
            cds = "".join(parts)
            if len(cds) % 3:
                count("transcripts_partial_cds_3prime")
            seq = "".join(CODON_TABLE[cds[i:i + 3]] for i in range(0, len(cds) - 2, 3))
            protein_cache[tid] = seq
        return seq

    n_variants = len(variants["variant_key"])
    for i in range(n_variants):
        key, chrom, pos, ref, alt = (variants["variant_key"][i], variants["chrom"][i],
                                     variants["pos"][i], variants["ref"][i], variants["alt"][i])
        if len(ref) != 1 or len(alt) != 1 or "," in alt:
            count("skipped_non_single_allele")
            continue
        tids = index.covering_transcripts(chrom, pos)
        if not tids:
            count("no_protein_coding_overlap")
            continue
        genomic_base = genome.fetch(chrom, pos - 1, pos)
        if genomic_base != ref:
            count("ref_mismatch", len(tids))
            continue
        for tid in tids:
            t = index.transcripts[tid]
            if t.first_frame != 0:
                count("transcript_incomplete_5prime")
                continue
            off = t.offset_of(pos)
            if off is None:
                count("internal_offset_miss")
                continue
            codon_idx, codon_pos = off // 3, off % 3
            if codon_idx * 3 + 2 >= 3 * t.n_codons:
                count("variant_in_partial_final_codon")
                continue
            positions = t.positions_for_offsets([codon_idx * 3 + j for j in range(3)])
            bases = [genome.fetch(chrom, p - 1, p) for p in positions]
            if t.strand == "-":
                bases = [complement(b) for b in bases]
            ref_codon = "".join(bases)
            if "N" in ref_codon:
                count("ref_codon_gapped")
                continue
            alt_base = alt if t.strand == "+" else complement(alt)
            alt_codon = ref_codon[:codon_pos] + alt_base + ref_codon[codon_pos + 1:]
            aa_ref, aa_alt = CODON_TABLE[ref_codon], CODON_TABLE[alt_codon]
            if aa_alt == aa_ref:
                count("synonymous")
                continue
            if aa_alt == "*":
                count("stop_gained")
                continue
            if aa_ref == "*":
                count("stop_in_reference_cds")
                continue
            transcripts_hit.add(tid)
            protein = protein_of(t, tid)
            protein_len = len(protein)
            mutant_protein = protein[:codon_idx] + aa_alt + protein[codon_idx + 1:]
            one_based = codon_idx + 1
            for k in peptide_lengths:
                s_lo = max(1, one_based - k + 1)
                s_hi = min(protein_len - k + 1, one_based)
                for s in range(s_lo, s_hi + 1):
                    ref_window = protein[s - 1:s - 1 + k]
                    mut_window = mutant_protein[s - 1:s - 1 + k]
                    if "*" in ref_window or "*" in mut_window:
                        count("window_contains_stop")
                        continue
                    if "X" in ref_window or "X" in mut_window:
                        count("window_contains_ambiguous_residue")
                        continue
                    cols["variant_key"].append(key)
                    cols["transcript_id"].append(tid)
                    cols["gene_id"].append(t.gene_id)
                    cols["strand"].append(t.strand)
                    cols["protein_pos"].append(one_based)
                    cols["ref_aa"].append(aa_ref)
                    cols["alt_aa"].append(aa_alt)
                    cols["peptide_length"].append(k)
                    cols["window_start"].append(s)
                    cols["mutant_offset"].append(one_based - s + 1)
                    cols["window_seq"].append(mut_window)
                    cols["ref_window_seq"].append(ref_window)

    count("variants_processed", n_variants)
    count("variants_with_peptides", len(set(cols["variant_key"])))
    count("windows_derived", len(cols["window_seq"]))
    run_params = {"peptide_lengths": list(peptide_lengths), "codon_table": "standard"}
    from neoforge.artifacts import StageRun

    run = StageRun(
        stage="translate",
        case=case,
        inputs={"variants": store.read_manifest(case)["stages"]["variants"]["artifact_id"]},
        parameters=run_params,
    )
    payload = {
        "parameters": run_params,
        "summary": dict(sorted(stats.items())) | {
            "transcripts_hit": len(transcripts_hit),
        },
        "derived": {
            "columns": list(cols),
            **cols,
            "uniqueness_key": ["variant_key", "transcript_id", "peptide_length", "window_start"],
        },
    }
    artifact_id, digest = store.write(run, "translate", payload)
    return {"artifact_id": artifact_id, "sha256": digest, "payload": payload}
