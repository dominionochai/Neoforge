# Evaluation Plan

Principle: every layer is evaluated only against evidence of its own kind.
What would falsify each layer is stated explicitly.

## L1 — Variant processing validation (input integrity)

**Case:** SEQC2 HCC1395/HCC1395BL (genomic truth-set case).
**Method:** our ingest must reproduce, losslessly and in sorted order, the
published truth VCF records within `High-Confidence_Regions_v1.2.bed`, and the
GDC open MAF for the context case must round-trip through our parser.
**Metrics:** record-coverage = 100% of truth records ingested; 0 silent drops;
every exclusion categorised with a machine-checkable reason.
**Falsifies:** any coordinate/reference/alt mismatch, lost indel representation,
or unparsed record classes would falsify the ingest contract.

## L2 — Peptide derivation validation

**Method:**
1. Synthetic mini-genome fixtures (3 exons, + and − strands, splice junctions,
   frameshift, stop-gain): exact expected peptide sets asserted in unit tests.
2. Real-coordinate spot tests: for hand-selected well-known mutations
   (e.g., BRAF V600E, KRAS G12D, TP53 R175H), the derived mutant 9-mer containing
   position 600/12/175 must match the published amino-acid context of the
   protein (UniProt cross-check; runs only in slow tests with network).
3. HCC1395 truth-set exercise: every truth SNV mapped to an expressed,
   protein-coding ENST must yield exactly the peptides a second independent
   derivation (codon method on the exonic CDS) produces. Reported metric:
   discrepancy count (target 0) and flagging (not discarding) of any edge cases
   (near-splice, UTR-overlap, frameshift boundaries).
**Falsifies:** any codon/strand/coordination mismatch, e.g. a derived mutant
peptide that does not contain the substituted amino acid at the asserted offset,
or a reference peptide differing from the reference protein at the same offset.

## L3 — Predicted MHC binding / presentation

**Method:** we validate the *use* of the predictor, not the predictor:
1. Allele-name normalisation round trips (2/4-field, serology, TWT swaps).
2. Known-ligands sanity benchmark: published HLA-A*02:01 ligands
   (e.g., `GILGFVFTL` influenza M1, `SLYNTVATL` HIV, `NLVPMVATV` CMV pp65)
   must rank in the strong-binder (< 500 nM) class through our wrapper; a
   random 50-mer background must show the expected strong-binder base rate.
   Lab-verified so far: GILGFVFTL → 25.8 nM (A*02:01), SLYNTVATL → 82.6 nM (C*07:01).
3. Documentation test: the registry pins MHCflurry version + model release; the
   version and model release are recorded in every report artifact.
**Stated limitation:** NetMHCpan 4.1 is the field's strongest reference
predictor but is not freely reproducible (registration + academic license).
MHCflurry (Cell Systems 2018; Cell Systems 2020) is used with full citation.
**Falsifies:** a known ligand scored weak or an unnormalised allele yielding
silently different results across runs.

## L4 — Immunogenicity: external benchmark of the ranking layer (executed — see status below)

**Data:** TESLA (Wells et al., *Cell* 2020) — file `mmc4.xlsx`
(`master-bindings-selected`): 608 candidate pMHC across 6 melanoma patients,
**37 experimentally validated immunogenic** (tetramer/T-cell assays; the
`VALIDATED` column); plus `mmc7.xlsx` (independent 310-peptide TIL cohort,
4 validated). Files and SHA-256 recorded in `DATA_PROVENANCE.md`. Ingest code +
licence note + evaluation harness are implemented: `neoforge/evaluation/tesla.py`,
CLI stage `tesla_eval` (artifact case `tesla`, kept strictly separate from any
project case). The harness refuses unsupported alleles rather than imputing
predictions, and records its falsification verdict either way.
**Metrics:** TESLA-standard TTIF (top-20 immunogenic fraction), fraction ranked,
AUPRC; our ranking (using only features our pipeline computes for its own cases:
predicted binding, available measured features) vs. two baselines:
(a) random ranking, (b) strongest published team baseline (reported in
`mmc5.xlsx` as reference context only — our pipeline must not be compared
favourably to itself by using features TESLA teams didn't have access to;
all features we use for TESLA ranking must be recomputable from public data).
**Expected limitation & honesty note:** n=37 positives → wide confidence
intervals; results are cohort-level evidence about the *ranking rule*, not
patient-level claims.
**Falsification:** if our ranking does not beat random ranking (e.g., TTIF at or
below the random baseline's expectation), the prioritisation design is
falsified and must be revised or presented as exploratory with the failure shown.

## Stage 7 — construct: design-layer verification (in-silico design, never a product claim)

**What the stage is:** a research-only design layer over the L1–L3 chain. It
generates a candidate polyepitope mRNA *research design* from Tier A
candidates, with every choice ledgered. It adds no biological prediction of
its own: a designed construct being selected, transcribed or translated in a
cell is not simulated and not measured. The output is stamped
"research design, not a manufacturable or validated sequence" in the artifact,
the FASTA header and the interface.

**What is objectively verifiable in a design layer** (the metrics below are
the full set — no expression/immunogenicity metric exists for a construct
that was never made):

- **Round-trip exactness:** translating the designed CDS must reproduce the
  designed polyepitope protein residue-for-residue with zero internal stop
  codons. The stage raises (`ArtifactChainError`) otherwise — the artifact can
  never record a construct whose own translation disagrees with its design.
- **Determinism:** identical chain artifacts + pinned inputs (codon-usage
  table sha256 recorded in the artifact; fixed MHCflurry version recorded for
  the junction audit) → hash-identical construct artifact.
- **Mechanical refusals:** no HLA resolved, or the chain stopping before
  `present`, must refuse — construction is never run on partial or imputed
  evidence. The junction predictor's allele-coverage check refuses rather
  than silently skipping unsupported alleles.
- **Ledger completeness:** every mutation in the tier pool appears in exactly
  one ledger side, each entry with a documented reason (`not_tier_A`,
  `no_mutant_vs_wild_type_discrimination`, `capacity_reached`,
  `overlaps_selected_window`, `transcript_offset_unresolved`).
- **Junction audit honesty:** every junction-spanning 8–11-mer is scored
  against every case allele; a new strong binder (≤ the shared 500 nM
  threshold) created by a linker join is *flagged, never averaged away*. The
  adjacent-swap minimisation is auditable via its history and must be
  monotone (an accepted ordering can never worsen the flag count).

**Design assumptions (declared, not validated):** AAY linkers assumed
proteolytically permissive (processing of linker-embedded epitopes is not
measured); signal peptide routes through the ER (tPA leader, UniProt P00750);
UTRs and poly(A) are labelled placeholders; N1-methylpseudouridine and lipid
formulation are not simulated; proteasomal processing and epitope transport
are not simulated. The full register is in the artifact (`design_assumptions`,
`uncertainty_register`).

**Falsification:** for the design layer itself, falsification is binary and
mechanical — a round-trip mismatch, an unflagged strong junction binder, an
ambiguous selection ledger, or a non-reproducible hash would each falsify the
stage's verification claims and are test-pinned. What *cannot* be claimed from
this stage, ever: expression, folding, trafficking, immunogenicity or safety of
a designed sequence. Those require wet-lab data we do not have and do not
simulate.

## Compute determinism

Full rerun-from-artifacts test: identical inputs → identical SHA-256s on all
artifact outputs (documented seeds, fixed tool versions).

## Status of this plan (prototype results recorded here)

Executed against the real HCC1395/SEQC2 case (branch state at this commit):

- **L1 ingestion/audit:** 39,447 SNV truth positions 100% contained in the v1.2
  high-confidence regions (its own release claim, independently re-verified);
  23/1,625 INDEL truth positions fall outside those regions — recorded as an
  honest `warn`, never hidden (the release README does not mention this).
- **L2 derivation:** all **41,189/41,189** derived reference windows plus mutant
  windows equal slices of the GENCODE v46 official protein translations
  (independent check in `test_real_case.py`), including exon-junction codons on
  both strands. Zero unmatched, zero unknown transcripts.
- **L3 prediction:** MHCflurry 2.4.0 (models release 2.3.0) run per-allele on
  41,189 windows x 5 HLA alleles = 205,945 paired mutant/reference prediction
  rows (~164 s). Wrapper-only use; no scores re-implemented.
- **L1+L3 ranking (prioritize):** 329 candidate mutations → 245 Tier B (strong
  predicted binder), 84 Tier C, **0 Tier A** — because the expression source is
  not wired into this run yet, the Tier-A gate records `missing`, which can
  never silently pass. VAF gate evaluated on measured SEQC2 `TVAF` values.
- **L4 immunogenicity (ranking-rule audit on TESLA):** executed on the real
  benchmark — 608 pMHC (37 VALIDATED) from 6 patients / 13 class-I alleles,
  scored by the same unmodified per-allele MHCflurry wrapper (608×13 pairs where
  applicable; ~14 s). Ranking = ascending predicted affinity. Results:
  **TTIF top-20% = 0.541** vs random expectation 0.200 (2.7× enrichment;
  oracle headroom 1.0); **AUPRC = 0.153** vs random expectation 0.061;
  median rank of validated = 116 of 608; `rank_rule_falsified = False`
  recorded inside the `tesla_eval` artifact (`7ef67c3674c2566f`, including rank@k curve and rule ablation). Honesty
  notes: n=37 positives → wide confidence intervals; this is cohort-level
  evidence about the *ranking rule*, never a patient-level or candidate-level
  immunogenicity claim; candidate-level L4 evidence remains absent by design.
  The report stage embeds this reference (`tesla_cohort_reference`) only when
  the tesla_eval artifact exists in the store, and its uncertainty register
  reports the audit as `info` — verdict included verbatim. Additional
  falsification-adjacent outputs, both on the SAME labels (no new score):
  * **rank@k curve** — fraction of validated found in the top k% of the
    ranking (k = 1, 2, 5, 10, 20, 30, 50) vs the analytic random expectation
    k/100; the k=20 point equals the TTIF headline. Enrichment exceeds random
    at every k ≥ 10 (e.g. 78.4% found in the top 30% vs 30% random).
  * **rule ablation** — the same candidates re-ranked by *descending
    presentation score*: TTIF 0.595, AUPRC 0.210, median rank 111 — i.e. on
    this benchmark the presentation-score ordering slightly OUTPERFORMS the
    published affinity rule. Recorded honestly in the artifact, CLI output and
    interface (`rule_ablation`); the affinity rule is not falsified, but this
    is a documented candidate improvement for future iterations, not a hidden
    result.
- **Second case (TCGA-SKCM `skcm_d9a3z1`), L1 ingestion + L2 derivation:**
 executed end-to-end on open GDC data (588 SNVs, 3 INDELs from the ensemble
 masked MAF; same-patient barcode gate = pass; hla.unavailable warns and the
 present stage then refuses mechanically — no imputed presentation input).
 Independent L1↔L2 cross-check: every GDC-annotated missense HGVSp for this
 patient compared against our genome+GENCODE-v46-derived windows — **338/353
 agree on the exact (ref AA, position, alt AA) on the MAF's own transcript**,
 10 more agree exactly on another derived transcript of the same gene, 5
 differ only in protein position with identical substitution (annotation-
 version transcript-model drift, recorded, not hidden) and **0 contradictions
 of the substitution itself**; encoded in `test_real_case.py::…skcm…hgvsp…`.
- **Stage 7 construct (design layer), executed on the real HCC1395 chain**
  (artifact `68a7ceb9cfe4eb31`, hammer-rerun verified hash-identical; SKCM
  refused mechanically — "chain stopped before present"): of 329 candidate
  mutations, 80 Tier A → **15 epitopes selected** (documented cap; slider
  defaults) with **314 ledgered exclusions (249 not_tier_A + 55
  capacity_reached + 10 no_mutant_vs_wild_type_discrimination)**; **5/5 case
  HLA alleles** covered (HLA-A\*29:02, B\*08:01, B\*45:01, C\*06:02, C\*07:01)
  across **15 distinct genes**; protein **205 aa**, back-translated CDS
  **618 nt**, GC **66.34%**, **0 internal stop codons**, round-trip translation
  exact (`translation_matches_design = true`); **1,535 junction-spanning
  (junction, peptide, allele) scoring rows** by MHCflurry 2.4.0 — **72 flagged
  new strong binders** after one accepted adjacent swap, down from 73 at the
  initial ordering (history: initial 73 → swapped positions 13/14 → 72 → "no
  improving adjacent swap"); **0 restriction sites** from the documented
  12-enzyme panel; 3 homopolymer runs of exactly 5 bases (within the
  documented ≤5 threshold, listed in the artifact); codon table derived from
  the pinned GENCODE v46 GRCh38 CDS (38,855,269 codons; stop TGA) with its
  sha256 embedded in the artifact. FASTA export stamped in the header.
- **Determinism:** six-stage CLI end-to-end test (`test_pipeline_e2e.py`) reruns
  the whole chain and requires hash-identical artifacts — passing on the real
  case, including a hash-identical report with its exportable Markdown.
