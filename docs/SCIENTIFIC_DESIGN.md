# Scientific Design

Status: foundation document. Every claim in this file is either a lab-verified fact
(sandbox run recorded in `DATA_PROVENANCE.md`) or explicitly marked **pending**.

## 1. Scope and non-claims

NEOFORGE is a computational research-support prototype for personalised neoantigen
research. It is **not**: a diagnostic tool, a clinical decision system, a claim of
therapeutic efficacy, or an autonomous vaccine design system. Nothing the system
outputs is experimentally demonstrated unless it is labelled as coming from a
published experimental dataset.

## 2. The four evidence layers — deliberately never merged

The pipeline produces four strictly separated layers. Every stage that touches a
candidate labels each number with its layer. The interface never presents a
higher layer as if it were validated by a lower one.

| Layer | Name | What it is | Validation status |
|---|---|---|---|
| **L1** | Variant processing | Somatic SNV/indel calls, filters, VAF where available | Input-integrity tests against SEQC2 truth artifacts (genomic truth-set validation case) |
| **L2** | Peptide derivation | Transcript-resolved, codon-exact mutant peptide windows (8–11 mers) | Unit tests + UniProt cross-checks + HCC1395 truth-set exercise |
| **L3** | Predicted MHC binding & predicted presentation | MHCflurry 2.4 affinity / processing / presentation scores (predictions, not measurements) | Correct-use tests (allele normalisation, check known-ligands sanity benchmark); the predictor itself is cited work, not our claim |
| **L4** | Experimental immunogenicity | T-cell-validated labels from published experiments (TESLA) | Used **only** as an external benchmark for the ranking rule; TESLA evaluation marked **pending** until fetch+licence finalisation is complete |

Rules enforced in code and UI:

1. Score namespaces are prefixed (`mhcflurry_`), never bare (`affinity`), so no
   number can be mistaken for a measured quantity.
2. `predicted` and `measured` are distinct enum values on every evidence item;
   the UI renders them with different visual treatments.
3. A candidate can never carry an L4 label. Validated immunogenicity lives only in
   the benchmark module operating on external data.

## 3. Pipeline stages

```
ingest → variants → translate → present → prioritize → report
```

Each stage is a pure function from artifacts to artifacts. Every artifact is
immutable, content-hashed (SHA-256), and stored with a provenance record
(inputs and their hashes, environment versions, parameters, execution timestamps).

- **ingest** — loads a documented case bundle (somatic VCF/MAF + HLA genotype +
  optional expression). Validates schema, genome build agreement, and HLA allele
  resolvability. Missing inputs (e.g. tumour VAF when the source truth set has
  none) are recorded as explicitly missing, never imputed.
- **variants** — applies documented somatic-variant criteria
  (Variant_Classification ∈ {Missense_Mutation, In_Frame_Ins, In_Frame_Del,
  Frame_Shift_*, Nonsense_Mutation} in protein-coding transcripts, transcript
  support level filter, normal-read checks when available). Every excluded
  mutation is retained with its exclusion reason so the filtering is inspectable.
- **translate** — derives mutant and reference peptide windows from the genome
  sequence and GENCODE gene model. Codon-exact: genomic coordinates → CDS
  position → codon → protein position (not string tricks on HGVSp).
- **present** — scores all mutant (and paired reference) peptide windows against
  the patient's HLA-I alleles with MHCflurry (affinity, processing,
  presentation). MHCflurry is used unmodified; we wrap, never re-implement it.
- **prioritize** — transparent rule-based tiers (see §4). No hidden weighted sum.
- **report** — the research summary: tiered candidates, evidence ledger,
  uncertainty register, provenance, and "what science would confirm/refute this".

## 4. Prioritisation: explainable, TESLA-grounded, no invented weights

Candidate evidence = for each (mutation, peptide length, allele):

- `mhcflurry_affinity` (nM), `mhcflurry_presentation_score`,
  `mhcflurry_presentation_percentile` — **L3 predictions**
- transcript expression (TPM) when available — measured cohort-level quantity
- tumour allele fraction when available — measured quantity
- antigen-binding reference peptide affinity (change vs reference — a mutated
  peptide must bind *better than* the reference fragment for a tumour-specific
  window, following pVACtools' documented practice)

Ranking rule (fixed, published in code and UI): tier assignment by documented
thresholds from the field (pVACtools, TESLA): strong binder < 500 nM;
expression TPM ≥ 1 when available; VAF ≥ 0.3 when available. Tiers

- **Tier A**: strong predicted binder (any allele) AND expressed AND (VAF ≥ 0.3 or VAF unavailable with explicit flag)
- **Tier B**: strong predicted binder, partial supporting evidence
- **Tier C**: candidate with all other mutations, ranked below, retained for inspection

Ordering inside a tier is by ascending best predicted affinity (a documented,
single, citable criterion — no composite score). The rationale: TESLA's validated
features are binding strength, tumour VAF, expression, and clonality; we use
binding as the ordering key and the others as tier gates with per-missing-feature
flags. This can never masquerade as an immunogenicity score.

## 5. Artifact and provenance system

- `work/<case>/<stage>.json` + `work/<case>/manifest.json`
- manifest: stage inputs (artifact id + sha256), parameters, software versions,
  wall time. A rerun with identical inputs must produce hash-identical outputs
  (checked by the deterministic run test).
- The web app reads artifacts read-only; both the Research Explorer and the
  Analysis Workbench bind to the same files.

## 6. HCC1395 facts the system must respect

- HCC1395: triple-negative breast-cancer cell line; HCC1395BL: matched
  EBV-immortalised B-lymphoblast line from the same donor (SEQC2 consortium).
- Published HLA-I genotype (OpenVax/TCLP): `HLA-A*29:02` (single A allele —
  HLA-A LOH documented in the cell line), `HLA-B*08:01`, `HLA-B*45:01`,
  `HLA-C*06:02`, `HLA-C*07:01`.
- SEQC2 truth VCFs carry caller-consensus INFO but **no tumour VAF**; VAF will be
  reported as unavailable for this case unless BAM pileups are added.
- HCC1395/SEQC2 is a **genomic truth-set validation case only** (L1/L2). It is
  never described as immunogenicity validation.
