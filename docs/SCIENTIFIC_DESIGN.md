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
construct = research-only design layer over the L1–L3 chain (no science added)
structure = residue-level propensity annotation of the design (pinned published parameters; not a 3-D model)
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
- **construct** — a **design layer, not a prediction layer**: builds a
  research-only in-silico polyepitope design from tier-A candidates (ledgered
  selection with allele-coverage and gene-diversity preferences, MHC-I linker
  joins, a junction audit that flags newly created strong binders, and
  back-translation with a pinned codon-usage table). It adds no biological
  claim; its artifact and FASTA are stamped "research design, not a
  manufacturable or validated sequence", and it mechanically refuses without
  the full L1–L3 chain, a resolved HLA genotype, or a surviving tier-A
  candidate. What a design layer can and cannot be evaluated on is specified
  in the construct section of `docs/EVALUATION_PLAN.md`.
- **structure** — an **annotation layer over the construct design**: labels
  every residue of the designed protein with a secondary-structure propensity
  (alpha-helix / beta-sheet / coil) computed from the published, pinned
  Chou–Fasman (1978) conformational parameters via a fixed, documented
  windowed-mean/argmax rule. It is derived arithmetic on the sequence — not a
  3-D model (no coordinates, contact maps or docking), not a measurement, and
  not compared to any experimental structure (validation pending, declared in
  the artifact). It refuses mechanically without a construct artifact.

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
- Correction (verified on the fetched files): the SEQC2 truth VCFs **do carry
  tumour VAF** in the `TVAF=` INFO field (alongside `NVAF` and `bwaDP` depths).
  The earlier note claiming "no tumour VAF" was wrong. The Tier-A VAF gate uses
  these measured values (and only these): never imputed, unavailable stays
  unavailable.
- Audited v1.2 release: all 39,447 SNV truth positions are contained in the
  high-confidence regions BED, but 23 of 1,625 INDEL positions fall outside it
  (e.g. chr1:48903925; first gap at [48903919, 48903940)) — recorded as an
  honest `warn` audit finding, not hidden.
- HCC1395/SEQC2 is a **genomic truth-set validation case only** (L1/L2). It is
  never described as immunogenicity validation.
