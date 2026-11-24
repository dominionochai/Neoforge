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

## L4 — Immunogenicity: external benchmark of the ranking layer (PENDING)

**Data:** TESLA (Wells et al., *Cell* 2020) — file `mmc4.xlsx`
(`master-bindings-selected`): 608 candidate pMHC across 6 melanoma patients,
**37 experimentally validated immunogenic** (tetramer/T-cell assays; the
`VALIDATED` column); plus `mmc7.xlsx` (independent 310-peptide TIL cohort,
4 validated). Files and SHA-256 recorded in `DATA_PROVENANCE.md`. Ingest code +
licence note + evaluation harness are **pending** — marked so everywhere.
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
- **L4 immunogenicity:** still **PENDING** (TESLA ingest + metric harness not
  built). Marked pending in the report artifact and its uncertainty register.
- **Determinism:** six-stage CLI end-to-end test (`test_pipeline_e2e.py`) reruns
  the whole chain and requires hash-identical artifacts — passing on the real
  case, including a hash-identical report with its exportable Markdown.
