# Judge-Facing Demonstration Script

One continuous research journey, every number on screen coming from the
computed artifacts (no disconnected mock data). Recommended 6-8 minutes.

## Step 0 - Orientation (30s)
"Every cancer is molecularly individual. NEOFORGE runs a real computational
pipeline on a real tumour-normal pair and shows the science of how
researchers go from sequencing data to investigated candidate targets."
Quick non-claims banner, then a pipeline map that doubles as dashboard nav.

## Step 1 - Tumour fingerprint (Research Explorer)
- Tumour vs normal molecular profile: mutation burden, class distribution,
  high-confidence-region coverage. Question answered: **what makes this
  tumour different from the matching normal?**

## Step 2 - Mutation explorer
- Interactive mutation map (positions harvested from the parsed truth VCF),
  filterable by consequence class, expressed genes first.
- Click a mutation - detail panel: transcript, affected protein, evidence tags
  (measured / predicted / missing).

## Step 3 - Protein view
- The transcript exons with codon-level zoom on the mutation; reference vs
  mutant peptide windows rendered from the codon-exact derivation stage.
- Label: conceptual diagram of exon structure + exact sequence track - clearly
  not a 3-D structure rendering.

## Step 4 - Neoantigen candidates
- All peptide windows for this mutation x patient HLA alleles x lengths
  8-11, sorted by MHCflurry prediction; predicted affinity, processing and
  presentation scores shown with the **prediction** tag and their allele.
- Reference-vs-mutant pairing: peptides whose reference counterpart also binds
  are flagged (tumour-specific windows highlighted).

## Step 5 - Prioritisation & evidence ledger
- Tiers A/B/C from the documented rule; each candidate shows the full evidence
  ledger: which features were measured, which predicted, which missing, and the
  exact threshold rule that placed it in its tier. Judges can toggle any rule
  and watch tier membership update - the ranking is inspectable, not a black box.

## Step 6 - Uncertainty & limitation register
- Per-candidate uncertainty: absent VAF, expression source, prediction-only
  status. Persistent non-claims banner, visually distinct from data.

## Step 7 - Analysis Workbench (any screen, one click away)
- Same data, full provenance: input files with checksums, stage parameters,
  software versions, per-stage timings, artifact hashes, determinism status.

## Step 8 - Research report
- Exportable Markdown research summary: candidates, evidence, limitations,
  provenance, and recommended next investigations for qualified researchers.

## Step 9 - Construct (design layer)
- A research-only in-silico polyepitope design generated from the tier-A
  shortlist: a proportional construct map (signal peptide, epitopes, AAY
  linkers, flagged junctions), the full selection ledger (included epitopes
  and every exclusion with its documented reason), the junction audit (which
  8-11-mers spanning the joins create newly flagged strong binders across
  case alleles, and how ordering improved), the back-transcription stats
  (codon table hash, GC content, homopolymer and restriction-site checks), and
  the secondary-structure propensity strip: every residue of the design
  labelled from pinned published Chou-Fasman (1978) parameters - a
  sequence-derived annotation, never a 3-D model, shown with its own
  "validation PENDING" state.
- Question answered: **what would the construct look like, reductively - and
  what does the software refuse to claim about it?** The artifact, FASTA and
  screen are all stamped "research design, not a manufacturable or validated
  sequence".

## Step 10 - Rosie pathway mapping (narrative close)
- The reported dog-vaccine case mapped, step by reported step, to NEOFORGE
  stages - each row stating what the software does, what it does not do, and
  where human expertise and the laboratory take over.
- Mentions the honest gaps: human-HLA vs canine-DLA (no canine predictor is
  integrated, and faithful canine cases are refused rather than imputed),
  protein structure handled only as sequence-derived propensity annotation
  (no 3-D modelling, no structure validation yet), and everything downstream of design
  (review, manufacture, administration) is outside the software by rule.

## What is deliberately NOT in the demo
- Any treatment or vaccine claim; any immunogenicity prediction; any
  cross-patient data blending. TESLA is shown only as audit evidence for the
  ranking rule at cohort level (L4), never as validation of the case; the
  benchmark page states its files, licensing and falsification clause verbatim.
- No claim that the designed construct is manufacturable, expressible,
  trafficked or immunogenic; the construct stage is a design layer whose
  outputs cannot be elevated to clinical statements by any toggle or export.
