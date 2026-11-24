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

## What is deliberately NOT in the demo
- Any treatment or vaccine claim; any immunogenicity prediction; any
  cross-patient data blending. TESLA is shown only as audit evidence for the
  ranking rule at cohort level (L4), never as validation of the case; the
  benchmark page states its files, licensing and falsification clause verbatim.
