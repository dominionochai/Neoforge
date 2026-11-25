# NEOFORGE

Personalised neoantigen **research-support** prototype: an original, reproducible
computational pipeline that takes a documented tumour-normal case, derives
mutant peptides from the reference genome and gene model, scores predicted MHC-I
presentation with a published open-source predictor, and prioritises candidates
with a transparent, evidence-backed rule set that judges can inspect.

## Non-claims (read first)

- NOT a diagnostic or clinical system; makes no treatment recommendations.
- Predicted MHC binding/presentation is **prediction**, not experiment.
- No immunogenicity is claimed for any candidate; experimental immunogenicity
  data (TESLA) is handled as an external cohort-level audit of the ranking rule
  only (executed: TTIF 0.541 vs random 0.200, AUPRC 0.153 vs prevalence 0.061;
  see `docs/EVALUATION_PLAN.md`).
- Missing evidence is shown as missing. Nothing is imputed or invented.

## Architecture

- `pipeline/` — `neoforge` Python package; stages `ingest → variants →
  translate → present → prioritize → report`; immutable content-hashed JSON
  artifacts with full provenance (`docs/SCIENTIFIC_DESIGN.md`). Includes
  `neoforge/maf.py` (GDC masked-MAF ingest with a single-patient barcode gate)
  and `neoforge/evaluation/tesla.py` (L4 benchmark harness).
- `frontend/` — zero-dependency static Research Explorer + Analysis Workbench
  (vanilla JS, hand-rolled SVG). It renders **only** projections of the hashed
  pipeline artifacts exported by `python -m neoforge.export_ui`; it computes no
  science of its own.
- `cases/` — per-case bundles pinning every input by SHA-256
  (`hcc1395` SEQC2 truth-set case, `skcm_d9a3z1` TCGA-SKCM open-tier case).
- `fetch/` — reference-data fetcher with pinned SHA-256 checksums
  (`fetch/manifests.csv`).
- `docs/` — scientific design, evaluation plan, demonstration script, data
  provenance.

## Quickstart

```bash
cd pipeline && pip install -e ".[dev]" && python -m pytest tests -q
python -m neoforge.cli --version
```

Case data is refreshed (with checksum verification) via
`fetch/fetch_manifests.py --out-dir <dir>` for sources whose checksums pin at
fetch time (GENCODE genome/GTF, HLA CSV).

End-to-end demo (pipeline stages → exports → re-export verified identical to
the tracked projections, modulo per-run wall-clock timings → serve the
interface):

```bash
./scripts/run_demo.sh          # uses /workspace/refs by default, override with NEOFORGE_REFS
# open http://localhost:8000
```

## The four evidence layers

| Layer | Content | Validation |
|---|---|---|
| L1 variant processing | somatic calls, filters, VAF | SEQC2 genomic truth-set case |
| L2 peptide derivation | codon-exact mutant peptides | unit tests + UniProt checks + truth-set exercise |
| L3 predicted presentation | MHCflurry scores | correct-use tests; predictor cited, not ours |
| L4 ranking-rule audit | TESLA experimental labels | executed — cohort level, never case-level |

`hcc1395` runs the full chain (L1+L2+L3). `skcm_d9a3z1` demonstrates the honest
boundary: no same-patient HLA genotype exists in the open GDC tier, so the
chain mechanically stops after L2 and the interface says so instead of
inventing predictions.

See `docs/EVALUATION_PLAN.md` for metrics, baselines and falsification
criteria, and `frontend/README.md` for the interface's data contract.
