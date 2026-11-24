# NEOFORGE

Personalised neoantigen **research-support** prototype: an original, repoducible
computational pipeline that takes a documented tumour-normal case, derives
mutant peptides from the reference genome and gene model, scores predicted MHC-I
presentation with a published open-source predictor, and prioritises candidates
with a transparent, evidence-backed rule set that judges can inspect.

## Non-claims (read first)

- NOT a diagnostic or clinical system; makes no treatment recommendations.
- Predicted MHC binding/presentation is **prediction**, not experiment.
- No immunogenicity is claimed for any candidate; experimental immunogenicity
  data (TESLA) is handled as an external benchmark only, currently **pending**.
- Missing evidence is shown as missing. Nothing is imputed or invented.

## Architecture

- `pipeline/` — `neoforge` Python package; stages `ingest → variants →
  translate → present → prioritize → report`; immutable content-hashed JSON
  artifacts with full provenance (`docs/SCIENTIFIC_DESIGN.md`).
- `app/` — FastAPI + artifact-driven single-page interface; Research Explorer
  and Analysis Workbench read the **same** artifacts.
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
`fetch/fetch_manifests.py --out-dir <dir> --allow-unpinned` for sources whose
checksums pin at fetch time (GENCODE genome/GTF, HLA CSV).

## The four evidence layers

| Layer | Content | Validation |
|---|---|---|
| L1 variant processing | somatic calls, filters, VAF | SEQC2 genomic truth-set case |
| L2 peptide derivation | codon-exact mutant peptides | unit tests + UniProt checks + truth-set exercise |
| L3 predicted presentation | MHCflurry scores | correct-use tests; predictor cited, not ours |
| L4 immunogenicity | T-cell-validated labels | TESLA external benchmark — **pending** |

See `docs/EVALUATION_PLAN.md` for metrics, baselines and falsification
criteria.
