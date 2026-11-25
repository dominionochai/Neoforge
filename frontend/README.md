# NEOFORGE frontend

Zero-dependency static interface (vanilla HTML/CSS/JS, hand-rolled SVG). It renders
**only** projections of the pipeline's content-addressed artifacts and performs no
scientific computation of its own — the single interactive recomputation (the
priority-rule playground on `#/prioritize`) restates the published artifact rule and
self-checks against the artifact's tier counts every render.

## Run

```bash
python -m neoforge.export_ui --work-dir ../work --cases-dir ../cases --out frontend/data
python -m http.server 8000 -d frontend
# open http://localhost:8000
```

`data/` holds the generated per-case projections (`case_<id>.json`) and the TESLA
benchmark export (`benchmark_tesla.json`). Every block inside them carries the
upstream artifact id + SHA-256 in its `source` field, which the Analysis Workbench
(`#/workbench`) displays.
