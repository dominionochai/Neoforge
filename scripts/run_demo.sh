#!/usr/bin/env bash
# NEOFORGE one-command demo: pipeline stages -> exports -> identity
# verification -> local server.
#
# Every number shown in the interface is recomputed from the reference data in
# REFS (deterministic, artefact-hashed); the export step then verifies that the
# projections are identical to the tracked frontend/data/ files (modulo
# per-run wall-clock timings, which legitimately vary), so what you browse is
# provably what the pipeline produced.
#
# Usage:
#   ./scripts/run_demo.sh              # run everything, then serve on :8000
#   ./scripts/run_demo.sh --no-serve   # run everything, skip the server
#
# Environment:
#   NEOFORGE_REFS   reference-data directory (default /workspace/refs;
#                   build it first with fetch/fetch_manifests.py)
#   PORT            server port (default 8000)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REFS="${NEOFORGE_REFS:-/workspace/refs}"
PORT="${PORT:-8000}"
SERVE=1
if [[ "${1:-}" == "--no-serve" ]]; then
    SERVE=0
fi

fail() {
    echo "run_demo: $*" >&2
    exit 1
}

[[ -d "$REFS" ]] || fail "reference data not found at '$REFS' (set NEOFORGE_REFS or build it with fetch/fetch_manifests.py --out-dir '$REFS')"
for req in gencode hla mhci2018_mmc4.xlsx tesla/mmc4.xlsx; do
    [[ -e "$REFS/$req" ]] || fail "missing $REFS/$req"
done

cd "$ROOT"
python3 -c "import neoforge" 2>/dev/null || pip install -e "./pipeline" >/dev/null

CASES=(hcc1395 skcm_d9a3z1)
stage() { # stage <case> <stage>
    echo "==> [$1] $2"
    python3 -m neoforge.cli "$2" --case "$1" \
        --work-dir "$ROOT/work" --cases-dir "$ROOT/cases" --refs-dir "$REFS"
}

# skcm_d9a3z1 runs through translate only: without HLA-resolved data the
# present stage refuses, and downstream stages construct nothing on partial
# chains — the refusal is the integrity mechanism, not a failure.
for st in ingest variants translate present prioritize report construct; do
    stage hcc1395 "$st"
done
for st in ingest variants translate; do
    stage skcm_d9a3z1 "$st"
done
echo "==> [tesla] tesla_eval"
python3 -m neoforge.cli tesla_eval \
    --tesla-xlsx "$REFS/tesla/mmc4.xlsx" \
    --work-dir "$ROOT/work" --cases-dir "$ROOT/cases" --refs-dir "$REFS"

echo "==> export projections"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
python3 -m neoforge.export_ui \
    --work-dir "$ROOT/work" --cases-dir "$ROOT/cases" --out "$TMP"

echo "==> verifying exports match frontend/data, modulo per-run wall-clock"
for f in case_hcc1395.json case_skcm_d9a3z1.json benchmark_tesla.json; do
    [[ -f "$TMP/$f" ]] || fail "export did not produce $f"
    if ! python3 - "$TMP/$f" "$ROOT/frontend/data/$f" <<'PYEOF'
import json, sys

def load(path):
    def strip_volatile(obj):
        if isinstance(obj, dict):
            return {k: strip_volatile(v) for k, v in obj.items() if k != "wall_seconds"}
        if isinstance(obj, list):
            return [strip_volatile(v) for v in obj]
        return obj
    return strip_volatile(json.loads(open(path).read()))

a, b = load(sys.argv[1]), load(sys.argv[2])
sys.exit(0 if a == b else "exports differ beyond per-run wall-clock timings")
PYEOF
    then
        fail "$f differs from the tracked export (beyond wall-clock timings) — regenerate with neoforge.export_ui before serving"
    fi
    echo "    ok: $f"
done

if [[ "$SERVE" -eq 1 ]]; then
    echo "==> serving $ROOT/frontend at http://localhost:$PORT (Ctrl-C to stop)"
    exec python3 -m http.server "$PORT" --bind 127.0.0.1 -d "$ROOT/frontend"
else
    echo "run_demo: done (--no-serve)"
fi
