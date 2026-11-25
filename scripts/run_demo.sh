#!/usr/bin/env bash
# Thin POSIX wrapper: all behaviour now lives in scripts/run_demo.py so the
# demo runs identically under PowerShell, bash and plain `python`.
#   ./scripts/run_demo.sh [--no-serve] [--fetch-refs] [--refs-dir DIR]
set -euo pipefail
exec python3 "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_demo.py" "$@"
