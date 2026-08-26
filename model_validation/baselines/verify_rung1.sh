#!/usr/bin/env bash
# Phase-1 acceptance gate: the rung-1 pinned run must stay byte-identical.
#
#   ./verify_rung1.sh            rebuild, rerun, compare against the manifest
#   ./verify_rung1.sh --record   rebuild, rerun, and (re)write the manifest
#
# --record is for establishing a new baseline. During phase 1 (extract
# TPhylogeny) it must not be used: the whole point is that the manifest does not
# move. Phase 2 changes the chain deliberately and gets a statistical gate
# instead -- see the ADR on canonical node ordering.
#
# acol.log is excluded: it carries a per-run ntfy topic UUID and wall-clock
# timings, neither of which say anything about the model.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$SCRIPT_DIR/../.."
OUT="$ROOT/model_validation/independent_y_s255_m255_seed42/rung1_pre_split"
MANIFEST="$SCRIPT_DIR/rung1_pre_split.sha256"

"$SCRIPT_DIR/run_rung1.sh" >/dev/null

cd "$OUT"
hashes() { find . -type f ! -name 'acol.log' | sort | xargs shasum -a 256; }

if [[ "${1-}" == "--record" ]]; then
    hashes > "$MANIFEST"
    echo "recorded $(wc -l < "$MANIFEST" | tr -d ' ') files -> ${MANIFEST#"$ROOT"/}"
    exit 0
fi

[ -f "$MANIFEST" ] || { echo "error: no manifest at $MANIFEST; run with --record" >&2; exit 1; }

if diff <(hashes) "$MANIFEST" > /tmp/acol_gate_diff.$$; then
    echo "PASS: all $(wc -l < "$MANIFEST" | tr -d ' ') output files byte-identical to the baseline"
    rm -f /tmp/acol_gate_diff.$$
else
    echo "FAIL: output differs from the baseline" >&2
    cat /tmp/acol_gate_diff.$$ >&2
    rm -f /tmp/acol_gate_diff.$$
    exit 1
fi
