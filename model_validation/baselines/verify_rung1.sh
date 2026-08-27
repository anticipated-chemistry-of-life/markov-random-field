#!/usr/bin/env bash
# Acceptance gate: the rung-1 pinned run must stay byte-identical.
#
#   ./verify_rung1.sh            rebuild, rerun, compare against the manifest
#   ./verify_rung1.sh --record   rebuild, rerun, and (re)write the manifest
#
# --record belongs to the one commit that deliberately moves the output, and to
# nothing else: during a behaviour-preserving stretch the whole point is that the
# manifest does not move. The manifest is a statement about a binary *and* a
# scenario -- the scenario is generated, so changing the harness that generates it
# invalidates it too. See README.md.
#
# acol.log is excluded: it carries a per-run ntfy topic UUID and wall-clock
# timings, neither of which say anything about the model.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$SCRIPT_DIR/../.."
OUT="$ROOT/model_validation/independent_y_s255_m255_seed42/rung1_gate"
MANIFEST="$SCRIPT_DIR/rung1_canonical_order.sha256"

"$SCRIPT_DIR/run_rung1.sh" >/dev/null

cd "$OUT"
# LC_ALL=C, not the ambient locale: collation differs between locales, so a manifest recorded
# under a UTF-8 locale and checked under LC_ALL=C (the default in most CI containers, and under
# `env -i`) diffs as FAIL with every hash identical -- the exact false alarm this gate exists to
# rule out. -print0/-0 for the same reason applied to whitespace in names.
hashes() { find . -type f ! -name 'acol.log' -print0 | LC_ALL=C sort -z | xargs -0 shasum -a 256; }

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
