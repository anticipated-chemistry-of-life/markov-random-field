#!/usr/bin/env bash
# Acceptance gate: the pinned rung-1 and rung-2 runs must stay byte-identical.
#
#   ./verify.sh            rebuild, rerun both rungs, compare against the manifests
#   ./verify.sh --record   rebuild, rerun both rungs, and (re)write the manifests
#
# --record belongs to the one commit that deliberately moves the output, and to
# nothing else: during a behaviour-preserving stretch the whole point is that the
# manifests do not move. A manifest is a statement about a binary *and* a
# scenario -- the scenario is generated, so changing the harness that generates it
# invalidates it too. See README.md.
#
# acol.log is excluded: it carries a per-run ntfy topic UUID and wall-clock
# timings, neither of which say anything about the model.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$SCRIPT_DIR/../.."
SCENARIO="$ROOT/model_validation/independent_y_s255_m255_seed42"

# LC_ALL=C, not the ambient locale: collation differs between locales, so a manifest recorded
# under a UTF-8 locale and checked under LC_ALL=C (the default in most CI containers, and under
# `env -i`) diffs as FAIL with every hash identical -- the exact false alarm this gate exists to
# rule out. -print0/-0 for the same reason applied to whitespace in names.
hashes() { find . -type f ! -name 'acol.log' -print0 | LC_ALL=C sort -z | xargs -0 shasum -a 256; }

status=0
for rung in rung1 rung2; do
    manifest="$SCRIPT_DIR/${rung}_canonical_order.sha256"
    "$SCRIPT_DIR/run_rung.sh" "$rung" >/dev/null
    cd "$SCENARIO/${rung}_gate"

    if [[ "${1-}" == "--record" ]]; then
        hashes > "$manifest"
        echo "$rung: recorded $(wc -l < "$manifest" | tr -d ' ') files -> ${manifest#"$ROOT"/}"
        continue
    fi

    [ -f "$manifest" ] || { echo "error: no manifest at $manifest; run with --record" >&2; exit 1; }
    diff_file="$(mktemp)"
    if diff <(hashes) "$manifest" > "$diff_file"; then
        echo "$rung PASS: all $(wc -l < "$manifest" | tr -d ' ') output files byte-identical"
    else
        echo "$rung FAIL: output differs from the baseline" >&2
        cat "$diff_file" >&2
        status=1
    fi
    rm -f "$diff_file"
done
exit "$status"
