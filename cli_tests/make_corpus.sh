#!/usr/bin/env bash
# make_corpus.sh -- the deck lists for the CLI load gate (PLAN_CONSOLIDATION sec. 2.1).
#   corpus_committed.txt  every .in tracked by git in macos + MACOS_resources
#                         (ZGD_test_files, the manual's examples, mmacos/pymacos
#                         tests/Rx, templates, examples, challenges, GMI/test_ff)
#   corpus_extended.txt   the same plus every .in under $MACOS_SANDBOX (default
#                         ~/dev/MACOS_sandbox), untracked legacy decks: informational
# Absolute paths, sorted, one per line.  Re-run whenever decks are added.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
macos=$(cd "$here/.." && pwd)
res=${MACOS_RESOURCES:-$(cd "$macos/../MACOS_resources" && pwd)}
sandbox=${MACOS_SANDBOX:-$HOME/dev/MACOS_sandbox}
# a prescription declares elements; autoconf templates (giza/*.in) and the like do not
is_rx() { while IFS= read -r f; do grep -q -m1 -E '^[[:space:]]*(nElt|Element)[[:space:]]*=' "$f" 2>/dev/null && echo "$f"; done; }
{
  (cd "$macos" && git ls-files -- '*.in' | sed "s|^|$macos/|")
  (cd "$res"   && git ls-files -- '*.in' | sed "s|^|$res/|")
} | grep -v '/build\|/segmirmaker/test_in/' | sort -u | is_rx > "$here/corpus_committed.txt"
{
  cat "$here/corpus_committed.txt"
  [ -d "$sandbox" ] && find "$sandbox" -name '*.in' -type f 2>/dev/null | is_rx
} | sort -u > "$here/corpus_extended.txt"
# the per-commit core: the fixture sets and the manual, not the challenges' ladders
# (segmirmaker/test_in holds that tool's PARENT fragments -- inputs to segmirmaker,
#  not loadable prescriptions -- so they are left out of both lists)
grep -E '/ZGD_test_files/|/docs/macos-manual/|/tests/Rx/|/GMI/test_ff/' \
  "$here/corpus_committed.txt" > "$here/corpus_core.txt"
echo "core: $(wc -l < "$here/corpus_core.txt")  committed: $(wc -l < "$here/corpus_committed.txt")  extended: $(wc -l < "$here/corpus_extended.txt")"
