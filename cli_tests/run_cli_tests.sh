#!/usr/bin/env bash
# run_cli_tests.sh -- the CLI regression suite (PLAN_CONSOLIDATION sec. 2).
#
#   ./run_cli_tests.sh load            the load gate on the CORE corpus, both compilers
#   ./run_cli_tests.sh corpus          the load gate on EVERY committed deck (weekly / A/B)
#   ./run_cli_tests.sh extended        ... plus the untracked sandbox decks (informational)
#   ./run_cli_tests.sh cmd             (sec. 2.2, not yet) one journal per cmdref entry
#   ./run_cli_tests.sh ab <pre> <post> [kind] [allow]   (sec. 2.3) the load-gate RECORDS at two
#                                      SHAs compared deck by deck, both compilers; kind = load
#                                      (default) | corpus; allow = a file naming the decks the
#                                      change was meant to move (an ab_*_moved record serves).
#                                      Exit 1 on an unnamed mover or an ok -> not-ok deck.  The
#                                      records must exist: run `load`/`corpus` from the clean
#                                      worktree at each SHA first (see README).
#   ./run_cli_tests.sh full            load + cmd (ab needs its two SHAs)
#
# Records: cli_tests/records/<kind>_<compiler>_<sha>[+dirty].csv -- the record IS the
# result; this script prints the summary, writes <kind>_<compiler>_latest.txt, and
# exits 1 on (a) a must-fail deck that loads, (b) any crash_*/timeout_*/load_fail on
# a deck the previous record of the same kind had as ok, (c) any round trip that was
# 'same' before and is not now.  New failures are listed; pre-existing ones are
# counted and named in the summary, not fatal (PLAN sec. 0 items are closed one at a
# time; the gate's job is to stop NEW ones).
#
# Binaries: by default from a CLEAN worktree of the engine at a committed SHA
# (MACOS_CLI_TREE, default ~/dev/macos_cli_base: `git worktree add ../macos_cli_base
# <sha>` then makems.sh release [gfortran] there) -- never the shared working tree,
# whose build dirs another lane may be rebuilding mid-edit.  MACOS_CLI_TREE=. uses
# this tree (the record is then tagged +dirty if the engine source is modified).
# MACOS_CLI_ONLY=ifx|gfortran limits the compilers; a missing binary is skipped.
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)
tree=${MACOS_CLI_TREE:-$HOME/dev/macos_cli_base}
[ "$tree" = "." ] && tree="$root"
tree=$(cd "$tree" && pwd)
kind=${1:-load}
mkdir -p "$here/records"
sha=$(cd "$tree" && git rev-parse --short HEAD)
(cd "$tree" && git diff --quiet -- macos_f90 CMakeLists.txt) || sha="${sha}+dirty"
echo "binaries from $tree at $sha"

case "$kind" in
  load)     list="$here/corpus_core.txt" ;;
  corpus)   list="$here/corpus_committed.txt" ;;
  extended) list="$here/corpus_extended.txt" ;;
  cmd)
    # sec. 2.2: one pty journal per command (cmd/*.jou), must_fail/ journals must NOT pass
    for comp in ifx gfortran; do
      [ -n "${MACOS_CLI_ONLY:-}" ] && [ "$MACOS_CLI_ONLY" != "$comp" ] && continue
      case $comp in ifx) bin="$tree/build_release/bin/macos" ;; gfortran) bin="$tree/build_release_gfortran/bin/macos" ;; esac
      [ -x "$bin" ] || { echo "== $comp: no binary (skipped)"; continue; }
      rec="$here/records/cmd_${comp}_${sha}.csv"; echo "== $comp: journals -> $rec"
      python3 "$here/cli_cmd_gate.py" "$bin" "$here/cmd" "$rec" | tee "$here/records/cmd_${comp}_latest.txt"
      awk -F, 'NR>1 && $2=="journal" && $3!="pass" {bad=1} NR>1 && $2=="must_fail" && $3=="pass" {bad=1} END {exit bad}' "$rec" || fail=1
    done
    exit ${fail:-0} ;;
  manual)
    # Dave 2026-10-09: load and re-emit the manual's example decks for documentation --
    # the load gate's FIRST SAVE of each, kept as docs/macos-manual/examples/emitted/<name>.in
    list="$here/corpus_manual.txt"; ls "$root"/docs/macos-manual/examples/*.in > "$list"
    bin="$tree/build_release_gfortran/bin/macos"; [ -x "$bin" ] || bin="$tree/build_release/bin/macos"
    wd="$here/records/manual_work"; rm -rf "$wd"; mkdir -p "$wd" "$root/docs/macos-manual/examples/emitted"
    python3 "$here/cli_load_gate.py" "$bin" "$list" "$here/records/manual_${sha}.csv" --workdir "$wd" | tee "$here/records/manual_latest.txt"
    i=0; while read -r d; do n=$(basename "$d" .in); f=$(printf '%s/rt%04d_a.in' "$wd" $i)
      [ -f "$f" ] && cp "$f" "$root/docs/macos-manual/examples/emitted/$n.in"; i=$((i+1)); done < "$list"
    rm -rf "$wd"; echo "emitted: $(ls "$root"/docs/macos-manual/examples/emitted | wc -l) decks"; exit 0 ;;
  ab)
    pre=${2:?pre sha}; post=${3:?post sha}; akind=${4:-load}; allow=${5:-}
    for comp in ifx gfortran; do
      [ -n "${MACOS_CLI_ONLY:-}" ] && [ "$MACOS_CLI_ONLY" != "$comp" ] && continue
      a="$here/records/${akind}_${comp}_${pre}.csv"; b="$here/records/${akind}_${comp}_${post}.csv"
      for f in "$a" "$b"; do [ -f "$f" ] || { echo "== $comp: no record $f -- run '$0 $akind' from the clean worktree at that SHA"; fail=1; continue 2; }; done
      echo "== $comp: $pre -> $post ($akind)"
      python3 "$here/ab_records.py" "$a" "$b" ${allow:+--allow "$allow"} --out "$here/records/ab_${akind}_${comp}_${pre}_${post}.txt" | tail -n +3 | tee "$here/records/ab_${comp}_latest.txt"
      [ "${PIPESTATUS[0]}" -eq 0 ] || fail=1
    done
    exit ${fail:-0} ;;
  full) "$0" load; r1=$?; "$0" cmd; r2=$?; exit $(( r1 || r2 )) ;;
  *) echo "usage: $0 [load|corpus|extended|cmd|manual|ab <pre> <post> [kind] [allow]|full]"; exit 2 ;;
esac
[ -s "${list:-}" ] || "$here/make_corpus.sh"
ls "$here"/must_fail/*.in > "$here/must_fail/list.txt"      # absolute paths, regenerated each run

fail=0
for comp in ifx gfortran; do
  [ -n "${MACOS_CLI_ONLY:-}" ] && [ "$MACOS_CLI_ONLY" != "$comp" ] && continue
  case $comp in ifx) bin="$tree/build_release/bin/macos" ;; gfortran) bin="$tree/build_release_gfortran/bin/macos" ;; esac
  if [ ! -x "$bin" ]; then echo "== $comp: no binary at $bin (skipped)"; continue; fi
  rec="$here/records/${kind}_${comp}_${sha}.csv"
  mf="$here/records/mustfail_${comp}_${sha}.csv"
  echo "== $comp: $(wc -l < "$list") decks -> $rec"
  python3 "$here/cli_load_gate.py" "$bin" "$here/must_fail/list.txt" "$mf" --deadline 30 > /dev/null
  python3 "$here/cli_load_gate.py" "$bin" "$list" "$rec" > "${rec%.csv}.log"
  prev=$(ls -t "$here/records/${kind}_${comp}_"*.csv 2>/dev/null | grep -v "$rec" | head -1)
  python3 "$here/judge_load.py" "$rec" "$mf" "${prev:-}" | tee "$here/records/${kind}_${comp}_latest.txt"
  [ "${PIPESTATUS[0]}" -eq 0 ] || fail=1
done
# the two compilers against each other (informational): a deck whose pass count or
# lost count differs, or whose RMS OPD differs by more than 1e-9 relative, is named
a="$here/records/${kind}_ifx_${sha}.csv"; b="$here/records/${kind}_gfortran_${sha}.csv"
if [ -f "$a" ] && [ -f "$b" ]; then python3 "$here/compare_compilers.py" "$a" "$b" | tee -a "$here/records/${kind}_compilers_latest.txt"; fi
exit $fail
