#!/usr/bin/env bash
# eof_gate.sh -- end of input at a prompt must stop a NO-READLINE macos, not loop.
#
# Windows builds never have readline, and neither does a build made with
# -DMACOS_FORCE_NO_READLINE=ON.  There PARSE_LOH reads the console with a bare
# READ; before 2026-11-20 end of input sent it straight back to the prompt, and
# a script that ran out at a sub-prompt printed the prompt without end (1.5 GB
# in three minutes on Windows).  Now it prints one line and stops with status 1,
# as the readline build and the model-size prompt already did.
#
# Usage:  tests/cli/tools/eof_gate.sh <path to a no-readline macos>
#   Build one with:
#     cmake -S . -B build_norl -G "Unix Makefiles" -DCMAKE_Fortran_COMPILER=gfortran \
#           -DCMAKE_C_COMPILER=gcc -DBUILD_MACOS=ON -DBUILD_SMACOS=OFF -DMACOS_FORCE_NO_READLINE=ON
#     cmake --build build_norl -j
# Exit status 0 = every leg passed.
BIN="${1:?usage: eof_gate.sh <no-readline macos binary>}"
[ -x "$BIN" ] || { echo "eof_gate: $BIN is not an executable"; exit 2; }
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
fail=0
# A 20 s limit on each leg, since the fault this gate looks for is a loop without end.
# macOS has no `timeout` command (CCMac, 2026-11-20: every leg returned status 127): use
# GNU timeout where it exists, Homebrew's gtimeout next, and perl's alarm otherwise.
# EOF_GATE_LIMIT=perl forces the last one (to test it where timeout exists).
limited() {
  if [ "${EOF_GATE_LIMIT:-}" != perl ] && command -v timeout > /dev/null 2>&1; then timeout 20 "$@"
  elif [ "${EOF_GATE_LIMIT:-}" != perl ] && command -v gtimeout > /dev/null 2>&1; then gtimeout 20 "$@"
  else perl -e 'alarm shift; exec @ARGV or exit 127' 20 "$@"
  fi
}
leg() {   # leg <name> <input> <want status> <want line count of the stop message> <max bytes>
  local name="$1" input="$2" want_rc="$3" want_n="$4" max="$5" rc n size
  printf "$input" | limited "$BIN" 2>&1 | head -c 2000000 > "$TMP/out.txt"
  rc=${PIPESTATUS[1]}
  n=$(grep -c 'no answer (end of input); stopping' "$TMP/out.txt")
  size=$(wc -c < "$TMP/out.txt")
  if [ "$rc" = "$want_rc" ] && [ "$n" = "$want_n" ] && [ "$size" -le "$max" ]; then
    echo "PASS  $name: status $rc, $size bytes, stop line x$n"
  else
    echo "FAIL  $name: status $rc (want $want_rc), $size bytes (limit $max), stop line x$n (want $want_n)"
    fail=1
  fi
}
leg "end of input at the STOP sub-prompt" '128\nstop\n' 1 1 20000
leg "end of input at the main prompt"     '128\n'       1 1 20000
leg "a script that ends with QUIT"        '128\nquit\n' 0 0 20000
exit $fail
