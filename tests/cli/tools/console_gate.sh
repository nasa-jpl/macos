#!/usr/bin/env bash
# console_gate.sh -- MACOS_CONSOLE at the compiler level: the engine's own routine
# (src/emsg_mod.F90, emsg_console_from_env) linked into a ten-line driver, run in every
# combination of the setting and of a console that is open or already closed.
#
# Why a gate below the bindings: the two compilers' runtimes differ here, and each
# difference was a crash or a stray file before it was handled (2026-11-20):
#   ifx       a failed console write ENDS the process (forrtl severe); WRITE(*) follows
#             unit 6 only with -assume noold_unit_star
#   gfortran  a failed console write is ignored; OPEN(UNIT=6) with the descriptor already
#             closed faulted (SIGSEGV) until unit 6 was closed first; an OPEN(UNIT=6) of a
#             file that cannot be opened left WRITE(*) going to a new file fort.6
#
# Usage:  tests/cli/tools/console_gate.sh [gfortran|ifx]      (default gfortran)
# Exit status 0 = every leg passed.  The pymacos tests (app/pymacos/tests/test_console.py)
# cover the same setting through the api's init.
FC="${1:-gfortran}"
here=$(cd "$(dirname "$0")" && pwd); root=$(cd "$here/../../.." && pwd)
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
cd "$TMP" || exit 2
cat > drv.f90 <<'F90'
program drv
  use emsg_mod, only: emsg_console_from_env
  integer :: k
  call emsg_console_from_env()
  call emsg_console_from_env()
  do k = 1, 3000
    write(*,*) 'a line to the console', k
    print *, 'and a print'
  end do
  write(0,'(a)') 'SURVIVED'
end program
F90
case "$FC" in
  ifx)      ifx -O1 -fpp -assume noold_unit_star -o drv "$root/src/emsg_mod.F90" drv.f90 > build.log 2>&1 ;;
  gfortran) gfortran -O1 -cpp -o drv "$root/src/emsg_mod.F90" drv.f90 > build.log 2>&1 ;;
  *) echo "console_gate: compiler must be gfortran or ifx"; exit 2 ;;
esac
[ -x drv ] || { echo "console_gate: the driver did not build with $FC"; tail -5 build.log; exit 2; }
fail=0
leg() {   # leg <setting or ''> <open|closed> <want: live|dead> <want stdout lines or -> <want file lines or ->
  local con="$1" cl="$2" want="$3" want_out="$4" want_file="$5" rc got out file stray
  rm -f engine.log fort.6 out.txt
  if [ "$cl" = closed ]; then ( [ -n "$con" ] && export MACOS_CONSOLE="$con" || unset MACOS_CONSOLE; ./drv >&- 2> err.txt ); rc=$?; out=-
  else ( [ -n "$con" ] && export MACOS_CONSOLE="$con" || unset MACOS_CONSOLE; ./drv > out.txt 2> err.txt ); rc=$?; out=$(wc -l < out.txt | tr -d ' '); fi
  if [ $rc -eq 0 ] && grep -q SURVIVED err.txt; then got=live; else got=dead; fi
  file=-; [ -f engine.log ] && file=$(wc -l < engine.log | tr -d ' ')
  stray=no; [ -f fort.6 ] && stray=yes
  if [ "$got" = "$want" ] && [ "$out" = "$want_out" ] && [ "$file" = "$want_file" ] && [ $stray = no ]; then
    echo "PASS  MACOS_CONSOLE=${con:-<unset>}, console $cl: $got, stdout lines $out, file lines $file"
  else
    echo "FAIL  MACOS_CONSOLE=${con:-<unset>}, console $cl: $got (want $want), stdout lines $out (want $want_out), file lines $file (want $want_file), fort.6 $stray; status $rc"
    fail=1
  fi
}
closed_on=dead; [ "$FC" = gfortran ] && closed_on=live       # the runtimes' own behavior with the console on
leg ""          open   live       6000 -
leg ""          closed $closed_on -    -
leg on          open   live       6000 -
leg on          closed $closed_on -    -
leg off         open   live       0    -
leg off         closed live       -    -
leg engine.log  open   live       0    6000
leg engine.log  closed live       -    6000
leg "$TMP/no_such_folder/engine.log" open   live       6001 -
leg "$TMP/no_such_folder/engine.log" closed $closed_on -    -
[ $fail -eq 0 ] && echo "console_gate ($FC): all legs pass" || echo "console_gate ($FC): FAILED"
exit $fail
