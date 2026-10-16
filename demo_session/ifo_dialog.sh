#!/usr/bin/env bash
# ifo_dialog.sh -- drive the polarization-PSI Twyman-Green demo ONE BEAT AT A
# TIME from outside MATLAB, so the demo can be a conversation instead of a
# script that runs away from you.
#
#   ./ifo_dialog.sh start [v1|v2]   launch a persistent MATLAB, run to the
#                                   point where beat 1 has finished and the
#                                   script is waiting
#   ./ifo_dialog.sh next            release ONE beat, wait for it, print it
#   ./ifo_dialog.sh show            reprint the last beat without advancing
#   ./ifo_dialog.sh log             the whole transcript so far
#   ./ifo_dialog.sh stop            shut the session down
#
# RUN IT PLAINLY -- do NOT pipe it into head/tail/grep.  The background
# MATLAB keeps the caller's pipe open, so the pipe never sees EOF and the
# command appears to hang while the session is in fact perfectly healthy.
# (setsid does not fix it, nor does detaching stdout; not worth more time --
# every subcommand already prints only what you want.)  Use `show` to
# reprint and `log` for the full transcript.
#
# HOW IT WORKS.  The demo script already pauses at every beat on an input()
# read from stdin.  We give MATLAB a FIFO as its stdin and push one newline
# per beat, so the thing being driven is the REAL demo -- the same file, the
# same beats, the same numbers Dave gets when he runs it by hand.  Nothing
# here is a demo-only reimplementation.
#
# Two details that are load-bearing:
#   * a FIFO delivers EOF as soon as the last writer closes, and MATLAB would
#     exit.  A holder process keeps the write end open for the whole session.
#   * MATLAB's stdout into a redirected file is BLOCK buffered, so "has the
#     beat finished" cannot be read from the log promptly.  The demo drops an
#     empty beat<N>.done file instead (TG_DEMO_MARKER) -- the filesystem is
#     not buffered.  Every wait below is bounded; none of them spins forever.

set -u
RES=/home/dcr/dev/MACOS_res_dev/mmacos
SESS=/tmp/ifo_dialog
FIFO=$SESS/cmd
LOG=$SESS/log
MAXWAIT=180          # seconds; beat 7 is the long one at ~20 s

case "${1:-}" in

start)
  # Children inherit this shell's stdout.  If the caller PIPES us
  # (`start | tail`), the background MATLAB and the FIFO holder keep that
  # pipe's write end open forever and the caller blocks -- the command looks
  # hung while the session is in fact perfectly healthy.  setsid does not fix
  # it; detaching the descriptor does.  Park the real stdout on fd 3, point
  # fd 1 at /dev/null so every child inherits THAT, and report on fd 3.
  exec 3>&1 1>/dev/null
  rig=${2:-v2}
  case "$rig" in
    v2) DIR=$RES/templates/90_polarization/tg_psi_dm_v2; SCRIPT=demo_tg_psi_v2.m ;;
    v1) DIR=$RES/templates/90_polarization/tg_psi_dm;    SCRIPT=demo_tg_psi.m ;;
    *)  echo "rig must be v1 or v2" >&3; exit 2 ;;
  esac
  "$0" stop >/dev/null 2>&1
  rm -rf "$SESS"; mkdir -p "$SESS"
  mkfifo "$FIFO"
  # holder: keeps the FIFO's write end open so MATLAB never sees EOF
  setsid sleep 100000 < /dev/null > "$FIFO" 2>/dev/null &
  echo $! > "$SESS/holder.pid"
  ( cd "$DIR" && \
    MACOS_HOME=/home/dcr/dev/macos/macos_f90 \
    TG_DEMO_PAUSE=1 TG_DEMO_MARKER="$SESS" \
    setsid matlab -nodisplay -nosplash < "$FIFO" > "$LOG" 2>&1 ) &
  echo $! > "$SESS/matlab.pid"
  echo "$rig" > "$SESS/rig"
  echo "$DIR" > "$SESS/figdir"
  if [ "$rig" = v2 ]; then echo demo2 > "$SESS/pfx"; else echo demo > "$SESS/pfx"; fi
  echo "0" > "$SESS/beat"
  # feed the two startup lines, then beat 1 runs unprompted (beat() only
  # pauses for n > 1) and the script blocks waiting for beat 2
  { echo "run('$RES/mmacos_setup.m')"; echo "run('$DIR/$SCRIPT')"; } > "$FIFO"
  # absorb MATLAB startup + beat 1 HERE, before anyone is watching, so the
  # first "run it" from the floor is instant instead of a 20-second boot
  for _ in $(seq 1 $((MAXWAIT*2))); do
    [ -f "$SESS/beat1.done" ] && { echo "$rig ready" >&3; exit 0; }
    sleep 0.5
  done
  { echo "!! $rig failed to start -- tail of the log:"; tail -25 "$LOG"; } >&3; exit 1
  ;;

next)
  n=$(( $(cat "$SESS/beat") + 1 ))
  [ "$n" -gt 1 ] && echo "" > "$FIFO"      # one Enter releases the next beat
  # wait for the marker the demo drops when beat n is finished and waiting
  for _ in $(seq 1 $((MAXWAIT*2))); do
    [ -f "$SESS/beat${n}.done" ] && break
    if ! kill -0 "$(cat "$SESS/matlab.pid")" 2>/dev/null; then
      echo "!! MATLAB exited -- tail of the log:" >&2; tail -25 "$LOG"; exit 1
    fi
    sleep 0.5
  done
  if [ ! -f "$SESS/beat${n}.done" ]; then
    echo "!! beat $n did not finish within ${MAXWAIT}s" >&2; tail -25 "$LOG"; exit 1
  fi
  echo "$n" > "$SESS/beat"
  # Four of the seven beats draw something.  A terminal cannot render an
  # image, so the picture goes on the DESKTOP -- one viewer window, opened
  # on cue, replacing the previous one so the screen never accumulates
  # clutter.  (Share the whole desktop, not just the terminal window.)
  case "$n" in
    1) fig=beat1_stack ;;   2) fig=beat2_layout ;;  3) fig=beat3_coating ;;
    4) fig=beat4_null ;;    5) fig=beat5_poke ;;
    6) fig=beat6_sweep ;;   7) fig=beat7_recovery ;;
    8) fig=beat8_calibration ;;
    *) fig= ;;
  esac
  if [ -n "$fig" ]; then
    f="$(cat "$SESS/figdir")/$(cat "$SESS/pfx")_$fig.png"
    g="${f%.png}.gif";  [ -f "$g" ] && f="$g"
    if [ -f "$f" ]; then
      pkill -x eog 2>/dev/null; sleep 0.3
      ( DISPLAY=:0 setsid eog "$f" < /dev/null > /dev/null 2>&1 & )
      echo "   [figure on screen: $(basename "$f")]"
    fi
  fi
  "$0" show
  ;;

show)
  n=$(cat "$SESS/beat")
  # the transcript between this beat's banner and its timing line, with the
  # engine's own trace chatter filtered out
  awk -v n="$n" '
    $0 ~ ("^#  BEAT " n " ") {on=1}
    on {print}
    $0 ~ ("\\[beat " n " took") {exit}
  ' "$LOG" \
  | sed 's/ *\[enter to run this beat\] *//' \
  | grep -vE '^ (Wavelength|u-v|x-y|Peak|Sum of|Tracing|OPD|RMS OPD|P-V OPD|Average OPD|Default|Compute)|^  (Elt|Number of|Geometric Prop|Polarization turned|Vector diff|Input file|Reading surface|mGridMat)|Grid surface array|Wavefront Prop|glass table|parameter file|Optical Train|Aperture grid|Tracing rays past|metrology beam lengths' \
  | grep -vE 'Graphics acceleration|performance might be diminished|mathworks.com|System Requirements|alternatePrintPath|^\[> In |^In [a-zA-Z]' \
  | cat -s
  ;;

ask)
  # Answer a question by COMPUTING it, in the session that just ran the demo.
  # After beat 7 the script has returned to the prompt with its whole
  # workspace intact, so a question about the result can be answered against
  # the actual data rather than from memory.  The body goes through a file to
  # dodge shell quoting entirely; diary captures the output and flushes on
  # close; a marker file signals completion (stdout to a file is buffered).
  rm -f "$SESS/ask.txt" "$SESS/ask.done"
  printf '%s\n' "$2" > "$SESS/ask_body.m"
  {
    echo "diary('$SESS/ask.txt');"
    echo "try, run('$SESS/ask_body.m'); catch e, fprintf('ERROR: %s\\n', e.message); end"
    echo "diary off;"
    echo "fid = fopen('$SESS/ask.done','w'); fclose(fid);"
  } > "$SESS/ask.m"
  echo "run('$SESS/ask.m')" > "$FIFO"
  for _ in $(seq 1 $((MAXWAIT*2))); do
    [ -f "$SESS/ask.done" ] && break
    sleep 0.5
  done
  [ -f "$SESS/ask.done" ] || { echo "!! no answer within ${MAXWAIT}s" >&2; exit 1; }
  grep -vE "^>> |^diary|ask_body|^$" "$SESS/ask.txt"
  ;;

log)  cat "$LOG" ;;

stop)
  for p in matlab holder; do
    [ -f "$SESS/$p.pid" ] && kill "$(cat "$SESS/$p.pid")" 2>/dev/null
  done
  # the matlab launcher forks the real process; take the session's children too
  pkill -x MATLAB 2>/dev/null
  pkill -x eog 2>/dev/null
  rm -rf "$SESS"
  echo "session stopped" >&2
  ;;

*) sed -n '2,25p' "$0"; exit 2 ;;
esac
