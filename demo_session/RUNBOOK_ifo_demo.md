# RUNBOOK — the live IFO demo (polarization-PSI Twyman-Green gauging a DM)

For 2026-09-01.  Two rigs: **v2** (`tg_psi_dm_v2`, a real MacNeille PBS cube,
eight beats) and **v1** (`tg_psi_dm`, the older plate rig, seven beats).  v2
is the default.  Every timing below was measured on this box.

---

## 0.  THE WAKE-UP LINE

Dave types, at the start of a **fresh** Claude Code session:

> **Let's run the IFO demo.  Runbook is `demo_session/RUNBOOK_ifo_demo.md` —
> arm it and wait for my cue.**

---

## FOR CC — read this section first, then do exactly this

You have no prior context.  Do not improvise; the room is watching the
conversation window.

1. **Arm the session, silently:**
   `cd /home/dcr/dev/macos/demo_session && ./ifo_dialog.sh start v2`
   Takes ~25 s (MATLAB boot + beat 1, absorbed before anyone is watching).
   It prints `v2 ready`.
2. **Reply in ONE line.**  Something like: *"Armed — v2, eight beats, beat 1
   computed and waiting.  Say the word."*  Nothing longer.  Your replies are
   projected; a wall of setup text is the failure mode here.
3. **Then follow the dialog protocol below and nothing else.**  Do not run
   ahead.  Do not run a beat until Dave says so.

### The protocol

| Dave says | you do |
|---|---|
| *"tell us about the first step"* / *"tell us about the next step"* | **Three or four lines**, plain language, no jargon. Say what the step does and what will appear. Then stop — do NOT run it. |
| *"run it"* | `./ifo_dialog.sh next`. **Paste the real output verbatim** in a fenced block (trim with a visible `…` if it is long; never reformat it). Then three or four short paragraphs pointing at the one thing that matters. |
| a question | see the three-bucket contract below |
| *"back up"* / *"wait"* | stop, answer, do not advance |

### Two rules that are load-bearing

- **A terminal renders markdown, not images.**  `SendUserFile` shows nothing.
  Figures go on the DESKTOP — the driver opens them automatically for beats
  2, 4, 5, 6, 7 and 8, one viewer window at a time.  Dave shares the whole
  desktop, not just the terminal.
- **Your reply IS the screen.**  Short. The output block plus a few
  paragraphs. No preamble, no recap of previous beats.

### The three-bucket contract for questions

Say which bucket, immediately, before doing anything:

- **answer now** — it is on screen or already measured; just answer.
- **one run** — say roughly how long, then compute it.  The MATLAB session
  stays alive at the prompt after beat 8 with the entire workspace loaded, so
  use it: `./ifo_dialog.sh ask '<matlab code>'`.  A real question gets a real
  number in well under a minute.  (Worked example: the calibration question
  that became beat 8.)
- **that's real work** — say so, note it, move on.  Do not start it live.

### If a beat stalls

Every drawing beat writes its PNG *before* it finishes.  Show the file and
keep talking; do not debug in front of the room.  `./ifo_dialog.sh show`
reprints the last beat without advancing.

---

## 1.  The eight beats

| # | beat | measured | draws? | the line that lands |
|---|---|---|---|---|
| 1 | BUILD | 0.3 s | — | The design condition *hands back the glass* — n = 1.6555, a dense flint. That is why real polarizing cubes are not BK7. |
| 2 | LAYOUT | 8.0 s | ✓ | Both returns leave by the **same** port. That is what a PBS buys, and it is 2.3× the light of a plate splitter. |
| 3 | COATING | 2.9 s | — | `R+T = 1.000000000000` measured across **two separate arm models**. Then: same textbook rule, one layer of termination different, and R_p goes from 4e-12 to 2.1%. |
| 4 | NULL | 3.6 s | ✓ | 73 picometres, with **nothing aligned**. And the leftover has a smooth saddle shape — a systematic you could chase, not a noise floor. |
| 5 | POKE | 2.5 s | ✓ | The bullseye does the arithmetic: half a fringe, double-passed, is ~158 nm. 146 recovered against 150 injected — and the demo *prints* the gap. |
| 6 | SWEEP | 2.4 s | ✓ (GIF) | 36 frames in 35 ms and **zero ray traces**. Nothing in the interferometer moved. |
| 7 | RECOVER | 13.0 s | ✓ | 6.26 nm measured against 6.35 truth, 0.183 nm residual. The truth panel is the thing you cannot get on a bench. |
| 8 | CALIBRATE | 12.0 s | ✓ | Gain 0.9912 ± 0.0020, linear to 0.00%. There is almost nothing to calibrate — and the two correlations say *where* the leftover error lives. |

Roughly 45 s of compute across the whole run; the rest is talking.

**Beat 8 is the "a customer asks a hard question" beat.**  Its conclusion is
deliberately *not* the one it was built expecting: over low-order aberrations
the gauge is already at its floor, so the answer is widen the basis, not the
gain.

---

## 2.  Cold-start checks (do these before the room fills)

Open a **terminal** (not the GNOME menu — MATLAB must inherit `MACOS_HOME`):

```bash
echo $MACOS_HOME          # /home/dcr/dev/macos/macos_f90
cd /home/dcr/dev/macos/demo_session
./ifo_dialog.sh start v2  # ~25 s, prints "v2 ready"
./ifo_dialog.sh next      # beat 1 should print instantly
./ifo_dialog.sh stop
```

If `start` fails it prints the tail of the MATLAB log and exits non-zero.

---

## 3.  The driver

```
./ifo_dialog.sh start [v1|v2]   boot MATLAB, run beat 1, hold
./ifo_dialog.sh next            release one beat, wait, print it
./ifo_dialog.sh show            reprint the last beat, no advance
./ifo_dialog.sh ask '<code>'    compute in the live session (see above)
./ifo_dialog.sh log             the whole transcript
./ifo_dialog.sh stop            tear it all down
```

**Run these plainly — never pipe them into `head`/`tail`/`grep`.**  The
background MATLAB holds the caller's pipe open, so it never sees EOF and
the command looks hung while the session is fine.  Every subcommand
already prints only what you want; use `show` to reprint a beat and `log`
for the whole transcript.

It drives the **real** demo script — same file Dave runs by hand, same beats,
same numbers.  It feeds one newline per beat into MATLAB's stdin through a
FIFO.  Two details it depends on: a holder process keeps the FIFO's write end
open (otherwise MATLAB sees EOF and exits), and each beat drops a
`beat<N>.done` file because stdout into a redirected file is block-buffered
while the filesystem is not.

Leftovers while a session is up: MATLAB (~300 MB), a `sleep 100000` holder,
and one `eog`.  `stop` takes all three.

---

## 4.  Fallback — no CC, drive it by hand

```matlab
>> run('/home/dcr/dev/MACOS_res_dev/mmacos/mmacos_setup.m')
>> cd /home/dcr/dev/MACOS_res_dev/mmacos/templates/90_polarization/tg_psi_dm_v2
>> demo_tg_psi_v2                 % Enter between beats; Ctrl-C aborts
```

Launch the **full desktop** (`matlab`, no `-batch`, no `-nodesktop`) or the
per-beat pauses do not happen — `usejava('desktop')` is what gates them, and
it returns 1 here.  `cd` matters: the decks reference their grid file by a
relative name.  A cold directory is fine; beat 1 writes what it needs.

v1 is the same, in `../tg_psi_dm`, as `demo_tg_psi` (seven beats, no
calibration beat).  Both demos are model size 256, so v1 and v2 can be run
back to back in one MATLAB session.  Do **not** run the test suite in that
session — `tTgPol2` is model 128, and it is one model size per process.

**Total-failure fallback:** the committed figures carry the result with no
live run at all — `demo2_beat*.png` in the v2 directory, and
`tg_psi_dm_v2_recovery.png` / `_sensitivity.png`.

---

## 5.  "How do you know it's right"

```matlab
% FRESH session -- model 128, do not mix with the demos
>> runtests('tests/tTgPol2.m')     % 9 gates, ~6 s   (v2)
>> runtests('tests/tTgPol.m')      % 9 gates, ~5 s   (v1)
```

The one to narrate is `test_engine_coated_diagonal_matches_the_macleod_analytic`:
the engine's coated 45° stack against Macleod's characteristic matrix,
written from the textbook and never transcribed from the Fortran — `R_s` to
9.6e-10.  Two gates in that file are deliberate counterexamples, so the suite
cannot pass vacuously.

---

## 6.  Where things live

| | |
|---|---|
| v2 rig, demo, README (every number) | `MACOS_res_dev/mmacos/templates/90_polarization/tg_psi_dm_v2/` |
| v1 rig | `…/tg_psi_dm/` |
| gates | `MACOS_res_dev/mmacos/tests/tTgPol.m`, `tTgPol2.m` |
| briefs + delivery logs | `macos/BRIEF_tg_demo.md`, `macos/BRIEF_tg_ifo_v2.md` |
| builder | `+macos/+design/twyman_green.m`, `pbs_macneille.m`, `thinfilm_rt.m`, `Bench.m` |

All LOCAL on `dev-candidate` in `MACOS_res_dev`.  Nothing pushed.
