#!/usr/bin/env python3
"""Regenerate the mmacos manual's listings, figures and quoted numbers from the chapter drivers
(BRIEF_to_mmacos_manual.md: generated, never typed -- the MACOS manual's rule, Dave 2026-10-19).

Each chapter has a driver tools/transcripts/<name>.m: a MATLAB script run by `matlab -batch` in a
fresh scratch directory (it copies the template it narrates there and runs the copy, so a committed
template's own outputs are never touched).  The driver speaks through four helpers
(tools/transcripts/mmlib/), each a tagged stdout line:

  mm_transcript(NAME, text) / mm_listing(NAME, file, rows)
        -> src/transcripts/NAME.md, an indented code block; a chapter pulls it in with the marker
           line <!-- transcript: NAME --> (the MACOS manual's include_transcripts.lua)
  mm_figure(NAME, png)   -> src/media/NAME.png
  mm_pin(NAME, value, rule, shown)
        -> a number the prose quotes: checked against tools/pins/<driver>.tsv (rule 'sig5' = five
           significant figures on a signal, 'abs:<bound>' = a bound on a numerical zero,
           'printed' = a number a template prints at its own precision, compared by the listing rule -- the
           cross-architecture rule, so the Mac gates the manual the way it gates the engine), and
           written to src/transcripts/pins.tsv, where tools/pins.lua puts SHOWN in place of
           %pin:NAME% in the prose at build time

Every listing is gated too: on a rerun a listing may differ from the committed one only in numbers
that agree by the same rule (five significant figures, or the printed precision when fewer are
printed); text must match exactly.  A listing within the rule keeps its committed text.

A driver that errors, a pin that does not hold, or a listing that moved writes NOTHING (a listing cannot rot silently);
the run's full output is kept in build/logs/<name>.log.  A driver with no pin file yet records its
pins (first run); --repin rewrites them (a deliberate re-pin, to be reviewed in the diff).

Removed from what is printed (and only this): the home directory prints as `~`, the scratch
directory as `<work>`, the e2e thread's directory (build/thread, MMMAN_THREAD: the stages of the
worked thread consume each other's products, so their drivers share it) as `<thread>`.

Usage (from docs/mmacos-manual):
  python3 tools/regen_mmacos_transcripts.py [--repin] [driver names...]
Environment: MATLAB (default `matlab`), MMACOS_ROOT (default ~/dev/MACOS_resources/mmacos, the
mmacos the drivers add to the path), MACOS_HOME (default ~/dev/macos/macos_f90).
"""
import os, re, shutil, subprocess, sys, tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
MAN = HERE.parent
OUT = MAN / "src" / "transcripts"
MEDIA = MAN / "src" / "media"
DRIVERS = HERE / "transcripts"
LIB = DRIVERS / "mmlib"
PINS = HERE / "pins"
LOGS = MAN / "build" / "logs"
HOME = str(Path.home())
TAG = re.compile(r"^@@MM\|([TEFP])\|([\w.-]+)\|?(.*)$")


def sig5_equal(a, b):
    return f"{float(a):.4e}" == f"{float(b):.4e}"


def holds(rule, value, pinned):
    if rule == "sig5":
        return sig5_equal(value, pinned)
    if rule.startswith("abs:"):
        return abs(float(value)) < float(rule[4:])
    raise SystemExit(f"unknown pin rule {rule!r}")


NUM = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?")


def numbers_agree(a, b):
    """Two printed numbers agree by the cross-architecture rule: to five significant figures, or to
    the printed precision when fewer are printed (one unit in the last printed digit)."""
    fa, fb = float(a.replace("D", "e").replace("d", "e")), float(b.replace("D", "e").replace("d", "e"))
    if fa == fb:
        return True
    if not re.search(r"[.eEdD]", a):
        return False                     # an integer (a count, an index) matches exactly
    digits = len(re.sub(r"[eEdD].*", "", a).lstrip("+-").replace(".", "").lstrip("0")) or 1
    if digits >= 5:
        return sig5_equal(fa, fb)
    mant = re.sub(r"[eEdD].*", "", a)
    dec = len(mant.split(".")[1]) if "." in mant else 0
    exp = int(re.split(r"[eEdD]", a)[1]) if re.search(r"[eEdD]", a) else 0
    return abs(fa - fb) <= 10.0 ** (exp - dec) * 1.0000001


def transcript_drift(old, new):
    """The lines where NEW differs from the committed OLD beyond the rule (text must match exactly)."""
    bad = []
    ol, nl = old.splitlines(), new.splitlines()
    if len(ol) != len(nl):
        return [f"line count {len(ol)} -> {len(nl)}"]
    for i, (x, y) in enumerate(zip(ol, nl)):
        if x == y:
            continue
        if NUM.sub("#", x) != NUM.sub("#", y):
            bad.append(f"{i+1}: text differs: {x.strip()!r} -> {y.strip()!r}")
            continue
        for u, v in zip(NUM.findall(x), NUM.findall(y)):
            if not numbers_agree(u, v):
                bad.append(f"{i+1}: {u} -> {v}")
    return bad


def read_pins(path):
    pins = {}
    if path.is_file():
        for line in path.read_text().splitlines():
            if line.strip() and not line.startswith("#"):
                name, value, rule, shown = line.split("\t")
                pins[name] = (value, rule, shown)
    return pins


def write_pins(path, pins):
    rows = ["# name\tvalue\trule\tshown  (generated by regen_mmacos_transcripts.py; a change is a re-pin)"]
    rows += [f"{n}\t{v}\t{r}\t{s}" for n, (v, r, s) in pins.items()]
    path.write_text("\n".join(rows) + "\n")


def clean(line, work):
    return (line.replace(work, "<work>").replace(str(MAN / "build" / "thread"), "<thread>")
            .replace(HOME, "~"))


def run_driver(name, repin):
    drv = DRIVERS / f"{name}.m"
    work = tempfile.mkdtemp(prefix=f"mmman_{name}_")
    env = dict(os.environ)
    env.setdefault("MMACOS_ROOT", str(Path.home() / "dev" / "MACOS_resources" / "mmacos"))
    env.setdefault("MACOS_HOME", str(Path.home() / "dev" / "macos" / "macos_f90"))
    env["MMMAN_WORK"] = work
    # the e2e thread: each stage consumes the one before, so its drivers share one persistent
    # directory (build/thread); the first stage's driver (ch01) copies the template there afresh
    thread = MAN / "build" / "thread"
    thread.mkdir(parents=True, exist_ok=True)
    env["MMMAN_THREAD"] = str(thread)
    cmd = [env.get("MATLAB", "matlab"), "-batch",
           f"addpath('{LIB}'); cd('{work}'); run('{drv}');"]
    print(f"== {name}: running in {work}", flush=True)
    p = subprocess.run(cmd, cwd=work, env=env, capture_output=True, text=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    (LOGS / f"{name}.log").write_text(p.stdout + "\n--- stderr ---\n" + p.stderr)
    if p.returncode != 0:
        print(f"   FAILED (rc {p.returncode}); nothing written; build/logs/{name}.log", flush=True)
        return False

    texts, closed, figs, got = {}, set(), {}, {}
    for raw in p.stdout.splitlines():
        m = TAG.match(raw)
        if not m:
            continue
        kind, key, rest = m.groups()
        if kind == "T":
            texts.setdefault(key, []).append(clean(rest, work))
        elif kind == "E":
            texts.setdefault(key, [])
            closed.add(key)
        elif kind == "F":
            figs[key] = rest
        elif kind == "P":
            value, rule, shown = rest.split("|", 2)
            got[key] = (value, rule, shown)
    bad = [k for k in texts if k not in closed]
    if bad:
        print(f"   FAILED: transcript(s) not closed: {bad}; nothing written", flush=True)
        return False

    pinfile = PINS / f"{name}.tsv"
    pinned = read_pins(pinfile)
    if pinned and not repin:
        fails = []
        for k, (v, r, s) in got.items():
            if k not in pinned:
                fails.append(f"{k}: new pin (run with --repin to record it)")
            elif r == "printed":
                if not numbers_agree(pinned[k][2], s):
                    fails.append(f"{k}: printed {s} vs pinned {pinned[k][2]}")
            elif not holds(r, v, pinned[k][0]):
                fails.append(f"{k}: {v} vs pinned {pinned[k][0]} ({r})")
        fails += [f"{k}: pinned but not produced" for k in pinned if k not in got]
        if fails:
            print("   FAILED pins; nothing written:\n     " + "\n     ".join(fails), flush=True)
            return False
    if not pinned or repin:
        PINS.mkdir(exist_ok=True)
        write_pins(pinfile, got)
        print(f"   pins recorded: tools/pins/{name}.tsv ({len(got)})", flush=True)

    bodies = {k: "\n".join("    " + l if l.strip() else "" for l in lines).rstrip() + "\n"
              for k, lines in texts.items()}
    if not repin:
        drift = []
        for k, body in bodies.items():
            f = OUT / f"{k}.md"
            if f.is_file():
                drift += [f"{k}.md {d}" for d in transcript_drift(f.read_text(), body)]
        if drift:
            print("   FAILED: a listing moved beyond the cross-architecture rule; nothing written:\n     "
                  + "\n     ".join(drift[:40]), flush=True)
            return False
    OUT.mkdir(parents=True, exist_ok=True)
    MEDIA.mkdir(parents=True, exist_ok=True)
    for k, body in bodies.items():
        f = OUT / f"{k}.md"
        # a listing within the rule keeps its committed text (no churn between machines)
        if repin or not f.is_file() or transcript_drift(f.read_text(), body):
            f.write_text(body)
    for k, src in figs.items():
        shutil.copyfile(src, MEDIA / f"{k}.png")
    print(f"   wrote {len(texts)} transcripts, {len(figs)} figures, {len(got)} pins", flush=True)
    shutil.rmtree(work, ignore_errors=True)
    return True


def merge_pins():
    rows = []
    for f in sorted(PINS.glob("*.tsv")):
        for n, (v, r, s) in read_pins(f).items():
            rows.append(f"{n}\t{s}")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "pins.tsv").write_text("\n".join(rows) + "\n")


def main(argv):
    repin = "--repin" in argv
    names = [a for a in argv if not a.startswith("--")]
    if not names:
        names = sorted(p.stem for p in DRIVERS.glob("*.m"))
    ok = all([run_driver(n, repin) for n in names])
    merge_pins()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
