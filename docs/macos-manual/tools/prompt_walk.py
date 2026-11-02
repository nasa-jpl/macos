#!/usr/bin/env python3
"""Walk a CLI dialog one prompt at a time and print PROMPT -> ANSWER pairs (an authoring aid for the
transcript journals in tools/transcripts/, chiefly the long NEW dialogs of Section 4).

  python3 tools/prompt_walk.py <answers file> [--binary B] [--model N] [--copy FILE ...] [--full]

The answers file is a journal's send lines: one answer per line, '#' lines ignored, $TMP replaced by a
fresh scratch directory (the CLI starts there, so NEW writes its deck there and OLD finds copies).
Each answer is sent only after the engine has printed a prompt and gone quiet, and the line it answers
is printed beside it, so an answer that lands on the wrong prompt shows at once.  At the end the
scratch directory's files are listed (the deck NEW wrote) and, with --full, the whole output printed.
Not a gate: the journal built from the walk is gated by cli_cmd_gate.py like every other.
"""
import argparse, os, re, select, shutil, sys, tempfile, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "cli_tests"))
from cli_load_gate import Cli  # noqa: E402

PROMPT = re.compile(r"(:|\?|\]|>)\s*$")


def quiet(cli, settle=0.25, tmax=60):
    """Read until the output ends in a prompt and nothing arrives for `settle` s."""
    end = time.time() + tmax
    last = time.time()
    while time.time() < end:
        r, _, _ = select.select([cli.fd], [], [], 0.05)
        if r:
            try:
                d = os.read(cli.fd, 65536)
            except OSError:
                return False
            if not d:
                return False
            cli.out += d
            last = time.time()
        elif time.time() - last > settle:
            tail = cli.out[-400:].decode(errors="replace").replace("\r", "")
            if PROMPT.search(tail.split("\n")[-1]):
                return True
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("answers")
    ap.add_argument("--binary", default=str(Path.home() / "dev/macos_cli_base/build_release/bin/macos"))
    ap.add_argument("--model", default="256")
    ap.add_argument("--copy", nargs="*", default=[])
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    tmp = tempfile.mkdtemp(prefix="pw_", dir="/tmp")
    for f in a.copy:
        shutil.copy(f, tmp)
    os.environ["CLI_PTY_COLS"] = "400"
    cwd = os.getcwd()
    os.chdir(tmp)
    try:
        cli = Cli(a.binary, a.model, 60)
        import fcntl, struct, termios
        fcntl.ioctl(cli.fd, termios.TIOCSWINSZ, struct.pack("HHHH", 50, 400, 0, 0))
        quiet(cli)
        for n, line in enumerate(open(os.path.join(cwd, a.answers)).read().splitlines(), 1):
            if line.startswith("#"):
                continue
            line = line.replace("$TMP", tmp)
            tail = cli.out.decode(errors="replace").replace("\r", "").rstrip("\n").split("\n")[-1]
            print(f"{n:4d}  {tail.strip()[-90:]:<90} <- {line!r}")
            cli.send(line)
            if not quiet(cli):
                print("  ** the CLI exited")
                break
        tail = cli.out.decode(errors="replace").replace("\r", "").rstrip("\n").split("\n")[-1]
        print(f"      {tail.strip()[-90:]}   (waiting here)")
        if a.full:
            print(cli.out.decode(errors="replace").replace("\r", ""))
        cli.close()
        print("files:", sorted(os.listdir(tmp)))
        for f in sorted(os.listdir(tmp)):
            if f.endswith(".in") and f not in [os.path.basename(c) for c in a.copy]:
                print(f"---- {f}\n" + open(f).read())
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
