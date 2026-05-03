"""
fix_examples.py — Normalize the extracted .in / .jou files so MACOS
can read them.

The PDF -> DOCX conversion gave us inconsistent line breaking inside
prescription bodies:

  - Some files run multiple `Key= value` pairs together on one
    physical line (Cassegrain.in style).
  - Other files put the key on one line and the values on the next
    (AOExample.in style).

MACOS's reader wants each prescription key + its values on a single
line.  This script normalizes every file in examples/ to that form by:

  1. Joining every line to its successors until it hits the next
     line that begins with a key (`\b\w+=`) -- so values on
     continuation lines flow back up.
  2. Splitting each resulting line at every NEW `\b\w+=` (after the
     first), so multi-key concatenated lines get one key per line.

Comments (lines starting with %) are preserved unchanged.
Empty lines are dropped.
"""
import re
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
EXAMPLES_DIR = PROJECT_ROOT / "examples"

# A "key" is an identifier followed by '='.  Identifiers are letters
# digits and underscore, must start with a letter.
KEY_RE = re.compile(r"\b([A-Za-z][\w]*)=")

# "Continuation" pattern: a key at the start of a line, with no value
# after it.  Such lines must absorb the next non-empty line.
CONT_RE = re.compile(r"^\s*[A-Za-z][\w]*=\s*$")


def is_comment(line: str) -> bool:
    return line.lstrip().startswith("%")


def is_blank(line: str) -> bool:
    return not line.strip()


def fix_one(text: str) -> str:
    raw_lines = text.splitlines()

    # Pass 1: join continuation lines.  If a line ends with `Key=` and
    # nothing else, pull the next non-empty non-comment line into it.
    joined = []
    i = 0
    while i < len(raw_lines):
        line = raw_lines[i]
        if is_comment(line) or is_blank(line):
            joined.append(line)
            i += 1
            continue
        # Greedily absorb successors that don't start with a fresh key
        # AND are not comments / blanks, up until either:
        #   - the line we just built contains text after at least one '='
        #     (i.e. has a value);
        #   - the next line *itself* starts with a `\w+=` (new key).
        while True:
            j = i + 1
            if j >= len(raw_lines):
                break
            nxt = raw_lines[j]
            if is_blank(nxt) or is_comment(nxt):
                break
            # Does the current line have a value after its last '=' ?
            tail = line.rsplit("=", 1)[-1].strip() if "=" in line else line.strip()
            current_has_value = bool(tail) and "=" not in tail
            # Does next start with a new key?
            nxt_starts_key = bool(KEY_RE.match(nxt.lstrip()))
            if current_has_value and nxt_starts_key:
                # Two distinct entries; stop joining.
                break
            if nxt_starts_key and not current_has_value:
                # Current line is `Key=` only AND next is also a key:
                # the current key has no value at all -- emit it standalone.
                break
            # Join.
            line = line.rstrip() + " " + nxt.lstrip()
            i += 1
        joined.append(line)
        i += 1

    # Pass 2: within each non-comment line, split at every fresh `Key=`
    # other than the first.
    out = []
    for line in joined:
        if is_comment(line) or is_blank(line):
            out.append(line)
            continue
        # Find positions of every key-equals in the line.
        positions = [m.start() for m in KEY_RE.finditer(line)]
        if len(positions) <= 1:
            out.append(line.rstrip())
            continue
        # Slice between successive key positions.
        chunks = []
        for k, pos in enumerate(positions):
            end = positions[k + 1] if k + 1 < len(positions) else len(line)
            chunks.append(line[pos:end].strip())
        # Preserve any pre-key prefix (e.g. leading whitespace before
        # the first key) by attaching it to the first chunk.
        prefix = line[:positions[0]].rstrip()
        if prefix:
            chunks[0] = (prefix + " " + chunks[0]).strip()
        out.extend(chunks)

    # Drop runs of consecutive blank lines (collapse to one).
    cleaned = []
    prev_blank = False
    for line in out:
        if is_blank(line):
            if not prev_blank:
                cleaned.append("")
            prev_blank = True
        else:
            cleaned.append(line.rstrip())
            prev_blank = False

    return "\n".join(cleaned).rstrip() + "\n"


def main():
    # Only the .in files use Key= value semantics and benefit from the
    # join+split normalization.  .jou files are one MACOS command per
    # line, and .test files are pure numeric data -- the Key= heuristic
    # would mangle both.
    files = sorted(EXAMPLES_DIR.glob("*.in"))
    if not files:
        sys.exit(f"No .in example files found under {EXAMPLES_DIR}/")

    for f in files:
        before = f.read_text()
        after = fix_one(before)
        if before == after:
            print(f"  unchanged  {f.name}")
            continue
        f.write_text(after)
        nb = before.count("\n")
        na = after.count("\n")
        print(f"  fixed      {f.name}   ({nb} -> {na} lines)")

    # Lightweight pass on .jou: only strip running-page-footer lines
    # (everything else is one command per line as-is).
    for f in sorted(EXAMPLES_DIR.glob("*.jou")):
        before = f.read_text()
        cleaned = []
        for line in before.splitlines():
            t = line.strip()
            if not t:
                cleaned.append("")
                continue
            # Drop "Modeling and Analysis... NN" footer lines
            if t.startswith("Modeling and Analysis"):
                continue
            cleaned.append(line.rstrip())
        # Collapse multiple blank lines.
        out = []
        prev_blank = False
        for line in cleaned:
            if not line.strip():
                if not prev_blank:
                    out.append("")
                prev_blank = True
            else:
                out.append(line)
                prev_blank = False
        after = "\n".join(out).rstrip() + "\n"
        if before != after:
            f.write_text(after)
            print(f"  trimmed    {f.name}   "
                  f"({before.count(chr(10))} -> {after.count(chr(10))} lines)")
        else:
            print(f"  unchanged  {f.name}")

    # .test files are left strictly alone.


if __name__ == "__main__":
    main()
