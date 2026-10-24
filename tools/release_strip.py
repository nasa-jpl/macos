#!/usr/bin/env python3
"""Apply a repo's release-exclude.txt to the index: the files it names leave `main`.

Usage (from the repo root, on the branch that will become main):
    python3 tools/release_strip.py            # dry run: prints what would go
    python3 tools/release_strip.py --apply    # git rm --cached -r the matches (then commit)

Pattern semantics (one per line, '#' comments): an exact tracked path; a directory
(every tracked file under it); a glob matched against the whole path, or -- when the
pattern has no '/' -- against the basename anywhere in the tree (so `BRIEF_*.md`
strips briefs in subdirectories too).  Idempotent: re-running after a merge from dev
strips only what came back.
"""
import fnmatch, os, subprocess, sys

def matches(path, pat):
    if path == pat or path.startswith(pat.rstrip('/') + '/'):
        return True
    if fnmatch.fnmatch(path, pat):
        return True
    return '/' not in pat and fnmatch.fnmatch(os.path.basename(path), pat)

def main():
    apply = '--apply' in sys.argv
    listfile = next((a for a in sys.argv[1:] if not a.startswith('--')), 'release-exclude.txt')
    pats = [l.strip() for l in open(listfile) if l.strip() and not l.lstrip().startswith('#')]
    files = subprocess.run(['git', 'ls-files', '-z'], capture_output=True, text=True, check=True).stdout.split('\0')
    files = [f for f in files if f]
    gone = sorted(f for f in files if any(matches(f, p) for p in pats))
    unused = [p for p in pats if not any(matches(f, p) for f in files)]
    print(f"{len(files)} tracked, {len(gone)} to strip by {len(pats)} patterns ({listfile})")
    for f in gone:
        print("  -", f)
    if unused:
        print("patterns matching nothing (fine if the file is already gone):", ", ".join(unused))
    if apply and gone:
        subprocess.run(['git', 'rm', '--cached', '-r', '-q', '--'] + gone, check=True)
        print("removed from the index; now: git commit -m 'Release strip (release-exclude.txt)'")
    elif apply:
        print("nothing to remove")

if __name__ == '__main__':
    main()
