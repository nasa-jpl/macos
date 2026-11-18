#!/usr/bin/env python3
"""Name the previous record of a kind: prev_record.py <records dir> <kind> <compiler> <this record.csv>
Prints the path of the most recently MADE <kind>_<compiler>_*.csv other than this one, or nothing.

A record's age is the date of the commit that added it; a record not yet committed is
aged by its file time.  File times alone cannot order the committed records: a clone
or a checkout stamps every file with the same moment, and the pick among them was then
arbitrary (2026-11-16, the first load gate on the re-cloned repository judged ifx
against a record five weeks old and named 16 movers that were not).  Outside a git
repository every record falls back to its file time.

The records directory has moved (cli_tests/records became tests/cli/records, 2026-11), and a
move makes git report every record as added by the commit that moved it.  So the history is
read for the record's NAME in any directory, with renames followed (-M): a rename is then
not an addition, and the date is still that of the commit that first made the record."""
import glob, os, subprocess, sys


def git(d, *a):
    try:
        r = subprocess.run(('git', '-C', d) + a, capture_output=True, text=True)
        return r.stdout.splitlines() if r.returncode == 0 else []
    except OSError:
        return []


def main(d, kind, comp, this):
    pat = '%s_%s_*.csv' % (kind, comp)
    this = os.path.realpath(this)
    cands = [p for p in glob.glob(os.path.join(d, pat)) if os.path.realpath(p) != this]
    if not cands:
        return ''
    tracked = set(git(d, 'ls-files', '--', pat))
    added, t = {}, 0
    for ln in git(d, 'log', '-M', '--diff-filter=A', '--format=@%ct', '--name-only', '--',
                  ':(top,glob)**/' + pat):
        if ln.startswith('@'):
            t = int(ln[1:])
        elif ln:
            n = os.path.basename(ln)                    # the record's name, whatever directory it was in then
            added[n] = max(t, added.get(n, 0))          # deleted and re-added: the later date

    def age(p):
        n = os.path.basename(p)
        if n not in tracked:
            return os.path.getmtime(p)                  # not committed yet: made on this machine
        if n not in added:                              # e.g. it arrived in a merge's own tree
            ln = git(d, 'log', '-1', '-M', '--follow', '--diff-filter=AM', '--format=%ct', '--', n)
            return int(ln[0]) if ln else 0              # its last change of content, not a move
        return added[n]

    return max(cands, key=lambda p: (age(p), os.path.basename(p)))


if __name__ == '__main__':
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    print(main(*sys.argv[1:5]))
