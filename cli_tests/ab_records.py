#!/usr/bin/env python3
"""The A/B gate (PLAN_CONSOLIDATION sec. 2.3): the load-gate record at a PRE commit
against the record at a POST commit, deck by deck.

    ab_records.py <pre.csv> <post.csv> [--allow <movers.txt>] [--out <record.txt>]

A deck MOVES when any of status / nElt / nPass / lost / RMS / P-V / avg / roundtrip
differs between the two records (numbers compared as printed -- the CLI prints
10 digits; a 1-ulp wobble in the SAVE round trip is the 'ulp' class, already
folded into the roundtrip column).  Every mover is printed with the columns that
moved and written to --out in the form of cli_tests/records/ab_*_moved_*.txt (the
header lines are yours to write: what the change was MEANT to move, and why).

Exit 1 when a deck moved that is NOT named in --allow (any file with the movers' names or
paths on its lines, '#' comments -- an ab_*_moved record itself serves) -- the commit must name every mover.  Exit 1 also when a
deck went ok -> not ok (a mover that an allow list cannot excuse).

Where the records come from: run_cli_tests.sh load|corpus from the clean worktree
(`~/dev/macos_cli_base`, `git checkout <sha>`, `makems.sh release [gfortran]`) at
each commit -- records/<kind>_<compiler>_<sha>.csv.  The builds stay a manual
step: one clean worktree holds one SHA at a time, and the record is what gets
committed.  Run both compilers' pairs.
"""
import argparse, csv, os, sys

COLS = ['status', 'nelt', 'npass', 'lost', 'rms', 'pv', 'avg', 'roundtrip']


def read(p):
    return {r['deck']: r for r in csv.DictReader(open(p))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pre'); ap.add_argument('post')
    ap.add_argument('--allow', default='')
    ap.add_argument('--out', default='')
    a = ap.parse_args()
    A, B = read(a.pre), read(a.post)
    allow = set()
    if a.allow:
        for l in open(a.allow):
            l = l.split('#', 1)[0].strip()
            for t in l.split():         # any deck token: a movers record's own lines serve
                if t.endswith('.in'):
                    allow.add(t); allow.add(os.path.basename(t))
    moved, broke, only = [], [], []
    for d in A:
        if d not in B:
            only.append(d); continue
        x, y = A[d], B[d]
        cols = [c for c in COLS if x.get(c, '') != y.get(c, '')]
        if cols:
            moved.append((d, cols, x, y))
            if x['status'] == 'ok' and y['status'] != 'ok':
                broke.append(d)
    lines = ['# A/B: pre %s  post %s  (%d decks in common, %d only in pre, %d only in post)' % (
        os.path.basename(a.pre), os.path.basename(a.post), len(A) - len(only), len(only), len([d for d in B if d not in A])),
        '# columns compared: ' + ' '.join(COLS), '# %d decks moved' % len(moved)]
    for d, cols, x, y in moved:
        lines.append('%s  %s' % (d, '  '.join('%s %s -> %s' % (c, x.get(c, ''), y.get(c, '')) for c in cols)))
    txt = '\n'.join(lines) + '\n'
    print(txt, end='')
    if a.out:
        open(a.out, 'w').write(txt)
    unexplained = [d for d, _, _, _ in moved if d not in allow and os.path.basename(d) not in allow]
    bad = 0
    if broke:
        print('FAIL: %d decks went ok -> not ok: %s' % (len(broke), ' '.join(os.path.basename(d) for d in broke))); bad = 1
    if unexplained:
        print('FAIL: %d movers not in the allow list: %s' % (len(unexplained), ' '.join(os.path.basename(d) for d in unexplained[:30]))); bad = 1
    if not bad:
        print('PASS: %d movers, all named' % len(moved) if moved else 'PASS: 0 movers')
    sys.exit(bad)


if __name__ == '__main__':
    main()
