#!/usr/bin/env python3
"""Judge a load-gate record: judge_load.py <record.csv> <mustfail.csv> [<previous.csv>]
Prints the summary; exit 1 on a new failure (see run_cli_tests.sh header)."""
import collections, csv, os, sys


def read(p):
    if not p or not os.path.exists(p):
        return {}
    return {r['deck']: r for r in csv.DictReader(open(p))}


rec, mf, prev = sys.argv[1], sys.argv[2], (sys.argv[3] if len(sys.argv) > 3 else '')
R, M, P = read(rec), read(mf), read(prev)
bad = 0
print('record  %s  (%d decks)' % (os.path.basename(rec), len(R)))

# (a) non-vacuity: every must-fail deck must NOT load
for d, r in M.items():
    ok = r['status'] != 'ok'
    print('  must-fail %-22s %-12s %s' % (os.path.basename(d), r['status'], 'OK' if ok else 'LOADED -- the gate is blind'))
    bad |= (not ok)

# status counts
c = collections.Counter(r['status'] for r in R.values())
rt = collections.Counter(r['roundtrip'] for r in R.values() if r['status'] == 'ok')
print('  status    ' + '  '.join('%s %d' % kv for kv in sorted(c.items())))
print('  roundtrip ' + '  '.join('%s %d' % kv for kv in sorted(rt.items())))
nan = [d for d, r in R.items() if r['status'] == 'ok' and 'NaN' in r['rms']]
if nan:
    print('  NaN OPD   %d: %s' % (len(nan), ' '.join(os.path.basename(d) for d in nan)))

# pre-existing failures: named, not fatal
for d, r in sorted(R.items()):
    if r['status'] != 'ok':
        was = P.get(d, {}).get('status', '')
        tag = 'NEW' if (P and was == 'ok') else ('pre-existing' if not P or was else 'new deck')
        print('  %-13s %-12s %s' % (r['status'], tag, d))
        if tag == 'NEW':
            bad = 1
    elif r['roundtrip'] not in ('same', 'ulp', 'n/a'):
        was = P.get(d, {}).get('roundtrip', '')
        tag = 'NEW' if (P and was in ('same', 'ulp')) else ('pre-existing' if not P or was else 'new deck')
        print('  rt-%-10s %-12s %s' % (r['roundtrip'], tag, d))
        if tag == 'NEW':
            bad = 1

# (c) numbers that moved against the previous record (informational unless asked)
if P:
    moved = [d for d, r in R.items() if d in P and r['status'] == 'ok' and P[d]['status'] == 'ok'
             and (r['rms'], r['pv'], r['npass'], r['lost']) != (P[d]['rms'], P[d]['pv'], P[d]['npass'], P[d]['lost'])]
    print('  vs %s: %d decks with a changed OPD/nPass/lost' % (os.path.basename(prev), len(moved)))
    for d in moved[:40]:
        print('    moved  %s  rms %s -> %s  pass %s -> %s' % (d, P[d]['rms'], R[d]['rms'], P[d]['npass'], R[d]['npass']))
print('  verdict   %s' % ('FAIL' if bad else 'PASS'))
sys.exit(1 if bad else 0)
