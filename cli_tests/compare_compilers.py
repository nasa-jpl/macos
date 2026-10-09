#!/usr/bin/env python3
"""compare_compilers.py <ifx.csv> <gfortran.csv>: decks the two compilers disagree on."""
import csv, sys


def f(x):
    try:
        return float(x.replace('D', 'E').replace('d', 'e'))
    except ValueError:
        return float('nan')


A = {r['deck']: r for r in csv.DictReader(open(sys.argv[1]))}
B = {r['deck']: r for r in csv.DictReader(open(sys.argv[2]))}
out = []
for d in A:
    if d not in B:
        continue
    a, b = A[d], B[d]
    if (a['status'], a['npass'], a['lost']) != (b['status'], b['npass'], b['lost']):
        out.append((d, 'status/pass/lost', a['status'], b['status'], a['npass'], b['npass'], a['lost'], b['lost']))
    elif a['status'] == 'ok' and a['rms'] and b['rms']:
        ra, rb = f(a['rms']), f(b['rms'])
        if ra == ra and rb == rb and abs(ra - rb) > 1e-9 * max(abs(ra), abs(rb), 1e-300):
            out.append((d, 'rms', a['rms'], b['rms']))
print('ifx vs gfortran: %d of %d decks differ' % (len(out), len(A)))
for o in out:
    print('  ' + '  '.join(str(x) for x in o))
