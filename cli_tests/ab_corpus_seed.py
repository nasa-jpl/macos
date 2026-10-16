#!/usr/bin/env python3
"""Corpus A/B for the FwdRoot change: for each deck, load it in the macos CLI
(model 256), run `opd nElt`, record nPassRays / RMS OPD / P-V / whether the
vertex-sheet note fired.  One CLI process per deck (a crash at load must not
take the rest down).  Prompt-driven: waits for 'MACOS>' instead of sleeping.
Usage: ab_corpus.py <binary> <decklist> <out.csv>"""
import os, pty, select, sys, time, signal, re, csv

binary, decklist, outcsv = sys.argv[1:4]
MODEL = os.environ.get('MODEL', '256')
DEADLINE = float(os.environ.get('DEADLINE', '60'))

def run_deck(deck):
    pid, fd = pty.fork()
    if pid == 0:
        os.execv(binary, [binary])
    out = b''
    t0 = time.time()
    def wait_for(pats, tmax):
        nonlocal out
        end = time.time() + tmax
        while time.time() < end:
            r, _, _ = select.select([fd], [], [], 0.05)
            if r:
                try:
                    d = os.read(fd, 65536)
                except OSError:
                    return None
                if not d:
                    return None
                out += d
            tail = out[-4000:].decode(errors='replace')
            for k, p in enumerate(pats):
                if re.search(p, tail):
                    return k
        return -1
    status = 'ok'
    try:
        if wait_for([r'size.*:'], 10) is None: raise RuntimeError('died0')
        os.write(fd, (MODEL + '\n').encode())
        if wait_for([r'MACOS>'], 20) != 0: raise RuntimeError('noprompt')
        n0 = len(out)
        os.write(fd, ('old ' + deck + '\n').encode())
        k = wait_for([r'MACOS>\s*$'], DEADLINE)
        if k is None: raise RuntimeError('crash_load')
        if k == -1: raise RuntimeError('timeout_load')
        txt = out[n0:].decode(errors='replace')
        m = re.search(r'Tracing rays past\s+(\d+)\s+elements', txt)
        if not m: raise RuntimeError('load_fail')
        nelt = int(m.group(1))
        n1 = len(out)
        os.write(fd, ('opd %d\n' % nelt).encode())
        k = wait_for([r'Graphics device.*:', r'MACOS>\s*$'], DEADLINE)
        if k is None: raise RuntimeError('crash_opd')
        if k == -1: raise RuntimeError('timeout_opd')
        if k == 0:
            os.write(fd, b'/null\n')
            k = wait_for([r'MACOS>\s*$'], 20)
        os.write(fd, b'quit\n')
        wait_for([r'NEVER'], 2)
    except RuntimeError as e:
        status = str(e)
    try:
        os.kill(pid, signal.SIGKILL)
    except OSError:
        pass
    try:
        os.waitpid(pid, 0)
    except OSError:
        pass
    txt = out.decode(errors='replace')
    g = lambda p: (re.search(p, txt).group(1) if re.search(p, txt) else '')
    return dict(deck=deck, status=status, nelt=g(r'Tracing rays past\s+(\d+)'),
                npass=g(r'nPassRays =\s*(\d+)'), rms=g(r'RMS OPD error is\s+(\S+)'),
                pv=g(r'P-V OPD error is\s+(\S+)'), lost=g(r'A total of\s+(\d+)\s+rays were lost'),
                note=('1' if 'VERTEX sheet' in txt else '0'),
                flips=g(r'differing from the legacy metric:\s*(\d+)'),
                dt='%.1f' % (time.time() - t0))

decks = [l.strip() for l in open(decklist) if l.strip()]
with open(outcsv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['deck','status','nelt','npass','rms','pv','lost','note','flips','dt'])
    w.writeheader()
    for i, d in enumerate(decks):
        r = run_deck(d)
        w.writerow(r); f.flush()
        print('%3d/%d %-12s n%-3s pass %-5s rms %-16s note %s %ss  %s' % (i+1, len(decks), r['status'], r['nelt'], r['npass'], r['rms'], r['note'], r['dt'], os.path.basename(d)), flush=True)
