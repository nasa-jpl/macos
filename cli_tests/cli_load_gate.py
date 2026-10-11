#!/usr/bin/env python3
"""CLI load gate (PLAN_CONSOLIDATION sec. 2.1): for every deck in a list, drive the
macos CLI through a pty (readline needs a tty) and record

    LOAD -> opd nElt -> SAVE s1 ;  new process: LOAD s1.in -> SAVE s2 ;
    new process: LOAD s2.in -> SAVE s3 ;   roundtrip = s2 vs s3,  settle = s1 vs s2

one CLI process per leg so a crash at load cannot take the rest down.  Prompt-
driven (waits for 'MACOS>'), never sleeps.  Each deck yields one CSV row:

  deck, status, nelt, npass, rms, pv, avg, lost, flips, rt_status, rt_nelt, roundtrip, settle, dt

status      ok | crash_load | timeout_load | load_fail | crash_opd | timeout_opd
            (crash_* = the process died: the host-killer class; load_fail includes
            the CLI's own "does not exist" re-prompt, answered with q)
roundtrip   same | ulp | eol | differ | n/a   -- s3 against s2: the SAVE fixed point
            the SAVE work asserts.  A hand-written deck's unit vectors are unitised
            ONCE at load, so s2 may differ from s1 by that settle (6a, 2026-10-10:
            102 of 138 'ulp' decks settle once; the 33 that DRIFTED were the zGrid
            defect) -- which is why the fixed point is judged on the THIRD save.
settle      same | ulp | eol | differ | n/a   -- s2 against s1, informational: 'ulp'
            here is the documented one-time settle; 'differ' is a SAVE that rewrites
            something on its first reload (name it)
flips       the FwdRoot 'differing from the legacy metric' count printed by WARN
avg         the CLI's 'Average OPD is' (DAvgl: the mean of ray-minus-chief OPL over the
            valid rays) -- the same under either OPD reference; it is exactly the
            constant by which the chief-referenced map differs from the mean-referenced
            one (2026-10-10), so a nonzero avg names a deck whose map moved

Usage:  cli_load_gate.py <binary> <decklist> <out.csv> [--model N] [--deadline S]
        [--workdir DIR]   (SAVE files go under DIR; default a temp dir, removed)
Exit code 0 always: the RECORD is the result; run_cli_tests.sh judges it.

Lineage: the FwdRoot corpus A/B driver (2026-10-04, scratchpad ab_corpus.py), kept
as cli_tests/ab_corpus_seed.py.
"""
import argparse, csv, os, pty, re, select, shutil, signal, tempfile, time


class Cli:
    """One macos CLI process on a pty."""

    def __init__(self, binary, model, deadline):
        self.deadline = deadline
        self.pid, self.fd = pty.fork()
        if self.pid == 0:
            os.execv(binary, [binary])
        self.out = b''
        if self.wait([r'size.*:'], 15) is None:
            raise RuntimeError('died0')
        self.send(model)
        if self.wait([r'MACOS>'], 30) != 0:
            raise RuntimeError('noprompt')

    def send(self, line):
        os.write(self.fd, (line + '\n').encode())

    def wait(self, pats, tmax):
        """Return the index of the first pattern seen in the tail, -1 on timeout,
        None if the process died."""
        end = time.time() + tmax
        while time.time() < end:
            r, _, _ = select.select([self.fd], [], [], 0.05)
            if r:
                try:
                    d = os.read(self.fd, 65536)
                except OSError:
                    return None
                if not d:
                    return None
                self.out += d
            tail = self.out[-6000:].decode(errors='replace')
            for k, p in enumerate(pats):
                if re.search(p, tail):
                    return k
        return -1

    def cmd(self, line, extra_pats=(), tmax=None):
        """Send a command and wait for the next prompt; returns (k, text since)."""
        n0 = len(self.out)
        self.send(line)
        k = self.wait([r'MACOS>\s*$'] + list(extra_pats), tmax or self.deadline)
        return k, self.out[n0:].decode(errors='replace')

    def text(self):
        return self.out.decode(errors='replace')

    def close(self):
        try:
            self.send('quit')
            self.wait([r'NEVER'], 1.5)
        except OSError:
            pass
        for sig in (signal.SIGTERM, signal.SIGKILL):
            try:
                os.kill(self.pid, sig)
            except OSError:
                pass
        try:
            os.waitpid(self.pid, 0)
        except OSError:
            pass
        try:
            os.close(self.fd)           # the pty fd: without this, select() fails at deck 511 (fd > 1023)
        except OSError:
            pass


def classify_diff(a, b):
    """'eol' when the two saves differ only in line endings; 'ulp' when only in numbers within 8 ulp (the known
    1-ulp unit-vector print wobble: psiElt / xGrid on tilted elements -- a SAVE ->
    load -> SAVE 2-cycle, PLAN sec. 0 candidate), else 'differ'."""
    la, lb = a.decode(errors='replace').splitlines(), b.decode(errors='replace').splitlines()
    if la == lb:
        return 'eol'                    # line endings only (CRLF vs LF: a Windows save)
    if len(la) != len(lb):
        return 'differ'
    num = re.compile(r'^[-+]?(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?$')
    for x, y in zip(la, lb):
        if x == y:
            continue
        tx, ty = x.split(), y.split()
        if len(tx) != len(ty):
            return 'differ'
        for u, v in zip(tx, ty):
            if u == v:
                continue
            if not (num.match(u) and num.match(v)):
                return 'differ'
            fu = float(u.replace('D', 'E').replace('d', 'e'))
            fv = float(v.replace('D', 'E').replace('d', 'e'))
            if abs(fu - fv) > 8 * 2.3e-16 * max(abs(fu), abs(fv), 1e-300):
                return 'differ'
    return 'ulp'


def g(txt, pat, default=''):
    m = re.search(pat, txt)
    return m.group(1) if m else default


def load(cli, deck):
    """LOAD a deck; returns (status, nelt, text).  A validator refusal re-prompts
    for a file name ('Pick a different file (or "q" to abort)'); answer q."""
    k, txt = cli.cmd('old ' + deck, extra_pats=[r'Pick a different file',
                                                  r'does not exist\. Check name/path'])
    if k in (1, 2):                     # a re-prompt would otherwise read as a timeout
        k2, txt2 = cli.cmd('q')
        txt += txt2
        return ('load_fail' if k2 is not None else 'crash_load'), '', txt
    if k is None:
        return 'crash_load', '', txt
    if k == -1:
        return 'timeout_load', '', txt
    m = re.search(r'Tracing rays past\s+(\d+)\s+elements', txt)
    if not m:
        return 'load_fail', '', txt
    return 'ok', m.group(1), txt


def save(cli, stem):
    """SAVE to <stem>.in (stem must not exist); returns (ok, text)."""
    k, txt = cli.cmd('save ' + stem, extra_pats=[r'Replace\?'])
    if k == 1:                      # should not happen with a fresh stem
        cli.send('yes')
        k = cli.wait([r'MACOS>\s*$'], cli.deadline)
    return (k == 0 and os.path.exists(stem + '.in')), txt


def run_deck(binary, deck, model, deadline, workdir, idx):
    t0 = time.time()
    row = dict(deck=deck, status='', nelt='', npass='', rms='', pv='', avg='', lost='',
               flips='', rt_status='', rt_nelt='', roundtrip='n/a', settle='n/a', dt='')
    if not os.path.isfile(deck):            # the CLI re-prompts forever on a missing file
        row['status'] = 'missing'; row['dt'] = '0.0'
        return row
    s1 = os.path.join(workdir, 'rt%04d_a' % idx)
    s2 = os.path.join(workdir, 'rt%04d_b' % idx)
    s3 = os.path.join(workdir, 'rt%04d_c' % idx)
    # ---- leg 1: load, opd, save
    try:
        cli = Cli(binary, model, deadline)
    except RuntimeError as e:
        row['status'] = str(e)
        row['dt'] = '%.1f' % (time.time() - t0)
        return row
    try:
        st, nelt, txt = load(cli, deck)
        row['status'], row['nelt'] = st, nelt
        if st == 'ok':
            k, txt = cli.cmd('opd ' + nelt, extra_pats=[r'Graphics device.*:'])
            if k == 1:
                k, txt2 = cli.cmd('/null')
                txt += txt2
            if k is None:
                row['status'] = 'crash_opd'
            elif k == -1:
                row['status'] = 'timeout_opd'
            else:
                row['npass'] = g(txt, r'nPassRays =\s*(\d+)')
                row['rms'] = g(txt, r'RMS OPD error is\s+(\S+)')
                row['pv'] = g(txt, r'P-V OPD error is\s+(\S+)')
                row['avg'] = g(txt, r'Average OPD is\s+(\S+)')
                row['lost'] = g(txt, r'A total of\s+(\d+)\s+rays were lost')
                row['flips'] = g(txt, r'differing from the legacy metric:\s*(\d+)')
                ok1, _ = save(cli, s1)
                if not ok1:
                    row['roundtrip'] = 'save1_fail'
    finally:
        cli.close()
    # ---- legs 2 and 3: reload the save, save again (twice); settle = s1 vs s2,
    #      roundtrip = s2 vs s3 (the fixed point after the one-time settle)
    def cmp(x, y):
        a = open(x + '.in', 'rb').read(); b = open(y + '.in', 'rb').read()
        return 'same' if a == b else classify_diff(a, b)
    for src, dst, col, failtag in ((s1, s2, 'settle', 'save2_fail'), (s2, s3, 'roundtrip', 'save3_fail')):
        if row['status'] != 'ok' or not os.path.exists(src + '.in'):
            break
        try:
            cli = Cli(binary, model, deadline)
            try:
                st2, nelt2, _ = load(cli, src + '.in')
                if col == 'settle':
                    row['rt_status'], row['rt_nelt'] = st2, nelt2
                elif st2 != 'ok':
                    row['rt_status'] = 'rt3_' + st2
                if st2 == 'ok':
                    ok2, _ = save(cli, dst)
                    row[col] = cmp(src, dst) if ok2 else failtag
            finally:
                cli.close()
        except RuntimeError as e:
            row['rt_status'] = str(e)
            break
        if row[col] in (failtag,):
            break
    row['dt'] = '%.1f' % (time.time() - t0)
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('binary'); ap.add_argument('decklist'); ap.add_argument('outcsv')
    ap.add_argument('--model', default='256')
    ap.add_argument('--deadline', type=float, default=90.0)
    ap.add_argument('--workdir', default='')
    a = ap.parse_args()
    decks = [l.strip() for l in open(a.decklist) if l.strip() and not l.startswith('#')]
    tmp = a.workdir or tempfile.mkdtemp(prefix='cli_load_')
    os.makedirs(tmp, exist_ok=True)
    fields = ['deck', 'status', 'nelt', 'npass', 'rms', 'pv', 'avg', 'lost', 'flips',
              'rt_status', 'rt_nelt', 'roundtrip', 'settle', 'dt']
    with open(a.outcsv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator='\n')
        w.writeheader()
        for i, d in enumerate(decks):
            r = run_deck(a.binary, d, a.model, a.deadline, tmp, i)
            w.writerow(r); f.flush()
            print('%3d/%d %-13s n%-3s pass %-5s rms %-16s rt %-6s settle %-6s %5ss  %s' % (
                i + 1, len(decks), r['status'], r['nelt'], r['npass'], r['rms'],
                r['roundtrip'], r['settle'], r['dt'], os.path.basename(d)), flush=True)
    if not a.workdir:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == '__main__':
    main()
