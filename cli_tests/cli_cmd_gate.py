#!/usr/bin/env python3
"""CLI command gate (PLAN_CONSOLIDATION sec. 2.2): one pty-driven journal per
command-reference entry.  A journal is a text file of CLI input lines, one per line,
sent in order, each after the CLI shows a prompt (the MACOS> prompt or any sub-prompt
ending in ':' / ']:' / '?'), plus directives:

  # MODEL: 256            model size answered to the first prompt (default 256)
  # DEADLINE: 90          seconds allowed per line (default 60)
  # EXPECT: <regex>       must match somewhere in the transcript (multi-line, re.M)
  # EXPECT-NOT: <regex>   must NOT match anywhere
  # COPY: <path>          copy that deck into the scratch directory; $DECK names the copy
                          (commands that WRITE beside the deck -- a TEXT spot, SAVE --
                          must never run on a deck inside the repos)
  # <anything else>       a comment
  $ROOT                   expands to the macos repo root; $RES to MACOS_resources;
  $TMP                    to a per-run scratch directory (SAVE targets)

Status per journal: pass | fail (an expectation missed) | crash (the CLI died) |
timeout.  The record is <out.csv> with the missed expectations; the transcript of
every journal is written beside it (<out_dir>/<journal>.txt) so a failure can be read.

Usage: cli_cmd_gate.py <binary> <journal_dir> <out.csv>
Non-vacuity: journal_dir/must_fail/*.jou are run too and must NOT pass (a journal
whose EXPECT can never match proves the runner reads expectations).
"""
import csv, glob, os, re, shutil, sys, tempfile, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cli_load_gate import Cli  # noqa: E402

PROMPT = r'(MACOS>\s*$|\]:\s*$|\?\s*$|:\s*$)'


def run_journal(binary, path, root, res, tmp):
    lines = open(path).read().splitlines()
    model, deadline, expects, nots, cmds, deck = '256', 60.0, [], [], [], ''
    for l in lines:
        if l.startswith('# COPY:'):
            src = l.split(':', 1)[1].strip().replace('$ROOT', root).replace('$RES', res)
            deck = os.path.join(tmp, os.path.basename(src))
            shutil.copy(src, deck)
        elif l.startswith('# MODEL:'):
            model = l.split(':', 1)[1].strip()
        elif l.startswith('# DEADLINE:'):
            deadline = float(l.split(':', 1)[1])
        elif l.startswith('# EXPECT-NOT:'):
            nots.append(l.split(':', 1)[1].strip())
        elif l.startswith('# EXPECT:'):
            expects.append(l.split(':', 1)[1].strip())
        elif l.startswith('#'):
            continue
        else:
            cmds.append(l.replace('$ROOT', root).replace('$RES', res).replace('$TMP', tmp).replace('$DECK', deck))
    t0 = time.time()
    status, missed = 'pass', []
    try:
        cli = Cli(binary, model, deadline)
    except RuntimeError as e:
        return str(e), [], '', time.time() - t0
    try:
        for c in cmds:
            n0 = len(cli.out)
            cli.send(c)
            k = cli.wait([PROMPT], deadline)
            if k is None:
                status = 'crash'
                break
            if k == -1:
                status = 'timeout'
                break
            since = cli.out[n0:].decode(errors='replace')
            if 'does not exist. Check name/path' in since or 'Pick a different file' in since:
                # a LOAD re-prompt would eat every following line: stop here, say why
                status = 'fail'; missed.append('LOAD refused or file missing: ' + c)
                cli.send('q')
                break
    finally:
        txt = cli.text()
        cli.close()
    if status == 'pass':
        for e in expects:
            if not re.search(e, txt, re.M):
                missed.append('EXPECT ' + e)
        for e in nots:
            if re.search(e, txt, re.M):
                missed.append('EXPECT-NOT ' + e)
        if missed:
            status = 'fail'
    return status, missed, txt, time.time() - t0


def main():
    binary, jdir, outcsv = sys.argv[1:4]
    root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    res = os.environ.get('MACOS_RESOURCES', os.path.abspath(os.path.join(root, '..', 'MACOS_resources')))
    tmp = tempfile.mkdtemp(prefix='cli_cmd_')
    tdir = os.path.splitext(outcsv)[0] + '_transcripts'
    os.makedirs(tdir, exist_ok=True)
    jous = sorted(glob.glob(os.path.join(jdir, '*.jou'))) + sorted(glob.glob(os.path.join(jdir, 'must_fail', '*.jou')))
    with open(outcsv, 'w', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['journal', 'kind', 'status', 'missed', 'dt'])
        for j in jous:
            kind = 'must_fail' if os.sep + 'must_fail' + os.sep in j else 'journal'
            st, missed, txt, dt = run_journal(binary, j, root, res, tmp)
            name = os.path.basename(j)
            open(os.path.join(tdir, name.replace('.jou', '.txt')), 'w').write(txt)
            w.writerow([name, kind, st, ' | '.join(missed), '%.1f' % dt]); f.flush()
            print('%-28s %-9s %-8s %5.1fs  %s' % (name, kind, st, dt, '; '.join(missed)[:90]), flush=True)
    shutil.rmtree(tmp, ignore_errors=True)


if __name__ == '__main__':
    main()
