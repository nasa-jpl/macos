#!/usr/bin/env python3
"""The prompt-and-default inventory behind the *Syntax:* / *Defaults:* lines of Part I
(2026-10-17): for every engine command, the command lines the CLI suite's journals issued and
the prompts (with their bracketed defaults) the CLI answered them with, read from one run's
transcripts.  Usage (from docs/macos-manual):
    python3 tools/cmd_prompt_inventory.py [../../cli_tests/records/cmd_ifx_<sha>_transcripts] > inventory.md
Rerun it after the engine changes a dialog, and compare the entry's lines."""
import re, sys, glob, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_cmdref as g
T = sys.argv[1] if len(sys.argv) > 1 else sorted(glob.glob('../../cli_tests/records/cmd_ifx_*_transcripts'))[-1]
cats = g.parse_help_inc()
entries = []
for cat, items in cats:
    for names, desc in items:
        forms = []
        for nm in names.split(','):
            w = nm.strip().split()[0] if nm.strip() else ''
            w = re.sub(r'[<>\[\]]', '', w)
            if not w:
                continue
            mm = re.match(r'^[A-Z!]+', w); mm = mm.group(0) if mm else w.upper()
            forms.append((mm.upper(), w.upper()))
        entries.append((cat, names, g._entry_key(names), forms))
def match(tok):
    tu = tok.upper(); hits = []
    for cat, names, key, forms in entries:
        for mm, full in forms:
            if len(tu) >= len(mm) and full.startswith(tu):
                hits.append(key); break
    return hits
inv = {}
for f in sorted(glob.glob(T + '/*.txt')):
    L = open(f, errors='replace').read().split('\n'); i = 0
    while i < len(L):
        m = re.match(r'^\s*MACOS>\s*(.*)$', L[i])
        if m and m.group(1).strip():
            cmdline = m.group(1).strip(); tok = cmdline.split()[0]
            if not re.match(r'^[A-Za-z!][A-Za-z]*$', tok):
                i += 1; continue
            prompts = []; j = i + 1
            while j < len(L) and not re.match(r'^\s*MACOS>', L[j]):
                pm = re.match(r'^\s*(.*?)\s*\[([^\]]*)\]:\s*(.*)$', L[j])
                if pm and ('Enter' in pm.group(1) or '?' in pm.group(1) or pm.group(1).strip().endswith(':') or len(pm.group(1)) < 70):
                    prompts.append((pm.group(1).strip().rstrip(':').strip(), pm.group(2)))
                else:
                    pm2 = re.match(r'^\s*(Enter [^\[]*?):\s*(\S.*)?$', L[j])
                    if pm2:
                        prompts.append((pm2.group(1).strip(), '(no default)'))
                j += 1
            for key in match(tok):
                e = inv.setdefault(key, {'lines': set(), 'seqs': {}})
                e['lines'].add(cmdline[:60])
                e['seqs'].setdefault(tuple(prompts), set()).add(os.path.basename(f)[:-4])
            i = j
        else:
            i += 1
print(f'# Prompt inventory from {T}')
for cat, items in cats:
    print(f'\n# {cat}')
    for names, desc in items:
        key = g._entry_key(names); e = inv.get(key)
        print(f'\n## {key}  ({names})')
        if not e:
            print('  (no invocation in the transcripts)'); continue
        print('  lines: ' + ' | '.join(sorted(e['lines'])[:8]))
        for sig, jous in sorted(e['seqs'].items(), key=lambda kv: -len(kv[1]))[:5]:
            print('  prompts (' + ','.join(sorted(jous))[:60] + '): ' + (' ; '.join(f'{p} [{d}]' for p, d in sig) if sig else '(none)'))
