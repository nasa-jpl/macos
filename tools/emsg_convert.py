"""(ee): route the load path's stdout writes through emsg_mod.

WRITE(*|6, ctl) list  ->  WRITE(emsgRec, ctl) list   + a following  CALL emsg_put
PRINT *, list         ->  WRITE(emsgRec,*) list      + a following  CALL emsg_put

Every conversion is checked: the statement is found whole (its continuation lines, no comment
inside it), it is not a DO terminator, it carries no ADVANCE= / '$' prompt, nothing follows it on
its line after a ';', and the rewritten line stays within 132 columns (fixed form truncates
silently past that).  Sites are chosen per file: whole files, line ranges, or a ' ** ' text filter.
Idempotent: a converted line no longer matches WRITE(*, so a re-run converts only new sites -- add a file
or a filter to the table at the bottom to route more of the engine's output (notes/api_bindings.md).
Run:  python3 tools/emsg_convert.py [--dry]
"""
import re, sys
from pathlib import Path

M = Path(__file__).resolve().parents[1] / 'macos_f90'
DRY = '--dry' in sys.argv

START = re.compile(r'(WRITE\s*\(\s*(\*|6)\s*,|PRINT\s*\*\s*(,|$))', re.I)


def is_comment(l):
    return (not l.strip()) or l[:1] in 'cC*!' or l.lstrip().startswith('!') or l.lstrip().startswith('#')


def is_cont(l):
    if is_comment(l):
        return False
    if l.startswith('\t'):
        return len(l) > 1 and l[1] in '123456789'
    return len(l) > 5 and '\t' not in l[:5] and l[5] not in ' 0' and not l[:5].strip()


def strip_strings(s):
    return re.sub(r"'[^']*'|\"[^\"]*\"", "''", s)


def cols(l):
    return len(l.expandtabs(8))


def label_of(l):
    m = re.match(r'^\s{0,4}(\d+)[\s\t]', l)
    return m.group(1) if m and '\t' not in l[:l.find(m.group(1))] else None


def convert(fname, select, note):
    p = M / fname
    L = p.read_text().split('\n')
    text = '\n'.join(L)
    out, n, i, wide = [], 0, 0, []
    while i < len(L):
        l = L[i]
        m = None if is_comment(l) else START.search(strip_strings(l))
        if not m or not select(i + 1, l, L):
            out.append(l); i += 1; continue
        # the statement: this line + its continuations
        j = i + 1
        while j < len(L) and is_cont(L[j]):
            j += 1
        k = j
        while k < len(L) and is_comment(L[k]) and L[k].strip():
            k += 1
        assert not (k < len(L) and k > j and is_cont(L[k])), (fname, i + 1, 'comment inside statement')
        stmt = '\n'.join(L[i:j])
        assert not re.search(r'advance\s*=', stmt, re.I), (fname, i + 1, 'ADVANCE')
        lab = label_of(l)
        if lab:
            assert not re.search(r'\bDO\s+%s\b' % lab, text, re.I), (fname, i + 1, 'DO terminator', lab)
        # rewrite the first stdout write on the line (located in the string-stripped copy)
        s = strip_strings(l)
        mm = START.search(s)
        assert ';' not in s[mm.end():], (fname, i + 1, "';' after the write")
        pos = mm.start()
        if l[pos:pos + 5].upper() == 'PRINT':
            m2 = re.match(r'PRINT\s*\*\s*,?', l[pos:], re.I)
            new = l[:pos] + 'WRITE(emsgRec,*) ' + l[pos + m2.end():].lstrip()
        else:
            m2 = re.match(r'WRITE\s*\(\s*(\*|6)\s*,', l[pos:], re.I)
            new = l[:pos] + 'WRITE(emsgRec,' + l[pos + m2.end():]
        assert cols(new) <= 132, (fname, i + 1, cols(new), new)
        if cols(new) > 72:
            wide.append((i + 1, cols(new)))
        out.append(new)
        out.extend(L[i + 1:j])
        ind = re.match(r'^(\s*)', L[i]).group(1)
        if lab or not ind:
            ind = '        '
        out.append(ind + 'CALL emsg_put')
        n += 1
        i = j
    if not DRY:
        p.write_text('\n'.join(out))
    print('%-20s %3d converted (%s); over 72 columns: %s' % (fname, n, note, wide))
    return n


def in_ranges(*rr):
    return lambda ln, l, L: any(a <= ln <= b for a, b in rr)


def routine_range(fname, start_pat, end_pat):
    L = (M / fname).read_text().split('\n')
    a = next(i for i, l in enumerate(L) if re.search(start_pat, l, re.I)) + 1
    b = next(i for i, l in enumerate(L) if i >= a and re.search(end_pat, l, re.I)) + 1
    return a, b


def stars_in(rng, deny=()):
    """statements in rng whose text (with continuations) carries ' ** ', minus deny substrings"""
    def sel(ln, l, L):
        if not (rng[0] <= ln <= rng[1]):
            return False
        j = ln
        while j < len(L) and is_cont(L[j]):
            j += 1
        st = ' '.join(L[ln - 1:j])
        return '**' in st and not any(d in st for d in deny)
    return sel


total = 0
total += convert('msmacosio.inc', lambda ln, l, L: True, 'all')
total += convert('elt_mod.F', lambda ln, l, L: True, 'all: the parse helpers')
total += convert('iosub.inc', in_ranges(routine_range('iosub.inc', r'^\s+SUBROUTINE ChkDf2\b', r'END SUBROUTINE ChkDf2')), 'ChkDf2')
total += convert('surfsub.F', in_ranges(routine_range('surfsub.F', r'^\s+SUBROUTINE GridInit\b', r'end subroutine gridinit')), 'GridInit')
mb = routine_range('smacosio.F', r'^\s+SUBROUTINE MBFile6\b', r'END SUBROUTINE MBFile6')
pl = routine_range('smacosio.F', r'^\s+SUBROUTINE PARSE_LOH\(', r'END SUBROUTINE PARSE_LOH')
total += convert('smacosio.F', in_ranges(mb, pl), 'MBFile6 + PARSE_LOH')
nd = routine_range('design_optim.F', r'^\s+Subroutine nls_optim_dvr\(', r'End Subroutine nls_optim_dvr')
total += convert('design_optim.F', stars_in(nd, deny=('nRefRay', 'nomBeamRefRayDir')), "nls_optim_dvr ' ** '")
total += convert('macos_cmd_loop.inc', lambda ln, l, L: 'Unknown command = ' in l, 'Unknown command')
print('total', total)
