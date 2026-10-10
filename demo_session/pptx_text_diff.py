import sys, difflib
from pptx import Presentation

# Compares the TEXT of two decks, shape by shape, slide by slide.
#
# TABLES ARE INCLUDED, and they must stay included: until 2026-09-18 this tool
# read only `sh.has_text_frame`, so an edit made inside a table cell produced an
# EMPTY diff.  sync_edit_deck.sh gates on this tool's output, so a table-only
# edit pass looked like "no unrecovered changes" and was silently overwritten by
# the next sync.  A table cell is emitted as one "r<row>c<col>: <text>" line so a
# changed cell shows up as a one-line diff naming its own coordinates.


def shape_text(sh):
    if sh.has_text_frame and sh.text_frame.text.strip():
        return sh.text_frame.text
    if getattr(sh, "has_table", False) and sh.has_table:
        rows = []
        for ri, row in enumerate(sh.table.rows):
            for ci, cell in enumerate(row.cells):
                t = cell.text.strip()
                if t:
                    rows.append(f"r{ri}c{ci}: {t}")
        if rows:
            return "\n".join(rows)
    return ""


def texts(path):
    prs = Presentation(path)
    out = []
    for sl in prs.slides:
        st = []
        for sh in sl.shapes:
            t = shape_text(sh)
            if t:
                st.append(t)
        out.append(st)
    return out


a = texts(sys.argv[1]); b = texts(sys.argv[2])
for si in range(max(len(a), len(b))):
    A = a[si] if si < len(a) else []
    B = b[si] if si < len(b) else []
    for i in range(max(len(A), len(B))):
        ta = A[i] if i < len(A) else ""
        tb = B[i] if i < len(B) else ""
        if ta != tb:
            print(f"\n=== SLIDE {si+1} shape {i}:")
            for l in difflib.unified_diff(ta.split("\n"), tb.split("\n"),
                                          lineterm="", n=0):
                if l[:3] in ("---", "+++"): continue
                print(" " + l)
