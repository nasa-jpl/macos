"""
audit.py — Read-only structural diagnostic for the manual.

Walks working/word/document.xml in order and produces a report on
pagination artifacts and heterogeneous-style buckets, so a human can
decide reclassification rules before running 00_/05_/08_ patches.

Usage:
    python scripts/unpack.py macosMan4.1beta.docx
    python scripts/audit.py                          # prints to stdout
    python scripts/audit.py audit_4.1beta.txt        # writes a file

What it shows:
    - Every <w:sectPr> location in document order, with heading context
    - Every column break with heading context
    - For each "messy" style (TableParagraph, Normal, FrameContents):
        up to 30 representative samples spread across the document
        (idx, font, size, bold, in-table?, snippet)
    - For TableParagraph: separate counts of in-tbl vs out-of-tbl uses
"""
import sys
from pathlib import Path
from collections import Counter, defaultdict

sys.path.insert(0, str(Path(__file__).parent))
from docx_helpers import (load_doc, wtag, get_style, get_text,
                          get_run_info, iter_paragraphs)

PROJECT_ROOT = Path(__file__).parent.parent
WORKING_DIR  = PROJECT_ROOT / "working"

MESSY_STYLES = ["TableParagraph", "Normal", "FrameContents"]
SAMPLE_LIMIT = 30
SNIPPET_LEN  = 90


def is_in_table(para, parent_map):
    """Walk parents; True iff some ancestor is <w:tbl>."""
    cur = parent_map.get(id(para))
    while cur is not None:
        if cur.tag == wtag("tbl"):
            return True
        cur = parent_map.get(id(cur))
    return False


def build_parent_map(root):
    """Map id(elem) -> parent elem.  ET.ElementTree has no .getparent()."""
    pm = {}
    for parent in root.iter():
        for child in parent:
            pm[id(child)] = parent
    return pm


def walk_in_order(root):
    """
    Yield document-body elements in document order, paragraph index
    incrementing only on <w:p> elements.  This lets us locate sectPr /
    column-break / heading paragraphs by the same index axis.
    """
    body = None
    for c in root:
        if c.tag == wtag("body"):
            body = c
            break
    if body is None:
        return
    pidx = -1
    for elem in body.iter():
        if elem.tag == wtag("p"):
            pidx += 1
        yield pidx, elem


def heading_context(idx, heading_track):
    """Return latest H2/H3/H4 seen at or before idx, as a string."""
    parts = []
    for level in ("H2", "H3", "H4"):
        last = heading_track.get(level)
        if last is not None and last[0] <= idx:
            parts.append(f"{level}@{last[0]}: {last[1][:60]}")
    return " | ".join(parts) if parts else "<no heading yet>"


def main(out=None):
    if not WORKING_DIR.exists():
        sys.exit("working/ not found — run: python scripts/unpack.py")

    print = (out.write if out else __builtins__.print)
    def w(line=""):
        if out:
            out.write(str(line) + "\n")
        else:
            __builtins__.print(line)

    tree, root = load_doc(WORKING_DIR)
    parent_map = build_parent_map(root)

    # ── Section breaks & column breaks (in document order) ────────────────
    sectprs = []     # list of (paragraph_idx, location_kind)
    column_breaks = []
    heading_track = {}     # level -> (idx, text)
    style_counter = Counter()
    style_idx_buckets = defaultdict(list)  # style -> list of (idx, p)

    # Body-final sectPr is the document's main page-setup; identify it.
    body = next(c for c in root if c.tag == wtag("body"))
    body_sectpr = next((c for c in reversed(list(body))
                        if c.tag == wtag("sectPr")), None)

    for idx, elem in walk_in_order(root):
        if elem.tag == wtag("sectPr"):
            kind = "body-final" if elem is body_sectpr else "intermediate"
            sectprs.append((idx, kind, elem))
        elif elem.tag == wtag("br"):
            if elem.get(wtag("type")) == "column":
                column_breaks.append((idx, elem))
        elif elem.tag == wtag("p"):
            style = get_style(elem)
            style_counter[style] += 1
            if style == "Heading1":
                heading_track["H2"] = (idx, get_text(elem).strip())
            elif style == "Heading2":
                heading_track["H2"] = (idx, get_text(elem).strip())
                heading_track.pop("H3", None)
                heading_track.pop("H4", None)
            elif style == "Heading3":
                heading_track["H3"] = (idx, get_text(elem).strip())
                heading_track.pop("H4", None)
            elif style == "Heading4":
                heading_track["H4"] = (idx, get_text(elem).strip())
            if style in MESSY_STYLES:
                style_idx_buckets[style].append((idx, elem))

    # We can't show context "from the past" with a single pass + single
    # heading_track, so re-pass while tracking.  Alternative: reuse the
    # tracker by rebuilding for each location.  Simpler: precompute a
    # heading-at-idx map.
    heading_at = {}     # paragraph_idx -> dict
    htrk = {}
    for idx, elem in walk_in_order(root):
        if elem.tag == wtag("p"):
            style = get_style(elem)
            text = get_text(elem).strip()
            if style == "Heading1":
                htrk["H2"] = (idx, text); htrk.pop("H3", None); htrk.pop("H4", None)
            elif style == "Heading2":
                htrk["H2"] = (idx, text); htrk.pop("H3", None); htrk.pop("H4", None)
            elif style == "Heading3":
                htrk["H3"] = (idx, text); htrk.pop("H4", None)
            elif style == "Heading4":
                htrk["H4"] = (idx, text)
            heading_at[idx] = dict(htrk)

    def context_at(idx):
        # Latest snapshot at or before idx
        keys = [k for k in heading_at if k <= idx]
        snap = heading_at[max(keys)] if keys else {}
        parts = []
        for lvl in ("H2", "H3", "H4"):
            v = snap.get(lvl)
            if v:
                parts.append(f"{lvl}#{v[0]}: {v[1][:55]}")
        return "  ".join(parts) if parts else "<no heading>"

    # ── REPORT ─────────────────────────────────────────────────────────────
    w("=" * 78)
    w("MACOS Manual — Structural Audit")
    w("=" * 78)
    w(f"Source: {WORKING_DIR}")
    w(f"Total paragraphs: {sum(style_counter.values())}")
    w()

    w("Style counts:")
    for s, n in style_counter.most_common():
        w(f"  {n:5d}  {s}")
    w()

    w("-" * 78)
    w(f"SECTION BREAKS (<w:sectPr>): {len(sectprs)} total")
    w("-" * 78)
    intermediates = [s for s in sectprs if s[1] == "intermediate"]
    final = [s for s in sectprs if s[1] == "body-final"]
    w(f"  body-final (keep):    {len(final)}")
    w(f"  intermediate (strip): {len(intermediates)}")
    w()
    w("First 20 intermediate sectPr locations (by paragraph index):")
    for idx, kind, _ in intermediates[:20]:
        w(f"  para#{idx:5d}  {context_at(idx)}")
    if len(intermediates) > 20:
        w(f"  ... and {len(intermediates)-20} more")
    w()

    w("-" * 78)
    w(f"COLUMN BREAKS: {len(column_breaks)} total")
    w("-" * 78)
    w("First 20 column-break locations:")
    for idx, _ in column_breaks[:20]:
        w(f"  para#{idx:5d}  {context_at(idx)}")
    if len(column_breaks) > 20:
        w(f"  ... and {len(column_breaks)-20} more")
    w()

    # ── Messy-style buckets ────────────────────────────────────────────────
    for style in MESSY_STYLES:
        bucket = style_idx_buckets.get(style, [])
        w("-" * 78)
        w(f"{style}: {len(bucket)} paragraphs")
        w("-" * 78)
        if not bucket:
            w("  (none)")
            w()
            continue

        if style == "TableParagraph":
            in_tbl  = sum(1 for _, p in bucket if is_in_table(p, parent_map))
            out_tbl = len(bucket) - in_tbl
            w(f"  inside <w:tbl>:  {in_tbl}   <- legitimate, keep as 'TableCell'")
            w(f"  outside tables: {out_tbl}   <- needs reclassification")
            w()

        # Sample evenly across the document
        n = len(bucket)
        step = max(1, n // SAMPLE_LIMIT)
        samples = bucket[::step][:SAMPLE_LIMIT]
        w(f"  {len(samples)} samples (every {step}-th occurrence):")
        for idx, p in samples:
            text = get_text(p).strip().replace("\n", " ")
            info = get_run_info(p)
            in_tbl = is_in_table(p, parent_map)
            tbl_flag = "T" if in_tbl else "."
            w(f"  [{tbl_flag}] para#{idx:5d}  font={info['font'] or '<none>':14s}"
              f" sz={info['size']:>3s}  b={int(info['bold'])}"
              f"  | {text[:SNIPPET_LEN]}")
        w()

    w("=" * 78)
    w("End of audit.")
    w("=" * 78)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        out_path = Path(sys.argv[1])
        with out_path.open("w") as f:
            main(out=f)
        sys.stderr.write(f"Audit written to {out_path}\n")
    else:
        main()
