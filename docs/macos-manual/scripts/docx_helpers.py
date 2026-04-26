"""
docx_helpers.py — Shared utilities for reading and editing document.xml.

Import this in any transform script:
    from scripts.docx_helpers import load_doc, save_doc, get_style, set_style,
                                      get_text, get_run_info, NAMESPACES, W, wtag
"""

import xml.etree.ElementTree as ET
from pathlib import Path
from collections import Counter
from typing import Optional

# ── Namespace registration ────────────────────────────────────────────────────
NAMESPACES = {
    "wpc": "http://schemas.microsoft.com/office/word/2010/wordprocessingCanvas",
    "r":   "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "m":   "http://schemas.openxmlformats.org/officeDocument/2006/math",
    "v":   "urn:schemas-microsoft-com:vml",
    "wp":  "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing",
    "w":   "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
    "w14": "http://schemas.microsoft.com/office/word/2010/wordml",
    "w15": "http://schemas.microsoft.com/office/word/2012/wordml",
    "wne": "http://schemas.microsoft.com/office/word/2006/wordml",
    "mc":  "http://schemas.openxmlformats.org/markup-compatibility/2006",
    "o":   "urn:schemas-microsoft-com:office:office",
    "wpg": "http://schemas.microsoft.com/office/word/2010/wordprocessingGroup",
    "wpi": "http://schemas.microsoft.com/office/word/2010/wordprocessingInk",
    "wps": "http://schemas.microsoft.com/office/word/2010/wordprocessingShape",
}
for prefix, uri in NAMESPACES.items():
    ET.register_namespace(prefix, uri)

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"

def wtag(name: str) -> str:
    """Return a fully-qualified w: tag name, e.g. wtag('p') → '{...}p'."""
    return f"{{{W}}}{name}"


# ── Load / save ───────────────────────────────────────────────────────────────

def load_doc(working_dir: Path) -> tuple[ET.ElementTree, ET.Element]:
    """Parse working/word/document.xml. Returns (tree, root)."""
    path = Path(working_dir) / "word" / "document.xml"
    tree = ET.parse(path)
    return tree, tree.getroot()


def save_doc(tree: ET.ElementTree, working_dir: Path) -> None:
    """Write tree back to working/word/document.xml."""
    path = Path(working_dir) / "word" / "document.xml"
    tree.write(str(path), xml_declaration=True, encoding="UTF-8")


def load_styles(working_dir: Path) -> tuple[ET.ElementTree, ET.Element]:
    """Parse working/word/styles.xml. Returns (tree, root)."""
    path = Path(working_dir) / "word" / "styles.xml"
    tree = ET.parse(path)
    return tree, tree.getroot()


def save_styles(tree: ET.ElementTree, working_dir: Path) -> None:
    path = Path(working_dir) / "word" / "styles.xml"
    tree.write(str(path), xml_declaration=True, encoding="UTF-8")


# ── Paragraph helpers ─────────────────────────────────────────────────────────

def get_style(para: ET.Element) -> str:
    """Return the paragraph style ID, or 'Normal' if unset."""
    pPr = para.find(wtag("pPr"))
    if pPr is not None:
        pStyle = pPr.find(wtag("pStyle"))
        if pStyle is not None:
            return pStyle.get(wtag("val"), "Normal")
    return "Normal"


def set_style(para: ET.Element, style_id: str) -> None:
    """Set (or replace) the paragraph style ID."""
    pPr = para.find(wtag("pPr"))
    if pPr is None:
        pPr = ET.Element(wtag("pPr"))
        para.insert(0, pPr)
    pStyle = pPr.find(wtag("pStyle"))
    if pStyle is None:
        pStyle = ET.SubElement(pPr, wtag("pStyle"))
        pPr.insert(0, pStyle)
    pStyle.set(wtag("val"), style_id)


def get_text(para: ET.Element) -> str:
    """Return all text content of a paragraph, concatenated."""
    return "".join(t.text or "" for t in para.iter(wtag("t")))


def get_run_info(para: ET.Element) -> dict:
    """
    Return formatting of the first content run as a dict:
        { font, size, bold, italic, mono }
    size is a string in half-points (e.g. '18' = 9pt, '24' = 12pt).
    """
    for r in para.iter(wtag("r")):
        t = r.find(wtag("t"))
        if t is not None and (t.text or "").strip():
            rPr = r.find(wtag("rPr"))
            info: dict = {"font": "", "size": "", "bold": False,
                          "italic": False, "mono": False}
            if rPr is not None:
                f = rPr.find(wtag("rFonts"))
                if f is not None:
                    info["font"] = (f.get(wtag("ascii"))
                                    or f.get(wtag("hAnsi")) or "")
                s = rPr.find(wtag("sz"))
                if s is not None:
                    info["size"] = s.get(wtag("val"), "")
                info["bold"]   = rPr.find(wtag("b"))  is not None
                info["italic"] = rPr.find(wtag("i"))  is not None
                info["mono"]   = any(k in info["font"]
                                     for k in ("Courier", "Mono", "mono",
                                               "Consolas", "Monaco"))
            return info
    return {"font": "", "size": "", "bold": False, "italic": False, "mono": False}


def has_image(para: ET.Element) -> bool:
    """Return True if the paragraph contains an inline image/drawing."""
    return bool(para.findall(f".//{wtag('drawing')}"))


# ── Document-level inspection ─────────────────────────────────────────────────

def style_counts(root: ET.Element) -> Counter:
    """Return a Counter of paragraph style IDs across the whole document."""
    return Counter(get_style(p) for p in root.iter(wtag("p")))


def iter_paragraphs(root: ET.Element):
    """Yield every <w:p> in document order (including inside tables)."""
    return root.iter(wtag("p"))


def print_style_summary(root: ET.Element) -> None:
    counts = style_counts(root)
    print("Paragraph style counts:")
    for style, count in counts.most_common():
        print(f"  {count:5d}  {style}")
    print(f"  -----")
    print(f"  {sum(counts.values()):5d}  total")


# ── Paragraph creation / mutation ─────────────────────────────────────────────

XML_SPACE = "{http://www.w3.org/XML/1998/namespace}space"


def make_para(style: str, text: str = "") -> ET.Element:
    """
    Create a fresh <w:p> with given paragraph style and (optional) text.
    Splits on '\\n' so callers can pass multi-line strings; line breaks
    become <w:br/> within a single paragraph.
    """
    p = ET.Element(wtag("p"))
    pPr = ET.SubElement(p, wtag("pPr"))
    pStyle = ET.SubElement(pPr, wtag("pStyle"))
    pStyle.set(wtag("val"), style)
    if text:
        r = ET.SubElement(p, wtag("r"))
        for i, line in enumerate(text.split("\n")):
            if i > 0:
                ET.SubElement(r, wtag("br"))
            t = ET.SubElement(r, wtag("t"))
            t.set(XML_SPACE, "preserve")
            t.text = line
    return p


def find_para_by_text(root: ET.Element,
                      substr: str,
                      style: Optional[str] = None,
                      start: int = 0) -> Optional[ET.Element]:
    """
    Return the first <w:p> whose concatenated text contains `substr`
    (case-sensitive). If `style` is given, restrict to paragraphs of that
    style. `start` skips that many matching candidates first.
    """
    skipped = 0
    for p in root.iter(wtag("p")):
        if style and get_style(p) != style:
            continue
        if substr in get_text(p):
            if skipped < start:
                skipped += 1
                continue
            return p
    return None


def find_para_parent(root: ET.Element,
                     target: ET.Element) -> Optional[ET.Element]:
    """Return the parent element of `target` somewhere under `root`."""
    for parent in root.iter():
        for child in list(parent):
            if child is target:
                return parent
    return None


def insert_after_para(root: ET.Element,
                      anchor: ET.Element,
                      new_paras: list) -> None:
    """Insert each <w:p> in new_paras directly after anchor in its parent."""
    parent = find_para_parent(root, anchor)
    if parent is None:
        raise RuntimeError("anchor paragraph not found in tree")
    children = list(parent)
    idx = children.index(anchor)
    for offset, np in enumerate(new_paras, start=1):
        parent.insert(idx + offset, np)


def replace_text_in_para(para: ET.Element, old: str, new: str) -> int:
    """Replace `old` with `new` in every <w:t> of para. Returns number of hits."""
    count = 0
    for t in para.iter(wtag("t")):
        if t.text and old in t.text:
            t.text = t.text.replace(old, new)
            count += 1
    return count


def replace_para_text(para: ET.Element, new_text: str, style: Optional[str] = None) -> None:
    """
    Drop all runs in `para` and replace with a single run holding `new_text`.
    Preserves the paragraph's pPr (style/indent/etc.) unless `style` is given,
    in which case the style is also rewritten.
    """
    pPr = para.find(wtag("pPr"))
    # Remove every child except pPr
    for child in list(para):
        if child is not pPr:
            para.remove(child)
    if style is not None:
        if pPr is None:
            pPr = ET.Element(wtag("pPr"))
            para.insert(0, pPr)
        pStyle = pPr.find(wtag("pStyle"))
        if pStyle is None:
            pStyle = ET.SubElement(pPr, wtag("pStyle"))
        pStyle.set(wtag("val"), style)
    r = ET.SubElement(para, wtag("r"))
    t = ET.SubElement(r, wtag("t"))
    t.set(XML_SPACE, "preserve")
    t.text = new_text
