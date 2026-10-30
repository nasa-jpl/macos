#!/usr/bin/env python3
"""Build reference.docx (the pandoc --reference-doc for both the manual
and the command reference) deterministically from pandoc's default.

Why: pandoc's default reference docx uses THEME font references
(asciiTheme="majorHAnsi" etc.) and leaves body styles with no explicit
font; different Word/LibreOffice versions resolve themes differently,
so output styling varied by machine.  It also lacks a "Source Code"
paragraph style, so code blocks rendered unpredictably.  This script
pins EXPLICIT fonts on every style we use and injects a proper
SourceCode style.  Fonts are chosen to substitute cleanly on Linux
(Liberation family) while native on Windows/Mac.

Usage (from docs/macos-manual):  python3 tools/make_reference_docx.py
Rebuilds ./reference.docx in place.
"""
import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "reference.docx"

HEAD_FONT = "Arial"
BODY_FONT = "Times New Roman"
MONO_FONT = "Courier New"

EXPLICIT = {
    # styleId: (font, halfpoint-size or None)
    "Normal":         (BODY_FONT, 22),
    "BodyText":       (BODY_FONT, None),
    "FirstParagraph": (BODY_FONT, None),
    "Compact":        (BODY_FONT, None),
    "Title":          (HEAD_FONT, 56),
    "Subtitle":       (HEAD_FONT, None),
    "Heading1":       (HEAD_FONT, 32),
    "Heading2":       (HEAD_FONT, 28),
    "Heading3":       (HEAD_FONT, 24),
    "Heading4":       (HEAD_FONT, 22),
    "Heading5":       (HEAD_FONT, None),
    "BlockText":      (BODY_FONT, None),
    "TOCHeading":     (HEAD_FONT, None),
}

RFONTS = ('<w:rFonts w:ascii="{f}" w:hAnsi="{f}" w:eastAsia="{f}"'
          ' w:cs="{f}" />')

SOURCECODE_STYLE = f"""<w:style w:type="paragraph" w:styleId="SourceCode">
<w:name w:val="Source Code" />
<w:basedOn w:val="Normal" />
<w:pPr>
<w:wordWrap w:val="0" />
<w:spacing w:before="60" w:after="60" w:line="240" w:lineRule="auto" />
<w:ind w:left="240" />
<w:shd w:val="clear" w:color="auto" w:fill="F5F5F5" />
</w:pPr>
<w:rPr>
<w:rFonts w:ascii="{MONO_FONT}" w:hAnsi="{MONO_FONT}" w:cs="{MONO_FONT}" />
<w:sz w:val="15" />
<w:szCs w:val="15" />
</w:rPr>
</w:style>"""

# centered figure paragraphs (tools/figure_captions.lua wraps the image and the FIGURE line)
FIGURE_STYLES = f"""<w:style w:type="paragraph" w:styleId="Figure">
<w:name w:val="Figure" />
<w:basedOn w:val="Normal" />
<w:pPr><w:keepNext /><w:spacing w:before="120" w:after="60" /><w:jc w:val="center" /></w:pPr>
</w:style>
<w:style w:type="paragraph" w:styleId="FigureCaption">
<w:name w:val="FigureCaption" />
<w:basedOn w:val="Normal" />
<w:pPr><w:spacing w:before="60" w:after="240" /><w:jc w:val="center" /></w:pPr>
<w:rPr><w:rFonts w:ascii="{BODY_FONT}" w:hAnsi="{BODY_FONT}" w:cs="{BODY_FONT}" /><w:sz w:val="20" /><w:szCs w:val="20" /></w:rPr>
</w:style>"""


def patch_style(xml, sid, font, size):
    m = re.search(r'<w:style [^>]*w:styleId="%s".*?</w:style>' % sid,
                  xml, re.S)
    if not m:
        print(f"  {sid:16s} not in reference — skipped")
        return xml
    block = m.group(0)
    nb = block
    rf = RFONTS.format(f=font)
    if "<w:rFonts" in nb:
        nb = re.sub(r"<w:rFonts[^/>]*/?>", rf, nb, count=1)
    elif "<w:rPr>" in nb:
        nb = nb.replace("<w:rPr>", "<w:rPr>" + rf, 1)
    else:
        nb = nb.replace("</w:style>",
                        "<w:rPr>" + rf + "</w:rPr></w:style>", 1)
    if size is not None:
        if "<w:sz " in nb:
            nb = re.sub(r'<w:sz w:val="\d+"\s*/>',
                        f'<w:sz w:val="{size}" />', nb, count=1)
        else:
            nb = nb.replace(rf, rf + f'<w:sz w:val="{size}" />', 1)
    print(f"  {sid:16s} -> {font}" + (f" {size/2:.0f}pt" if size else ""))
    return xml.replace(block, nb)


def main():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        base = td / "base.docx"
        base.write_bytes(subprocess.run(
            ["pandoc", "--print-default-data-file", "reference.docx"],
            capture_output=True, check=True).stdout)
        work = td / "unpacked"
        with zipfile.ZipFile(base) as z:
            z.extractall(work)
        sf = work / "word" / "styles.xml"
        xml = sf.read_text(encoding="utf8")

        for sid, (font, size) in EXPLICIT.items():
            xml = patch_style(xml, sid, font, size)

        # character styles used for inline/verbatim code
        # inline code: mono at 9pt (it inherited the 11pt body); code
        # BLOCKS (SourceCode below) at 7.5pt, indent 240 (a 100-column transcript line then fits the 6.5 in text width) -- Dave 2026-10-19, the
        # transcripts were too large beside the 11pt prose
        for sid in ("VerbatimChar",):
            xml = patch_style(xml, sid, MONO_FONT, 15)

        # inject SourceCode paragraph style if absent
        if 'w:styleId="SourceCode"' not in xml:
            xml = xml.replace("</w:styles>",
                              SOURCECODE_STYLE + "</w:styles>")
            print("  SourceCode       injected (mono, 7.5pt, shaded)")

        # page numbers: a footer part with a centered PAGE field, referenced from the
        # document's section properties; the count starts at 0 so the cover is page 0
        # (Dave 2026-10-19).  pandoc keeps the reference document's sectPr.
        FOOTER = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                  '<w:ftr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
                  'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
                  '<w:p><w:pPr><w:jc w:val="center"/></w:pPr>'
                  '<w:r><w:fldChar w:fldCharType="begin"/></w:r>'
                  '<w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>'
                  '<w:r><w:fldChar w:fldCharType="separate"/></w:r>'
                  '<w:r><w:t>0</w:t></w:r>'
                  '<w:r><w:fldChar w:fldCharType="end"/></w:r></w:p></w:ftr>')
        (work / "word" / "footer1.xml").write_text(FOOTER, encoding="utf8")
        rels = work / "word" / "_rels" / "document.xml.rels"
        rx = rels.read_text(encoding="utf8")
        if "footer1.xml" not in rx:
            rx = rx.replace("</Relationships>",
                            '<Relationship Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/footer" '
                            'Id="rIdFooter1" Target="footer1.xml" /></Relationships>')
            rels.write_text(rx, encoding="utf8")
        ct = work / "[Content_Types].xml"
        cx = ct.read_text(encoding="utf8")
        if "footer1.xml" not in cx:
            cx = cx.replace("</Types>",
                            '<Override PartName="/word/footer1.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml" /></Types>')
            ct.write_text(cx, encoding="utf8")
        doc = work / "word" / "document.xml"
        dx = doc.read_text(encoding="utf8")
        SECT = ('<w:sectPr><w:footerReference w:type="default" r:id="rIdFooter1"/>'
                '<w:pgSz w:w="12240" w:h="15840"/>'
                '<w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="720" w:footer="720" w:gutter="0"/>'
                '<w:pgNumType w:start="0"/></w:sectPr>')
        import re as _re
        dx2, n = _re.subn(r"<w:sectPr\s*/>|<w:sectPr>.*?</w:sectPr>", SECT, dx, count=1, flags=_re.S)
        assert n == 1, "no sectPr in the reference document"
        doc.write_text(dx2, encoding="utf8")
        print("  footer1 (PAGE field, centered) + sectPr (letter, 1in, numbering from 0) injected")

        if 'w:styleId="FigureCaption"' not in xml:
            xml = xml.replace("</w:styles>", FIGURE_STYLES + "</w:styles>")
            print("  Figure/FigureCaption injected (centered)")

        sf.write_text(xml, encoding="utf8")

        # repack (styles.xml first is not required; zip whole tree)
        if OUT.exists():
            OUT.unlink()
        with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as z:
            for f in sorted(work.rglob("*")):
                if f.is_file():
                    z.write(f, f.relative_to(work))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
