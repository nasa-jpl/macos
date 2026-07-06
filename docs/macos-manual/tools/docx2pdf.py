#!/usr/bin/env python3
"""docx -> PDF via headless LibreOffice, refreshing the Table of
Contents field first (a plain `soffice --convert-to pdf` exports the
TOC unpopulated).

Usage: python3 tools/docx2pdf.py input.docx output.pdf
"""
import os
import subprocess
import sys
import time
import uno  # noqa: F401  (ships with LibreOffice; import must succeed)
from com.sun.star.beans import PropertyValue


def prop(name, value):
    p = PropertyValue()
    p.Name = name
    p.Value = value
    return p


def main(docx, pdf):
    port = 21713
    soffice = subprocess.Popen(
        ["soffice", "--headless", "--invisible", "--norestore",
         f"--accept=socket,host=localhost,port={port};urp;",
         "--nologo", "--nodefault"],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        local = uno.getComponentContext()
        resolver = local.ServiceManager.createInstanceWithContext(
            "com.sun.star.bridge.UnoUrlResolver", local)
        ctx = None
        for _ in range(30):
            try:
                ctx = resolver.resolve(
                    f"uno:socket,host=localhost,port={port};urp;"
                    "StarOffice.ComponentContext")
                break
            except Exception:
                time.sleep(0.5)
        if ctx is None:
            raise RuntimeError("could not connect to soffice listener")

        desktop = ctx.ServiceManager.createInstanceWithContext(
            "com.sun.star.frame.Desktop", ctx)
        url_in = uno.systemPathToFileUrl(os.path.abspath(docx))
        doc = desktop.loadComponentFromURL(
            url_in, "_blank", 0, (prop("Hidden", True),))

        # refresh all document indexes (TOC) twice: first pass paginates,
        # second pass fixes page numbers shifted by the TOC's own length
        for _ in range(2):
            indexes = doc.getDocumentIndexes()
            for i in range(indexes.getCount()):
                indexes.getByIndex(i).update()
            doc.refresh()

        url_out = uno.systemPathToFileUrl(os.path.abspath(pdf))
        doc.storeToURL(url_out, (
            prop("FilterName", "writer_pdf_Export"),))
        doc.close(False)
        print(f"wrote {pdf}")
    finally:
        soffice.terminate()
        soffice.wait()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
