# MACOS manual — build PDF / HTML / DOCX from the Markdown source in src/.
#
#   make            -> build/macosMan.docx + build/macosMan.html
#   make pdf        -> build/macosMan.pdf   (docx -> PDF via LibreOffice)
#   make appendix-c -> regenerate src/93_appendix_c_commands.md from
#                      macos_f90/macos_help.inc
#   make clean
#
# Section order is the lexical order of src/*.md (numbered prefixes).
# Figures live in src/media/ and are referenced relative to src/.

PANDOC     ?= pandoc
SOFFICE    ?= soffice
BUILD      := build
SRC        := $(sort $(wildcard src/[0-9]*.md))
REFDOC     := reference.docx
TITLE      := MACOS Manual

PANDOC_OPTS := --resource-path=src --toc --toc-depth=3 \
               --metadata title="$(TITLE)"

CMDSRC     := $(sort $(wildcard cmdref/[0-9]*.md))
CMDREF_OPTS := --resource-path=cmdref --toc --toc-depth=2 \
               --metadata title="MACOS Command Reference"

.PHONY: all docx html pdf appendix-c clean \
        cmdref cmdref-pdf cmdref-regen

all: docx html

docx: $(BUILD)/macosMan.docx
html: $(BUILD)/macosMan.html
pdf:  $(BUILD)/macosMan.pdf

$(BUILD)/macosMan.docx: $(SRC) | $(BUILD)
	$(PANDOC) $(PANDOC_OPTS) --reference-doc=$(REFDOC) -o $@ $(SRC)

$(BUILD)/macosMan.html: $(SRC) | $(BUILD)
	$(PANDOC) $(PANDOC_OPTS) --standalone --embed-resources \
	  -o $@ $(SRC)

# No LaTeX on this box: go docx -> PDF through headless LibreOffice.
# tools/docx2pdf.py refreshes the TOC field (page numbers) before export;
# a plain `soffice --convert-to pdf` would leave the TOC empty.
$(BUILD)/macosMan.pdf: $(BUILD)/macosMan.docx
	python3 tools/docx2pdf.py $< $@

appendix-c:
	python3 tools/gen_appendix_c.py

# ---- MACOS Command Reference (cmdref/) ------------------------------
cmdref: $(BUILD)/macosCmdRef.docx $(BUILD)/macosCmdRef.html
cmdref-pdf: $(BUILD)/macosCmdRef.pdf

# regenerate the auto-skeleton catalogs from engine + binding sources
cmdref-regen:
	python3 tools/gen_cmdref.py

$(BUILD)/macosCmdRef.docx: $(CMDSRC) | $(BUILD)
	$(PANDOC) $(CMDREF_OPTS) --reference-doc=$(REFDOC) -o $@ $(CMDSRC)

$(BUILD)/macosCmdRef.html: $(CMDSRC) | $(BUILD)
	$(PANDOC) $(CMDREF_OPTS) --standalone --embed-resources \
	  -o $@ $(CMDSRC)

$(BUILD)/macosCmdRef.pdf: $(BUILD)/macosCmdRef.docx
	python3 tools/docx2pdf.py $< $@

$(BUILD):
	mkdir -p $(BUILD)

clean:
	rm -rf $(BUILD)
