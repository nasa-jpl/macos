# mmacos manual -- MATLAB + MACOS with an AI assistant (PLAN_MMACOS_MANUAL.md, BRIEF_to_mmacos_manual.md).
#
#   make             -> build/mmacosMan.docx + build/mmacosMan.html
#   make pdf         -> build/mmacosMan.pdf   (docx -> PDF via LibreOffice)
#   make transcripts -> regenerate src/transcripts/ and src/media/ from the chapter drivers
#                       tools/transcripts/<chapter>.m (tools/regen_mmacos_transcripts.py; one MATLAB
#                       -batch run per driver; CH=<name> runs one)
#   make clean
#
# The pipeline is the MACOS manual's, reused BY REFERENCE (MANTOOLS): the Lua filters (transcript
# markers, section numbering + cross-references, figure captions), reference.docx, the CSS and the
# docx -> PDF step.  Only the transcript generator is this manual's own (MATLAB drivers, not CLI
# journals).  Section order is the lexical order of src/*.md; figures live in src/media/.

PANDOC     ?= pandoc
MANUAL     := ../macos-manual
MANTOOLS   := $(MANUAL)/tools
REFDOC     := $(MANUAL)/reference.docx
BUILD      := build
SRC        := $(sort $(wildcard src/[0-9]*.md))
TRANSCRIPTS := $(wildcard src/transcripts/*.md)
MEDIA      := $(wildcard src/media/*)
TITLE      := mmacos Manual

PANDOC_OPTS := --resource-path=src --toc --toc-depth=3 --metadata title="$(TITLE)"
DOCX_OPTS  := --resource-path=src

# tools/pins.lua (this manual's own): %pin:NAME% in the prose -> the number the chapter's driver produced
NUMBER     := --lua-filter=$(MANTOOLS)/include_transcripts.lua --lua-filter=tools/pins.lua \
              --lua-filter=$(MANTOOLS)/number_sections.lua --lua-filter=$(MANTOOLS)/figure_captions.lua
FILTERS    := $(MANTOOLS)/include_transcripts.lua tools/pins.lua $(MANTOOLS)/number_sections.lua \
              $(MANTOOLS)/figure_captions.lua src/transcripts/pins.tsv

.PHONY: all docx html pdf transcripts clean

all: docx html

docx: $(BUILD)/mmacosMan.docx
html: $(BUILD)/mmacosMan.html
pdf:  $(BUILD)/mmacosMan.pdf

$(BUILD)/mmacosMan.docx: $(SRC) $(TRANSCRIPTS) $(MEDIA) $(REFDOC) $(FILTERS) | $(BUILD)
	$(PANDOC) $(DOCX_OPTS) $(NUMBER) --reference-doc=$(REFDOC) -o $@ $(SRC)

$(BUILD)/mmacosMan.html: $(SRC) $(TRANSCRIPTS) $(MEDIA) $(FILTERS) $(MANTOOLS)/manual.css | $(BUILD)
	$(PANDOC) $(PANDOC_OPTS) $(NUMBER) --standalone --embed-resources --mathml \
	  --css $(MANTOOLS)/manual.css -o $@ $(SRC)

$(BUILD)/mmacosMan.pdf: $(BUILD)/mmacosMan.docx
	python3 $(MANTOOLS)/docx2pdf.py $< $@

transcripts:
	python3 tools/regen_mmacos_transcripts.py $(CH)

$(BUILD):
	mkdir -p $@

clean:
	rm -rf $(BUILD)
