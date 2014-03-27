SHELL := /bin/bash
BUILDDIR = _build
TARGET = kidney
INPUT = "../_style:"
TEX2PDF := cd $(BUILDDIR) && TEXINPUTS=$(INPUT) pdflatex #-interaction=batchmode
BIBTEX := cd $(BUILDDIR) && BSTINPUTS=$(INPUT) bibtex

all: $(TARGET)

$(BUILDDIR):
	mkdir -p $@

$(BUILDDIR)/%.pdf: %.tex
	($(TEX2PDF) $(<F) 1>/dev/null)


%: %.tex $(BUILDDIR)
	cp $(@).{tex,bib} $(BUILDDIR)
#	($(TEX2PDF) $(@).tex)
#	($(BIBTEX) $(@))
	($(TEX2PDF) $(@).tex)
	($(TEX2PDF) $(@).tex)
	cp $(BUILDDIR)/$(@).pdf $(@).pdf

clean:
	rm -rf $(BUILDDIR)/* $(TARGET).pdf
