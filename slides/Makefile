FILE  = talk
THEME = Copenhagen
BUILD = _build
#FLAGS = --theme=$(THEME) --stylesheet=preamble.txt --centerfigs=True --overlaybullets=False
FLAGS = --theme=$(THEME) --stylesheet=preamble.txt --overlaybullets=False

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  latex        to make LaTeX files"
	@echo "  latex-notes  to make LaTeX files with notes"
	@echo "  pdf          to make PDF slides"
	@echo "  notes        to make PDF slides with notes"
	@echo "  clean        to delete generated files"

latex:
	cd _build && rst2beamer $(FLAGS) ../$(FILE).rst $(FILE).tex
	cd _build && rst2beamer $(FLAGS) ../$(FILE).rst $(FILE).tex
	@echo
	@echo "Build finished. The TEX file is in _build."

latex-notes:
	cd _build && rst2beamer $(FLAGS) --shownotes bottom ../$(FILE).rst $(FILE)-notes.tex
	cd _build && rst2beamer $(FLAGS) --shownotes bottom ../$(FILE).rst $(FILE)-notes.tex
	@echo
	@echo "Build finished. The TEX file is in _build."

pdf: latex
	cd _build && pdflatex $(FILE).tex
	@echo
	@echo "Build finished. The PDF is _build/$(FILE).pdf."

notes: latex-notes
	cd _build && pdflatex $(FILE)-notes.tex
	@echo
	@echo "Build finished. The PDF is _build/$(FILE)-notes.pdf."

clean:
	rm $(BUILD)/$(FILE){,-notes}.{aux,log,nav,out,pdf,snm,tex,toc,vrb}
