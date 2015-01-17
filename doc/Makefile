texfiles = $(shell ls *.tex)
mainfile = YAFEL_documentation.tex
pdffile = $(mainfile:.tex=.pdf)

all: $(pdffile)

$(pdffile): $(texfiles)	
	pdflatex $(mainfile)
	pdflatex $(mainfile)

clean:
	rm *.aux *.log *.out *.nav

open: $(pdffile)
	gnome-open $<

.PHONY: clean