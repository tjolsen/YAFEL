texfiles = $(wildcard *tex)
output = YAFEL_overview.pdf


all: $(output)

$(output): main.pdf
	cp $< $@

main.pdf: $(texfiles)
	pdflatex main.tex
	pdflatex main.tex
	pdflatex main.tex

clean:
	rm *.aux *.log *.out *.nav *.toc *.snm *.vrb $(output) main.pdf

open: $(output)
	gnome-open $<

.PHONY: all clean open