TEXT = anderson-etal-blackswan-timeseries
MS = ms
SOM = supp
VERSION = 4.0

all: text

text:
	latexmk $(TEXT)

quick:
	pdflatex $(TEXT)

dropbox: text rtf
	cp $(TEXT).pdf ~/Dropbox/heavy-tails/$(TEXT)-v$(VERSION).pdf
	mv $(TEXT).rtf ~/Dropbox/heavy-tails/$(TEXT)-v$(VERSION).rtf
	open ~/Dropbox/heavy-tails/

rtf: text
	latex2rtf -E0 $(TEXT).tex

extractbib: text
	bibtool -x $(MS).aux -o $(MS).bib -- 'expand.macros = ON'
	bibtool -x $(SOM).aux -o $(SOM).bib -- 'expand.macros = ON'

clean:
	latexmk -c
