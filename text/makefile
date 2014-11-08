TEXT = anderson-etal-blackswan-timeseries
MS = ms
SOM = supp
VERSION = 3.2

all: text

text:
	latexmk $(TEXT)

quick:
	pdflatex $(TEXT)

dropbox: text rtf
	#rm ~/Dropbox/heavy-tails/*
	cp $(TEXT).pdf ~/Dropbox/heavy-tails/$(TEXT)-v$(VERSION).pdf
	mv $(TEXT).rtf ~/Dropbox/heavy-tails/$(TEXT)-v$(VERSION).rtf
	open ~/Dropbox/heavy-tails/

rtf: text
	latex2rtf -E0 anderson-etal-blackswan-timeseries.tex

extractbib: text
	bibtool -x $(MS).aux -o $(MS).bib -- 'expand.macros = ON'
	bibtool -x $(SOM).aux -o $(SOM).bib -- 'expand.macros = ON'

clean:
	latexmk -c
