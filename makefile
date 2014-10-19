TEXT = anderson-etal-blackswan-timeseries
MS = ms
SOM = supp

all: quick

quick:
	pdflatex $(TEXT)

text: 
	latexmk $(TEXT)

dropbox: text
	cp $(TEXT).pdf ~/Dropbox/Public/$(TEXT)-v3.pdf

rtf: text
	latex2rtf -E0 anderson-etal-blackswan-timeseries.tex

extractbib: text
	bibtool -x $(MS).aux -o $(MS).bib -- 'expand.macros = ON'
	bibtool -x $(SOM).aux -o $(SOM).bib -- 'expand.macros = ON'

clean:
	latexmk -c
