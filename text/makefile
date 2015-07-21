TEXT = anderson-etal-blackswan-timeseries
MS = ms
SOM = supp
VERSION = 4.1

all: text

text:
	pandoc -S --no-wrap --bibliography=$(REFS).bib --natbib -fmarkdown-all_symbols_escapable ms.md -o ms.tex
	perl -p -i -e "s/\\\citep/\\\cite/g" ms.tex
	perl -p -i -e "s/\\\citet/\\\cite/g" ms.tex
	perl -p -i -e "s/Figure /Figure~/g" ms.tex
	perl -p -i -e "s/Fig\. /Fig.~/g" ms.tex
	perl -p -i -e "s/vs\. /vs.~/g" ms.tex
	perl -p -i -e "s/Table /Table~/g" ms.tex
	perl -p -i -e "s/Figs\. /Figs.~/g" ms.tex
	perl -p -i -e "s/\\\citep/\\\cite/g" ../analysis/cause-table.tex
	perl -p -i -e "s/\\\citet/\\\cite/g" ../analysis/cause-table.tex
	latexmk $(TEXT)
	#pdflatex $(TEXT)
	latexmk $(SOM)

supp:
	latexmk $(SOM)

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

continuous:
	latexmk -pvc $(TEXT)
