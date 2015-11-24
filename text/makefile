TEXT = anderson-etal-blackswan-timeseries
MS = ms
SOM = supp
VERSION = 6.0

all: text

text:
# 	pandoc -S --no-wrap --bibliography=$(REFS).bib --natbib -fmarkdown-all_symbols_escapable ms.md -o ms.tex
# 	perl -p -i -e "s/\\\citep/\\\cite/g" ms.tex
# 	perl -p -i -e "s/\\\citet/\\\cite/g" ms.tex
# 	perl -p -i -e "s/Figure /Figure~/g" ms.tex
# 	perl -p -i -e "s/Fig\. /Fig.~/g" ms.tex
# 	perl -p -i -e "s/vs\. /vs.~/g" ms.tex
# 	perl -p -i -e "s/Table /Table~/g" ms.tex
# 	perl -p -i -e "s/Figs\. /Figs.~/g" ms.tex
	perl -p -i -e "s/\\\citep/\\\cite/g" ../analysis/cause-table.tex
# 	perl -0777p -i -e "s/\n\\\begin{align}/\\\begin{align}/gm" ms.tex
# 	perl -0777p -i -e "s/\\\end{align}\n/\\\end{align}/gm" ms.tex
	perl -p -i -e "s/\\\citet/\\\cite/g" ../analysis/cause-table.tex
	latexmk $(TEXT)
	#pdflatex $(TEXT)
	latexmk $(SOM)

supp:
	latexmk -pdf $(SOM)

quick:
	pdflatex $(TEXT)

dropbox: text
	cp $(TEXT).pdf ~/Dropbox/heavy-tails/$(TEXT)-v$(VERSION).pdf
	mv $(TEXT).rtf ~/Dropbox/heavy-tails/$(TEXT)-v$(VERSION).rtf
	cp $(SOM).pdf ~/Dropbox/heavy-tails/$(SOM)-v$(VERSION).pdf
	cp $(MS).md ~/Dropbox/heavy-tails/$(TEXT)-v$(VERSION).md
	open ~/Dropbox/heavy-tails/

rtf: text
	latex2rtf -M4 -t2 -E0 $(TEXT).tex

extractbib: text
	bibtool -x $(MS).aux -o $(MS).bib -- 'expand.macros = ON'
	bibtool -x $(SOM).aux -o $(SOM).bib -- 'expand.macros = ON'

clean:
	latexmk -c

#continuous:
	#latexmk -pvc $(TEXT)

continuous:
	while true; do make --silent; sleep 1; done
