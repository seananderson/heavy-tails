TEXT = anderson-etal-blackswan-timeseries

all:
	cat ms.md supp.md > temp.md
	extract_bib temp.md > refs.bib
	rm temp.md
	pandoc -S --no-wrap --bibliography=refs.bib --natbib ms.md -o ms.tex
	pandoc -S --no-wrap --bibliography=refs.bib --natbib supp.md -o supp.tex
	perl -p -i -e "s/Fig. /Fig.~/g" ms.tex
	perl -p -i -e "s/Fig. /Fig.~/g" supp.tex
	perl -p -i -e "s/Figs. /Figs.~/g" ms.tex
	perl -p -i -e "s/Figs. /Figs.~/g" supp.tex
	pdflatex anderson-etal-blackswan-timeseries.tex
	cp anderson-etal-blackswan-timeseries.pdf ~/Dropbox/Public/anderson-etal-blackswan-timeseries-v2.pdf

bib:
	pdflatex $(TEXT)
	bibtex ms
	bibtex supp
	pdflatex $(TEXT)
	pdflatex $(TEXT)

clean:
	latexmk -c
