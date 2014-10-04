TEXT = anderson-etal-blackswan-timeseries

all:
	cat ms.md supp.md > temp.md
	extract_bib temp.md > refs.bib
	rm temp.md
	pandoc -S --no-wrap --bibliography=refs.bib --natbib ms.md -o ms.tex
	pandoc -S --no-wrap --bibliography=refs.bib --natbib supp.md -o supp.tex
	perl -p -i -e "s/Fig. /Fig.~/g" ms.tex
	perl -p -i -e "s/Fig. /Fig.~/g" som.tex
	perl -p -i -e "s/Figs. /Figs.~/g" ms.tex
	perl -p -i -e "s/Figs. /Figs.~/g" som.tex
	pdflatex anderson-etal-blackswan-timeseries.tex

bib:
	pdflatex $(TEXT).tex
	bibtex $(TEXT)
	pdflatex $(TEXT).tex
	pdflatex $(TEXT).tex

clean:
	latexmk -c
