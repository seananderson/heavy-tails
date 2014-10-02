TEXT = anderson-etal-blackswan-timeseries
                                 

all:
	cat ms.md supp.md > temp.md
	extract_bib temp.md > refs.bib
	rm temp.md
	#extract_bib supp.md >> refs.bib
	pandoc -S --no-wrap --bibliography=refs.bib --natbib ms.md -o ms.tex
	pandoc -S --no-wrap --bibliography=refs.bib --natbib supp.md -o supp.tex
	latexmk anderson-etal-blackswan-timeseries.tex
	#pandoc --nowrap supp.md -o supp.tex
	#pdflatex

#all: $(TEXT).pdf

#$(TEXT).pdf: $(TEXT).tex
	#pdflatex $(TEXT).tex
#all:
	#pdflatex $(TEXT).tex
	#latexmk -pdf $(TEXT).tex

bib:
	pdflatex $(TEXT).tex
	bibtex $(TEXT)
	pdflatex $(TEXT).tex
	pdflatex $(TEXT).tex
	#latexmk -pdf $(TEXT).tex

clean:
	latexmk -c
