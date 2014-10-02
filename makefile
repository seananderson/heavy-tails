TEXT = anderson-etal-blackswan-timeseries
                                 
#all: $(TEXT).pdf

#$(TEXT).pdf: $(TEXT).tex
	#pdflatex $(TEXT).tex
all:
	pdflatex $(TEXT).tex
	#latexmk -pdf $(TEXT).tex

bib:
	pdflatex $(TEXT).tex
	bibtex $(TEXT)
	pdflatex $(TEXT).tex
	pdflatex $(TEXT).tex
	#latexmk -pdf $(TEXT).tex

clean:
	latexmk -c
