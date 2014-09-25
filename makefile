TEXT = anderson-etal-blackswan-timeseries
                                 
all: $(TEXT).pdf

$(TEXT).pdf: $(TEXT).tex
	latexmk -pdf $(TEXT).tex

clean:
	latexmk -c
