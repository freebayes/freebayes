all: paper

paper:
	pdflatex main
	bibtex main
	pdflatex main

clean:
	rm -f main.{aux,bbl,blg,log,pdf}

