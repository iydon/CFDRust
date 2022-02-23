.PHONY: clean

main.pdf: main.tex config/* part/*
	latexmk -xelatex -shell-escape main

clean:
	latexmk -C
	rm main.bbl main.listing main.run.xml
