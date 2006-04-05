all:	vign

vign:	tr_2006_01.tex
	pdflatex tr_2006_01
	pdflatex tr_2006_01
	make clean

pdf:	tr_2006_01.tex
	echo "Sweave('tr_2006_01.Rnw')" | nohup R-rel --vanilla > out.txt
	pdflatex tr_2006_01
	pdflatex tr_2006_01
	make clean
	acroread tr_2006_01.pdf &

clean:
	rm -f *.out *.bbl *.log *.aux *.blg *.toc  *~
	rm -f *.dvi out.txt tr_2006_01-* Rplots.ps
