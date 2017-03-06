LATEX       = pdflatex
BASH        = bash -c
ECHO        = echo
RM          = rm -rf
TMP_SUFFS   = pdf aux bbl blg log dvi ps eps out toc lof lot
CHECK_RERUN =

all: ms.pdf

ms.pdf: *.tex *.bib
	${LATEX} ms
	bibtex ms 
	${LATEX} ms
	( grep Rerun ms.log && ${LATEX} ms ) || echo "Done."
	( grep Rerun ms.log && ${LATEX} ms ) || echo "Done."

clean:
	${RM} $(foreach suff, ${TMP_SUFFS}, ms.${suff})
