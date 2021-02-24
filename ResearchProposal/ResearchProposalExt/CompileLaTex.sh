#!/bin/bash
# Author: Amy Solman amy.solman@imperial.ac.uk
# Script: CompileLaTex.sh
# Desc: Compiles LaTex document with citation and produces pdf file
# Date: Oct 2019
pdflatex $1.tex
pdflatex $1.tex
bibtex $1
pdflatex $1.tex
pdflatex $1.tex
evince $1.pdf &

## Cleanup
rm *~
rm *.aux
rm *.dvi
rm *.log
rm *.nav
rm *.out
rm *.out
rm *.snm
rm *.toc
rm *.blg
rm *.bbl