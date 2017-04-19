#!/usr/bin/env bash

# run from ./docs
rm -rf _build
# generating the html pages
make html
mkdir ../website/docs
cp -r _build/html/* ../website/docs/
rm -rf _build

# # generating the pdf
# make latex
# # manual edit to remove unecessary ines in the index page
# # \begin{sphinxadmonition}{note}{Note:}
# # Some of the links might not work properly in the pdf version...
# # \end{sphinxadmonition}
# cp _build/latex/AACON.pdf ../website/docs/aacon_manual.pdf

# for multiple versions of this manual we need to implement some solution
# like this one https://robpol86.github.io/sphinxcontrib-versioning/index.html
# in the future
