#!/usr/bin/env bash

# run from ./docs
rm -rf _build
# generating the html pages
make html
mkdir ../../website/docs/
cp -r _build/html/* ../../website/docs/
rm -rf _build

# # generating the pdf
# make latex
# # manual edit to remove unecessary ines in the index page

## add the following to the latex and change the date to match the release date

# \author{Agnieszka Golicz \and Peter V. Troshin \and Fábio Madeira \and David M. A. Martin \and James B. Procter \and Geoffrey J. Barton}
# \authoraddress{
# The Barton Group
#
# Division of Computational Biology
#
# School of Life Sciences
#
# University of Dundee
#
# Dow Street
#
# Dundee DD1 5EH
#
# Scotland, UK
# }

# add the note.
# # \begin{sphinxadmonition}{note}{Note:}
# # Some of the hyper-links found in this pdf might not work properly.
# # \end{sphinxadmonition}

# cp _build/latex/aacon.pdf ../../website/docs/aacon_manual.pdf
# rm -rf _build

# for multiple versions of this manual we need to implement some solution
# like this one https://robpol86.github.io/sphinxcontrib-versioning/index.html
# in the future
