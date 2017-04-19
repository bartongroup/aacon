# New 2017
# Autogenerate the Manual pages with Sphinx v1.5.1 (using the readthedocs theme)

# I installed sphinx
pip install sphinx

# then ran the command below, choosing appropriate configurations
sphinx-quickstart

# then converted the current pages to *.rst

# then generated the html pages and pdf (via latex and local build with pdflatex)
make html
make latex

# then local build pdf using (Sublime's (MacOS app) Latex build tools)
# if pdflatex is installed the command should work
make latexpdf
