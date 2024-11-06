#!/bin/bash

#   Pipeline
Rscript code/coordinates.R
Rscript code/cellnumber.R
Rscript code/spat_cdftest.R
Rscript -e 'library(rmarkdown); rmarkdown::render("docs/02_celldensity.Rmd", "html_document")'
echo "Tasks complete!"