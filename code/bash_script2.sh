#!/bin/bash

#   Pipeline
Rscript code/coordinates.R
Rscript code/cellnumber.R
Rscript code/spat_cdftest.R
Rscript -e 'library(rmarkdown); rmarkdown::render("docs/results2_celldensity_communitycomplexity.Rmd", "html_document")'
echo "Tasks complete!"