#!/bin/bash

#   Pipeline
Rscript code/hyperframe.R
Rscript code/spat_K.R
Rscript code/spat_K_fractions.R
Rscript -e 'library(rmarkdown); rmarkdown::render("docs/03_Intraspecific_spatial.Rmd", "html_document")'
echo "Tasks complete!"