#!/bin/bash

#   Pipeline
Rscript code/spat_pcf.R
Rscript code/spat_pcf_fractions.R
Rscript -e 'library(rmarkdown); rmarkdown::render("docs/04_Interspecific_spatial.Rmd", "html_document")'
echo "Tasks complete!"