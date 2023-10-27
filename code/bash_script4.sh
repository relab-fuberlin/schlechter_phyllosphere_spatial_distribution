#!/bin/bash

#   Pipeline
Rscript code/spat_pcf.R
Rscript code/spat_pcf_fractions.R
Rscript -e 'library(rmarkdown); rmarkdown::render("docs/results4_PCF_communitycomplexity.Rmd", "html_document")'
echo "Tasks complete!"