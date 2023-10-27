#!/bin/bash

#   Pipeline
Rscript code/hyperframe.R
Rscript code/spat_K.R
Rscript code/spat_K_fractions.R
Rscript -e 'library(rmarkdown); rmarkdown::render("docs/results3_Kest_communitycomplexity.Rmd", "html_document")'
echo "Tasks complete!"