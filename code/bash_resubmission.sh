#!/bin/bash

#   Create an output directory
mkdir -p results

#   Pipeline
Rscript -e 'library(rmarkdown); rmarkdown::render("docs/05resubmission.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("docs/06resubmission.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("docs/07resubmission.Rmd", "html_document")'