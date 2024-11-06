#!/bin/bash

#   Input directory (data)
input=$1

#   Create an output directory
outdir=results
mkdir -p $outdir

#   Pipeline
Rscript code/data_cfu.R $input $outdir
Rscript -e 'library(rmarkdown); rmarkdown::render("docs/01_CFU.Rmd", "html_document")'
echo "Tasks complete!"