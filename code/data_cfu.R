#!/usr/bin/env Rscript

# Data processing CFU data
# Load the required libraries
library(here)
library(tidyverse)

# Check command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("You must provide three command line arguments: [1] input directory and 
         [2] output directory")
}

# Check the output directory and CFU directory from command line arguments
input <- args[1]
outdir <- args[2]

# Read data
cfu_data = read.csv(here(input, "cfu.csv")) %>% 
    mutate(cfu_log = log10(cfu),
           sample = factor(sample))

# Export data
write.csv(cfu_data, here(outdir, "data_processed.csv"), row.names = FALSE)

# Summary data
cfu_data %>% 
    group_by(synID, syncom, strain, dpi, taxa) %>% 
    summarise(mean = mean(cfu_log),
              sd = sd(cfu_log),
              cv = sd/mean*100,
              .groups = "drop") %>% 
    write.csv(here(outdir, "cfu_data_summary.csv"), row.names=FALSE)

# Total CFU per sample
cfu_data %>% 
    group_by(exp, synID, comID, taxa, dpi, sample) %>% 
    summarise(cfu_total = sum(cfu),
              log_total = log10(cfu_total),
              .groups = "drop") %>% 
    write.csv(here(outdir, "cfu_total.csv"), row.names=FALSE)