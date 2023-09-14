#!/usr/bin/env Rscript

# Data processing CFU data
source("code/libraries_syncom.R")

cfu_data = read.csv(here("data/cfu.csv")) %>% 
    mutate(cfu_log = log10(cfu),
           sample = factor(sample))

write.csv(cfu_data, here("results/cfu_data_processed.csv"), row.names=FALSE)

# Summary data
cfu_data %>% 
    group_by(synID, syncom, strain, dpi, taxa) %>% 
    summarise(mean = mean(cfu_log),
              sd = sd(cfu_log),
              cv = sd/mean*100,
              .groups = "drop") %>% 
    write.csv(here("results/cfu_data_summary.csv"), row.names=FALSE)

# Total CFU per sample
cfu_data %>% 
    group_by(exp, synID, comID, taxa, dpi, sample) %>% 
    summarise(cfu_total = sum(cfu),
              log_total = log10(cfu_total),
              .groups = "drop") %>% 
    write.csv(here("results/cfu_total.csv"), row.names=FALSE)