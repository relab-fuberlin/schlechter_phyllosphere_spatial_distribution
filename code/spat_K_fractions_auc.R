#!/usr/bin/env Rscript 

##  Area under the curve spatial patterns frequency plots from K-estimates
#   Load required libraries
library(here)
library(tidyverse)
library(DescTools)

# Load data
fractions <- read.csv(here('results', 'stat_K_fractions.csv'), header = TRUE, row.names = "X") %>% 
    tibble %>% 
    separate(col='syncom', into = c('synID', 'synC'), sep = '\\.', remove = FALSE) %>% 
    select(-synC) %>% 
    mutate(taxa = case_when(
        strain == "meL85" ~ "Methylobacterium",
        strain == "meL92" ~ "Methylobacterium",
        strain == "mr01" ~ "Methylobacterium",
        TRUE ~ "Sphingomonas"))

auc_fractions <- fractions %>% 
    group_by(syncom, synID, taxa, dpi, type, strain) %>% 
    summarise(auc = AUC(r, fraction, absolutearea = FALSE),
              .groups = 'drop')

auc_C <- auc_fractions %>% 
    filter(synID == 'C') %>% 
    select(-syncom:-taxa)

auc_fold_change <- auc_fractions %>% 
    filter(synID != 'C') %>% 
    inner_join(., auc_C, by = c('dpi', 'type', 'strain'), suffix = c('.inter', '.intra')) %>% 
    mutate(fractional_change = (auc.inter-auc.intra)/auc.intra)

