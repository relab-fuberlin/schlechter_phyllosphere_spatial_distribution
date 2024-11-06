#!/usr/bin/env Rscript 

##  Area under the curve spatial patterns frequency plots from PCF
#   Load required libraries
library(here)
library(tidyverse)
library(DescTools)

# Load data
G_fractions <- read.csv(here('results', 'stat_pcf_fractions.csv'), header = TRUE, row.names = "X") %>% 
    tibble %>% 
    separate(col='syncom', into = c('synID', 'synC'), sep = '\\.', remove = FALSE) %>% 
    dplyr::select(-synC) %>% 
    mutate(taxa_pair = case_when(
        strain_pair == "meL85.meL92" | strain_pair == "meL85.mr01" | strain_pair == "meL92.mr01" ~ "MM",
        strain_pair == "smfr1.spfa2" ~ "SS",
        TRUE ~ "MS"))

G_auc_fractions <- G_fractions %>% 
    group_by(syncom, synID, dpi, type, taxa_pair, strain_pair) %>% 
    summarise(auc = AUC(r, fraction, absolutearea = FALSE),
              .groups = 'drop')

G_auc_S2 <- G_auc_fractions %>% 
    filter(synID == 'S2') %>% 
    dplyr::select(dpi, type, strain_pair, auc)

G_auc_fold_change <- G_auc_fractions %>% 
    filter(synID != 'S2') %>% 
    inner_join(., G_auc_S2, by = c('dpi', 'type', 'strain_pair'), suffix = c('.S3', '.S2')) %>% 
    mutate(fractional_change = (auc.S3-auc.S2)/auc.S2) %>% 
    dplyr::select(-auc.S3:-auc.S2)
