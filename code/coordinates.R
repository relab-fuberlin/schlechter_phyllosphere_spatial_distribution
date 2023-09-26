#!/usr/bin/env Rscript

# Data processing cell coordinates
library(tidyverse)

metadata = read_csv(here('data', 'comm_id.csv'))
data.S2 = read.csv(here('data', 'coordinates_S2.csv'), header = TRUE, row.names = "X")
data.S3 = read.csv(here('data', 'coordinates_S3.csv'), header = TRUE, row.names = "X")
data.all = rbind(data.S2, data.S3)

data.all %>%
    mutate(synID = case_when(
        syn == "s2" & comm == "comm01" ~ "C",
        syn == "s2" & comm == "comm04" ~ "C",
        syn == "s2" & comm == "comm07" ~ "C",
        syn == "s2" & comm == "comm10" ~ "C",
        syn == "s2" & comm == "comm13" ~ "C",
        TRUE ~ syn),
        synID = str_replace(synID, "s", "S"),
        dpi = str_replace(dpi, "7d", "07dpi"),
        dpi = str_replace(dpi, "14d", "14dpi"),
        comID = str_replace(comm, "comm", "Com"),
        comm_id = str_replace(comID, "Com", ""),
        syncom = paste(synID, comm_id, sep=".")) %>%
    rename(channel = ch) %>% 
    left_join(., metadata, by = c("exp", "dpi", "synID", "comID", "syncom", "channel")) %>% 
    select(exp, dpi, synID, comID, syncom, strain, channel, img, obj_n, area, x, y) %>% 
    saveRDS(file = here('results', 'coordinates.rds'))
