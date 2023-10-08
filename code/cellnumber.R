#!/usr/bin/env Rscript

library(tidyverse)

# Area in cm2
area = 124.94*100.24*1e-8

# Coordinates data
# Each entry represents one individual cell
data <- readRDS(here('results', 'coordinates.rds'))

# Metadata includes biological replicates (rep) and number of images taken
metadata <- read.csv(here('data', 'metadata.csv')) %>% 
    pivot_longer(cols=c(C0, C1, C2), names_to = "channel", values_to = "strain") %>% 
    na.omit

# Combine data with metadata
combined_data <- data %>% 
    group_by(channel, exp, synID, comID, syncom, img) %>% 
    tally %>% 
    left_join(metadata, ., by = c("channel", "exp", "synID", "comID", "syncom", "img")) %>% 
    na.omit

# Number of cells per field of views
combined_data %>% 
    pivot_wider(id_cols = c(exp, dpi, synID, comID, syncom, img, rep), 
                names_from = channel, values_from = n, values_fill = 0) %>% 
    group_by(exp, dpi, synID, comID, syncom, rep) %>% 
    summarise(nFOV = max(img),
              sumCell = sum(C0) + sum(C1) + sum(C2),
              avCell = sumCell/nFOV,
              .groups = "drop")  %>% 
    write.csv(here('results', 'nfov.csv'), row.names = FALSE)

# Cell coverage (number of cells per area)
coverage <- combined_data %>% 
    group_by(exp, dpi, synID, comID, syncom, strain, rep, channel) %>%
    summarise(nFOV = max(img),
              cell = sum(n),
              cm2 = area * nFOV,
              cell_density = cell/cm2,
              logCell = log10(cell_density),
              .groups = "drop")

write.csv(coverage, here('results', 'cell_density.csv'), row.names = FALSE)

# Summary statistic of cell coverage
coverage %>% 
    group_by(dpi, synID, comID, syncom, strain, channel) %>% 
    summarise(mean = mean(logCell),
              sd = sd(logCell),
              cv = sd/mean*100,
              .groups = "drop") %>% 
    write.csv(here('results', 'cell_density_summary.csv'), row.names = FALSE)
