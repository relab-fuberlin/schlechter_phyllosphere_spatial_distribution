##  Pair correlation analysis

#   Set seed
set.seed(19900725)

#   Dependencies
library(here)
library(spatstat)
library(tidyverse)
source(here('code', 'function_pcf.R')) # Custom PCF function
source(here('code', 'palette_syncom.R'))

#   Data
H <- readRDS(here('results', 'hyperframe_syncom.rds')) # Read hyperframe

#   PCF for S2
H_S2 <- H[grep("S2", H$syncom)] # Filter by S2 communities. C and S3 are excluded
H_S2 <- H_S2[H_S2$C0 > 10 & H_S2$C1 > 10] # Filter by FOV with at least 10 cells per channel
pcf_S2_C0C1 <- pcf_function(H_S2, i = "C0", j = "C1")
pcf_S2_C0C1 <- pcf_S2_C0C1 %>% mutate(pair = "C0.C1", pair_inv = "C1.C0")

#   PCF for S3
H_S3 <- H[grep("S3", H$syncom)] # Filter by S3 
pcf_S3_C0C1 <- pcf_function(H_S3[H_S3$C0 > 10 & H_S3$C1 >10,], i = "C0", j = "C1")
pcf_S3_C0C1 <- pcf_S3_C0C1 %>% mutate(pair = "C0.C1", pair_inv = "C1.C0")

pcf_S3_C0C2 <- pcf_function(H_S3[H_S3$C0 > 10 & H_S3$C2 >10,], i = "C0", j = "C2")
pcf_S3_C0C2 <- pcf_S3_C0C2 %>% mutate(pair = "C0.C2", pair_inv = "C2.C0")

pcf_S3_C1C2 <- pcf_function(H_S3[H_S3$C1 > 10 & H_S3$C2 >10,], i = "C1", j = "C2")
pcf_S3_C1C2 <- pcf_S3_C1C2 %>% mutate(pair = "C1.C2", pair_inv = "C2.C1")


#   Combine results of PCF for S2 and S3 communities
pcf <- rbind(pcf_S2_C0C1, pcf_S3_C0C1, pcf_S3_C0C2, pcf_S3_C1C2) %>% 
    pivot_longer(pair:pair_inv, names_to="type_pair", values_to="pair") %>% 
    separate(rep, into=c("syncom", "dpi", "exp", "img"), sep="_")

#   Add metadata
metadata <- read.csv(here('data', 'comm_id.csv')) %>% 
    group_by(exp, dpi, syncom, strain, channel) %>%
    unique() %>% 
    pivot_wider(id_cols = exp:syncom, names_from=channel, values_from = strain) %>% 
    unite("C0.C1", C0, C1, sep = ".", remove = FALSE) %>% 
    unite("C1.C0", C1, C0, sep = ".", remove = FALSE) %>% 
    unite("C0.C2", C0, C2, sep = ".", remove = FALSE) %>% 
    unite("C2.C0", C2, C0, sep = ".", remove = FALSE) %>% 
    unite("C1.C2", C1, C2, sep = ".", remove = FALSE) %>% 
    unite("C2.C1", C2, C1, sep = ".", remove = TRUE) %>% 
    select(-C0) %>% 
    pivot_longer(cols="C0.C1":"C2.C1", names_to = 'pair', values_to = "strain_pair") %>% 
    filter(!grepl("NA", strain_pair))

stat_pcf <- left_join(pcf, metadata, by=c("exp", "dpi", "syncom", "pair")) %>% 
    select(syncom, synID, comID, dpi, exp, img, r, obs, theo, lo, hi, strain_pair) %>% 
    filter(strain_pair %in% unique_pair)

#   Save data as RDS
write_rds(stat_pcf, here('results', 'stat_pcf_inhom_table.rds'))

