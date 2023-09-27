## Spatial fractions

#   Open file
data <- readRDS(here('results', 'stat_K_inhom_table.rds'))

data %>%
    mutate(
        aggregate = case_when(obs > hi ~ 1, TRUE ~ 0),
        regular = case_when(obs < lo ~ 1, TRUE ~ 0),
        random = case_when(obs < hi & obs > lo ~ 1, TRUE ~ 0)) %>% 
    group_by(syncom, strain, dpi, r) %>% 
    summarise(
        aggregate = sum(aggregate),
        regular = sum(regular),
        random = sum(random),
        .groups = 'drop') %>% 
    mutate(
        total = aggregate + regular + random,
        aggregate_fraction = if_else(is.na(aggregate/total), 0, aggregate/total),
        regular_fraction = if_else(is.na(regular/total), 0, regular/total),
        random_fraction = if_else(is.na(random/total), 0, random/total)) %>% 
    pivot_longer(cols = c("aggregate_fraction", "regular_fraction", "random_fraction"), 
                 names_to = "type", 
                 values_to = "fraction", 
                 values_drop_na = TRUE) %>% 
    select(syncom, strain, dpi, r, type, fraction) %>% 
    mutate(strain = factor(strain),
           dpi = factor(dpi),
           type = factor(type, levels = c('regular_fraction', 'random_fraction', 'aggregate_fraction'))) %>% 
    write.csv(here('results', 'stat_K_fractions.csv'))