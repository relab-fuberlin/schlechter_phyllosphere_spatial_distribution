

area = 124.94*100.24*1e-12 #m^2

data <- readRDS(here('results', 'coordinates.rds'))

metadata <- read.csv(here('data', 'syncom_metadata.csv')) %>% 
    pivot_longer(cols=c(C0, C1, C2), names_to = "channel", values_to = "strain") %>% 
    na.omit

combined_data <- data %>% 
    group_by(channel, exp, synID, comID, syncom, img) %>% 
    tally %>% 
    left_join(metadata, ., by = c("channel", "exp", "synID", "comID", "syncom", "img")) %>% 
    na.omit

nFOV <- combined_data %>% 
    pivot_wider(id_cols = c(exp, dpi, synID, comID, syncom, img), 
                names_from = channel, values_from = n, values_fill = 0, values_fn = mean) %>% 
    group_by(exp, dpi, synID, comID, syncom) %>% 
    summarise(nFOV = max(img),
              nC0 = sum(C0),
              nC1 = sum(C1),
              nC2 = sum(C2),
              .groups = "drop") %>% 
    mutate(sumCell = nC0 + nC1 + nC2,
           avCell = sumCell/nFOV)

coverage = combined_data %>% 
    mutate(cm2 = area) %>% 
    group_by(exp, dpi, synID, comID, syncom, strain, rep) %>% 
    summarise(cell = sum(n),
              total_area = sum(cm2),
              .groups = "drop") %>% 
    mutate(cell_area = cell/total_area,
           logCell = log10(cell_area))

cell_summary = coverage %>% 
    group_by(dpi, synID, comID, syncom, strain) %>% 
    summarise(mean = mean(logCell),
              sd = sd(logCell),
              cv = sd/mean*100,
              .groups = "drop")
