# Define a custom function to calculate the ECDF for a dataset
ecdf_func <- function(data) {
    # Use the ecdf() function to compute the ECDF
    ecdf_object <- ecdf(data)
    value <- ecdf_object(data)
    return(value)
}

#   Open file
data <- readRDS(here('results', 'stat_K_inhom_table.rds')) %>% 
    mutate(taxa = case_when(
        strain == "meL85" ~ "Methylobacterium",
        strain == "meL92" ~ "Methylobacterium",
        strain == "mr01" ~ "Methylobacterium",
        TRUE ~ "Sphingomonas"))

## Empirical cumulative distribution
ecdf_all <- data %>% 
    filter(obs > hi) %>% 
    group_by(dpi, syncom, taxa, strain, exp, rep) %>% 
    mutate(ecdfFun = ecdf_func(r)) %>% 
    filter(ecdfFun > 0.95) %>% 
    summarise(r_min = min(r), .groups="drop") %>% 
    group_by(dpi, syncom, taxa, strain) %>% 
    summarise(r = mean(r_min), .groups="drop") %>% 
    separate(col="syncom", into=c("synID", "com"), sep="\\.", remove = FALSE)
    
r_c <- ecdf_all %>% 
    filter(synID == "C") %>% 
    select(-syncom:-com)

r_fold <- ecdf_all %>% 
    filter(synID != "C") %>% 
    left_join(., r_c, by=c("dpi", "strain", "taxa"), suffix = c(".inter", ".intra")) %>% 
    mutate(fold = r.inter/r.intra,
           log2FC = log2(fold))

r_fold_s2 <- r_fold %>% 
    filter(synID=="S2") %>% 
    select(syncom, dpi, strain, r.inter, log2FC) %>% 
    group_by(syncom, dpi) %>% 
    mutate(to=rev(strain)) %>% 
    arrange(to) %>% 
    pivot_wider(id_cols = c(strain, dpi), names_from=to, values_from=c(r.inter,log2FC)) %>% 
    arrange(strain)

r_fold_s3 <- r_fold %>% 
    filter(synID=="S3") %>% 
    pivot_wider(id_cols = c(strain, dpi), names_from=syncom, values_from=c(r.inter,log2FC)) %>% 
    arrange(strain)

r_fold_back <- cbind(r_fold_s2, r_fold_s3[,c(-1:-2)]) %>% 
    pivot_longer(cols=-c(strain,dpi), names_to = c("type", "strain2"), names_sep="_", ) %>% 
    pivot_wider(id_cols = c(strain, dpi, strain2), names_from = type, values_from = value) %>% 
    ungroup %>% 
    arrange(dpi,strain) %>% 
    mutate(index = seq(1, nrow(.)))