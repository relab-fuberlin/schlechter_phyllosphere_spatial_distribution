# Define a custom function to calculate the ECDF for a dataset
ecdf_func <- function(data) {
    # Use the ecdf() function to compute the ECDF
    ecdf_object <- ecdf(data)
    value <- ecdf_object(data)
    return(value)
}

#   Open file
data <- readRDS(here('results', 'stat_K_inhom_table.rds'))
data %>% head

## Empirical cumulative distribution
ecdf_all <- data %>% 
    filter(obs > hi) %>% 
    group_by(dpi, syncom, strain, exp, rep) %>% 
    mutate(ecdfFun = ecdf_func(r)) %>% 
    filter(ecdfFun > 0.95) %>% 
    summarise(r_min = min(r), .groups="drop") %>% 
    group_by(dpi, syncom, strain) %>% 
    summarise(r = mean(r_min), .groups="drop") %>% 
    separate(col="syncom", into=c("synID", "com"), sep="\\.", remove = FALSE)
    
r_c <- ecdf_all %>% 
    filter(synID == "C") %>% 
    select(-syncom:-com)
r_fold <- ecdf_all %>% 
    filter(synID != "C") %>% 
    left_join(., r_c, by=c("dpi", "strain"), suffix = c(".inter", ".intra")) %>% 
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
    

r_fold_back %>% 
    ggplot(aes(x=fct_reorder(strain2, index), y=strain))+
    facet_wrap(~dpi, ncol=1)+
    geom_tile(colour='black', fill='white')+
    geom_point(aes(fill=log2FC, size = r.inter), shape = 21)+
    scale_fill_gradientn(colours = wes_palette("Zissou1")[c(1,1,3,5,5)], values = c(0,0.5,1), limits=c(-0.4,0.4))+
    coord_fixed()+
    scale_size_continuous(range = c(2,12), limits=c(8,15))

ecdf_all %>% 
    ggplot(aes(x=strain, y=r))+
    geom_jitter(width=0.1)+
    geom_boxplot(alpha = 0.5, width=0.2)+
    scale_y_continuous(limits = c(6,16), breaks = seq(8,16,4), expand = c(0,0))+
    theme_classic()

ecdf_all %>% 
    kruskal_test(r ~ strain)

ecdf_all %>% 
    group_by(dpi) %>% 
    dunn_test(r ~ synID, p.adjust.method='holm')

