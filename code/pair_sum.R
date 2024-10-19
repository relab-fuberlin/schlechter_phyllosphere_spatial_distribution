source(here("code", "function_triad.R"))


# Step 1: Split the data
split_data <- split(data_cfu, 
                    list(data_cfu$synID, data_cfu$syncom, data_cfu$exp, data_cfu$sample, data_cfu$dpi), 
                    drop = TRUE,
                    sep = "_")

# Step 2: Apply the function to each group
result_list <- lapply(split_data, process_triad)

# Step 3: Combine the results into a single dataframe
final_result <- do.call(rbind, result_list)

# Print the final result
cfu_pair <- final_result %>% 
    rownames_to_column(var = "label") %>% 
    separate(label, into = c("synID", "syncom", "exp", "sample", "dpi"), sep = "_") %>% 
    separate(dpi, into = c("dpi", "combn"), convert = TRUE) %>% 
    unite(strain1:strain2, col = "pair") %>% 
    group_by(synID, syncom, dpi, pair, excluded) %>% 
    summarise(cfu_median = median(sum), .groups = "drop") %>%
    filter(synID != "C") %>% 
    arrange(pair)


pair_S2 <- cfu_pair %>% 
    filter(synID == "S2") %>% 
    select(dpi, pair, cfu_median)

pair_S3 <- cfu_pair %>% 
    filter(synID == "S3" & excluded != 'NA') %>% 
    left_join(., pair_S2, by = c("dpi", "pair"), suffix = c(".S3", ".S2")) %>% 
    mutate(fc_pair = log2(cfu_median.S3/cfu_median.S2))

pair_S3 %>% 
    ggplot(., aes(x = pair, y = excluded))+
    facet_grid(rows = vars(dpi), scales = "free_x")+
    geom_tile(aes(fill = fc_pair))+
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
    scale_fill_gradientn(name = bquote(Log[2]~"FC"), 
                         colours = wes_palette("Zissou1")[c(1,2,3,4,5)], 
                         values = c(0,0.57,1), limits=c(-5,0), 
                         breaks = seq(-5,0,1), na.value = 'grey90')

pair_S3 %>% 
    mutate(
        taxa_pair = case_when(
            pair == "meL85_meL92" | pair == "meL85_meL92" | pair == "meL85_mr01" | pair == "meL92_mr01" ~ "MM",
            pair == "meL85_smfr1" | pair == "meL85_spfa2" | pair == "meL92_smfr1" | pair == "meL92_spfa2" | pair == "mr01_smfr1" | pair == "mr01_spfa2" ~ "MS",
            pair == "smfr1_spfa2" ~ "SS"),
        taxa_excluded = case_when(
            excluded == "smfr1" | excluded == "spfa2" ~ "S",
            excluded == "meL85" | excluded == "meL92" | excluded == "mr01" ~ "M"
        )) %>%
    ggplot(., aes(taxa_pair, fc_pair, fill = taxa_excluded))+
    facet_grid(cols = vars(dpi))+
    geom_boxplot(position = position_dodge(width = 0.75, preserve = "single"))+
    geom_point(aes(group = taxa_excluded), position = position_dodge(width = 0.75))

pair_S3 %>% 
    mutate(
        taxa_pair = case_when(
            pair == "meL85_meL92" | pair == "meL85_meL92" | pair == "meL85_mr01" | pair == "meL92_mr01" ~ "MM",
            pair == "meL85_smfr1" | pair == "meL85_spfa2" | pair == "meL92_smfr1" | pair == "meL92_spfa2" | pair == "mr01_smfr1" | pair == "mr01_spfa2" ~ "MS",
            pair == "smfr1_spfa2" ~ "SS"),
        taxa_excluded = case_when(
            excluded == "smfr1" | excluded == "spfa2" ~ "S",
            excluded == "meL85" | excluded == "meL92" | excluded == "mr01" ~ "M"
        )) %>%
    group_by(synID, dpi, taxa_pair, excluded) %>% 
    summarise(fc_pair = mean(fc_pair), .groups = "drop") %>% 
    ggplot(., aes(taxa_pair, excluded))+
    facet_grid(rows = vars(dpi))+
    geom_tile(aes(fill = fc_pair))+
    geom_text(aes(label = sprintf(fc_pair, fmt = '%#.2f')))+
    scale_fill_gradientn(name = bquote(Log[2]~"FC"), 
                         colours = wes_palette("Zissou1")[c(1,2,3,4,5)], 
                         values = c(0,0.57,1), limits=c(-5.1,0), 
                         breaks = seq(-5,0,1), na.value = 'grey90')




subset_df <- split_data[["S3_S3.10_e2_4_14dpi"]]

process_triad(subset_df)

strain <- subset_df$strain
cfu <- subset_df$cfu
    
# Generate all combinations of 3 objects
comb3 <- combn(strain, 3, simplify = FALSE)
        
# Create a dataframe to store results
results <- data.frame(strain1 = character(),
                      strain2 = character(),
                      excluded = character(),
                      sum = numeric())
        
# Loop through each combination of 3 objects
for (combo in comb3) {
        # Generate all pairs from the combination
    pair_combinations <- as.data.frame(t(combn(combo, 2)))
    colnames(pair_combinations) <- c("strain1", "strain2")
        
        # Calculate the sum of the values for each pair
    pair_combinations$sum <- apply(pair_combinations, 1, function(row) {
        sum(cfu[strain %in% row])
        })
            
        # Add the excluded object (the one not included in the pair)
    pair_combinations$excluded <- sapply(1:nrow(pair_combinations), function(i) {
        setdiff(combo, pair_combinations[i, c("strain1", "strain2")])
        })
            
        # Combine results
    results <- rbind(results, pair_combinations) 
}
        
    } else if (length(strain) == 2) {
        # If only two objects, handle as a pair
        pair_combinations <- data.frame(strain1 = strain[1], strain2 = strain[2], 
                                        excluded = NA,
                                        sum = sum(cfu))
        results <- pair_combinations
    }
    return(results)
}
