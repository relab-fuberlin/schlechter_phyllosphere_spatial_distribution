
process_triad <- function(subset_df, cell) {
    strain <- subset_df$strain
    cfu <- subset_df[[cell]]
    
    if (length(strain) == 3) {
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
        
    } else if (length(strain) == 1) {
        # If only two objects, handle as a pair
        pair_combinations <- data.frame(strain1 = strain[1], strain2 = NA, 
                                        excluded = NA,
                                        sum = sum(cfu))
        results <- pair_combinations
    }
    return(results)
}
