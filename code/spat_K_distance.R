#   Create the function.
getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}

#   Open file
data <- readRDS(here('results', 'stat_K_inhom_table.rds'))

## total empirical cumulative distribution
agg_cdf = data[data$obs > data$hi,]
ecdf_k = ecdf(agg_cdf$r)
plot(ecdf_k)
agg_cdf$ecdf = ecdf_k(agg_cdf$r)
