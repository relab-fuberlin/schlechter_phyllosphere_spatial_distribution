fractional_change_median = function(df, var="fractional_change", group=NULL){
    summary <- df %>% 
        group_by(!!!syms(group)) %>% 
        summarise(median = median(!!!syms(var)),
                  q1 = format(round(quantile(!!!syms(var), 0.25), 2), nsmall = 2),
                  q3 = format(round(quantile(!!!syms(var), 0.75), 2), nsmall = 2),
                  IQR = paste0(q1,"-(",q3,")", sep=''),
                  percentage = abs(100*median),
                  .groups = "drop")
    
    return(summary)
}