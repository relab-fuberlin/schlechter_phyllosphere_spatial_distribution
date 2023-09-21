# Perform Dunnet test
dun_func1 <- function(dataframe, cell_log, group_by_variable, explanatory_variable){
    
    #.This function works for dataframe of the structure of data_cfu (first argument) for performing Dunnet test
    # with different groups and explanatory variable
    # 1st arg: dataframe
    # 2nd arg: str character of name response variable
    # 3nd arg: vector of column names as grouping factor(s)
    # 4th arg: str character of name explanatory variable
    
    col_index_cell <- which(colnames(dataframe)==cell_log)
    col_index <- which(colnames(dataframe)==explanatory_variable)
    
    formula_str <- paste(colnames(dataframe)[col_index_cell], "~", colnames(dataframe)[col_index])
    
    result <- dataframe %>% 
        group_by( !!!syms(group_by_variable)) %>% 
        dunn_test(as.formula(formula_str), p.adjust.method = 'holm', detailed = TRUE) %>% 
        select(!!!syms(group_by_variable), group1, group2, estimate, statistic, p.adj) %>% 
        group_by(!!!syms(group_by_variable)) %>% 
        mutate(p.adj_new = case_when(p.adj < 0.05 ~ 0.05, TRUE ~ p.adj))
    
    return(result)
}


# Calculates log2 fold change (Log2FC)
fold_func1 <- function(dataframe, cell, group_by_variable, explanatory_variable){
    
    # This function calculates the fold change of the median of cell numbers
    # with different groups and explanatory variable
    # 1st arg: dataframe
    # 2nd arg: str character of name response variable
    # 3nd arg: vector of column names as grouping factor(s)
    # 4th arg: str character of name explanatory variable
    
    median_control = data_cfu %>% 
        filter(synID == "C") %>% 
        group_by(!!!syms(group_by_variable)) %>% 
        summarise(m_mono = median(!!!syms(cell)), .groups="drop")
    
    result <- dataframe %>% 
        group_by(!!!syms(group_by_variable), !!!syms(explanatory_variable)) %>% 
        summarise(m = median(!!!syms(cell)), .groups="drop") %>% 
        inner_join(., median_control, by = c(group_by_variable)) %>% 
        mutate(fold = log2(m/m_mono))
    
    return(median_control)
}