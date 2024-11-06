# Perform Dunnet test
dunnet_function <- function(dataframe, cell_log, group_by_variable, explanatory_variable){
    
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
        group_by(!!!syms(group_by_variable)) %>% 
        dunn_test(as.formula(formula_str), p.adjust.method = 'holm', detailed = TRUE) %>% 
        dplyr::select(!!!syms(group_by_variable), group1, group2, estimate, statistic, p.adj) %>% 
        group_by(!!!syms(group_by_variable)) %>% 
        mutate(p_label = case_when(p.adj < 0.05 ~ "< 0.05", TRUE ~ as.character(p.adj)),
               p_size = case_when(p.adj < 0.05 ~ 0.05, TRUE ~ p.adj))
    
    return(result)
}
