# Calculates log2 fold change (Log2FC)

fold_change <- function(df, value, ref_variable, ref_group, group_by_variable = NULL, additional_variable = NULL){
    
    # This function calculates the fold change of the median of cell numbers
    # with different groups and explanatory variable
    # 1st arg: dataframe
    # 2nd arg: str character of name response variable
    # 3nd arg: vector of column names as grouping factor(s)
    # 4th arg: str character of name explanatory variable
    
    if(!is.null(group_by_variable) && length(group_by_variable) > 0) {
        # Calculate the median for the control group
        median_control = df %>% 
            filter(!!sym(ref_variable) == ref_group) %>% 
            group_by(across(all_of(group_by_variable))) %>% 
            summarise(m_mono = median(!!!syms(value), na.rm = TRUE), .groups="drop")
        
        # Calculate the median for all groups and join with control median
        result <- df %>% 
            group_by(across(all_of(group_by_variable)), 
                     !!!syms(ref_variable),
                     !!!syms(additional_variable)) %>% 
            summarise(m = median(!!!syms(value)), .groups="drop") %>% 
            inner_join(median_control, by = group_by_variable) %>% 
            mutate(FC = m/m_mono,
                   log2FC= ifelse(FC != 0, log2(FC), NA))
        
    } else {
        # Calculate the median for the control group
        median_control = df %>% 
            filter(!!sym(ref_variable) == ref_group) %>% 
            summarise(m_mono = median(!!!syms(value), na.rm = TRUE), .groups="drop")
        
        # Calculate the median for all groups and join with control median
        result <- df %>% 
            group_by(!!!syms(ref_variable),
                     !!!syms(additional_variable)) %>% 
            summarise(m = median(!!!syms(value)), .groups="drop") %>% 
            inner_join(median_control, by = character()) %>% 
            mutate(FC = m/m_mono,
                   log2FC= ifelse(FC != 0, log2(FC), NA))
        
    }
    
    return(result)
}
