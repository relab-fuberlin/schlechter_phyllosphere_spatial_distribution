# Load necessary library
library(stats)

# Sample data: 3 dependent groups
group1 <- c(15, 20, 25, 30, 35)
group2 <- c(10, 25, 30, 28, 32)
group3 <- c(20, 24, 28, 29, 35)

# Combine data into a matrix
data_matrix <- cbind(group1, group2, group3)

# Friedman test
friedman.test(data_matrix)

list_full <- cfu_pair_full %>% 
    filter(synID == "S3") %>% 
    select(syncom, exp, sample, dpi, pair, sum) %>% 
    #unite(c(syncom, exp), col = "label", sep="_") %>% 
    split(., list(.$pair, .$dpi), drop = TRUE, sep = "-")

friedmantest <- lapply(list_full, pivot_wider, names_from = "syncom", values_from = "sum") %>% 
    lapply(., unite, c(exp, sample), col = "label", sep = "_") %>% 
    lapply(., column_to_rownames, var = "label") %>% 
    lapply(., select, -pair, -dpi) %>% 
    lapply(., as.matrix) %>% 
    lapply(., friedman.test) %>% 
    do.call(rbind, .) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "label") %>% 
    separate(label, into = c("pair", "dpi"), sep = "-")

list_full <- cfu_pair_full %>% 
    filter(synID == "S3") %>% 
    select(exp, sample, dpi, syncom, pair, sum) %>%
    split(., .$dpi, drop = TRUE)

friedman.test(sum ~ pair | syncom, data = list_full[[1]])

friedmantest <- lapply(list_full, pivot_wider, names_from = "pair", values_from = "sum") %>% 
    lapply(., select, -exp:-dpi) %>% 
    lapply(., as.matrix) %>% 
    lapply(., friedman.test, ) %>% 
    do.call(rbind, .) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "label") %>% 
    separate(label, into = c("pair", "dpi"), sep = "-")
    

lmm <- lmer(log10(sum) ~ pair * dpi * synID + (1 | syncom), data = cfu_pair_full %>% filter(synID!="C"), na.action = na.omit)
summary(lmm)

# Post-hoc comparisons
emm <- emmeans(lmm, ~ pair | dpi | synID)
emm_contrasts<- pairs(emm)

# Display the CLD result
print(cld_result)

model_pair <- cfu_pair_full %>% 
    filter(synID != "C") %>% 
    split(., list(.$dpi, .$synID), sep = "_") %>% 
    lapply(., lmer, formula = log10(sum) ~ pair + (1 | syncom))

emm_list_pair <- lapply(model_pair, emmeans, ~pair)
contrast_list_pair <- lapply(emm_list_pair, pairs, adjust = "BH") %>% 
    lapply(., data.frame)

do.call(rbind, contrast_list_pair) %>% 
    rownames_to_column(var = "label") %>% 
    separate(label, into = c("dpi", "syncom"), sep = "_") %>% 
    filter(p.value < 0.05)


##
models <- cfu_pair_full %>% 
    filter(synID != "C") %>% 
    split(., list(.$dpi, .$synID), sep = "_") %>% 
    lapply(., lmer, formula = log10(sum) ~ taxa_pair + (1 | syncom))

emm_list <- lapply(models, emmeans, specs = "taxa_pair")
contrast_list <- lapply(emm_list, cld, alpha = 0.05, Letters = letters, adjust = "BH") %>% 
    lapply(., data.frame) %>% 
    do.call(rbind, .) %>% 
    rownames_to_column(var = "label") %>% 
    separate(label, into = c("dpi", "syncom"), sep = "_") %>% 
    mutate(across(where(is.character), ~ str_remove_all(., "\\s+")),
           synID = case_when(
               str_detect(syncom, "S2") ~ "S2",
               TRUE ~ "S3"))

do.call(rbind, contrast_list) %>% 
    rownames_to_column(var = "label") %>% 
    separate(label, into = c("dpi", "syncom"), sep = "_") %>% 
    filter(p.value < 0.05)
    


percentage <- function(){
    
}

datatest <- split(data_cfu, list(data_cfu$dpi, data_cfu$syncom, data_cfu$exp, data_cfu$sample))


relative_fractions_list <- lapply(datatest, function(df) {
    df %>%
        mutate(Relative_Fraction = cfu / sum(cfu))
    })
