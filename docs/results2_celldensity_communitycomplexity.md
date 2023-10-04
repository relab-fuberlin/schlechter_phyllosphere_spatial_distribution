---
title: "Spatial distribution paper - Section 2"
author: "Rudolf Schlechter"
output:
  html_document:
    df_print: paged
    keep_md: yes
  pdf_document: default
editor_options: 
  chunk_output_type: console
---



## Spatial distribution of individual strains depend on their community context




```r
data_cell %>% head
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["exp"],"name":[1],"type":["chr"],"align":["left"]},{"label":["dpi"],"name":[2],"type":["fct"],"align":["left"]},{"label":["synID"],"name":[3],"type":["fct"],"align":["left"]},{"label":["comID"],"name":[4],"type":["chr"],"align":["left"]},{"label":["syncom"],"name":[5],"type":["chr"],"align":["left"]},{"label":["strain"],"name":[6],"type":["fct"],"align":["left"]},{"label":["sample"],"name":[7],"type":["int"],"align":["right"]},{"label":["cell"],"name":[8],"type":["int"],"align":["right"]},{"label":["total_area"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["cell_area"],"name":[10],"type":["dbl"],"align":["right"]},{"label":["logCell"],"name":[11],"type":["dbl"],"align":["right"]},{"label":["taxa"],"name":[12],"type":["fct"],"align":["left"]},{"label":["channel"],"name":[13],"type":["chr"],"align":["left"]},{"label":["cfu"],"name":[14],"type":["dbl"],"align":["right"]},{"label":["cfu_log"],"name":[15],"type":["dbl"],"align":["right"]}],"data":[{"1":"e1","2":"07dpi","3":"C","4":"Com01","5":"C.01","6":"meL85","7":"1","8":"34761","9":"0.0063","10":"5551108","11":"6.7","12":"Methylobacterium","13":"C0","14":"8800000","15":"6.9","_rn_":"1"},{"1":"e1","2":"07dpi","3":"C","4":"Com01","5":"C.01","6":"meL85","7":"1","8":"34761","9":"0.0063","10":"5551108","11":"6.7","12":"Methylobacterium","13":"C1","14":"8350000","15":"6.9","_rn_":"2"},{"1":"e1","2":"07dpi","3":"C","4":"Com01","5":"C.01","6":"meL85","7":"2","8":"34587","9":"0.0103","10":"3367879","11":"6.5","12":"Methylobacterium","13":"C0","14":"15400000","15":"7.2","_rn_":"3"},{"1":"e1","2":"07dpi","3":"C","4":"Com01","5":"C.01","6":"meL85","7":"2","8":"34587","9":"0.0103","10":"3367879","11":"6.5","12":"Methylobacterium","13":"C1","14":"5760000","15":"6.8","_rn_":"4"},{"1":"e1","2":"07dpi","3":"C","4":"Com01","5":"C.01","6":"meL85","7":"3","8":"30764","9":"0.0060","10":"5117514","11":"6.7","12":"Methylobacterium","13":"C0","14":"19900000","15":"7.3","_rn_":"5"},{"1":"e1","2":"07dpi","3":"C","4":"Com01","5":"C.01","6":"meL85","7":"3","8":"30764","9":"0.0060","10":"5117514","11":"6.7","12":"Methylobacterium","13":"C1","14":"23900000","15":"7.4","_rn_":"6"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>


```r
# Correlation between bacterial densities determined by CFU or by cell counts
data_cell %>% 
    cor_test(
        vars = c("logCell", "cfu_log"),
        method = "pearson",
        conf.level = 0.95)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["var1"],"name":[1],"type":["chr"],"align":["left"]},{"label":["var2"],"name":[2],"type":["chr"],"align":["left"]},{"label":["cor"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["statistic"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["p"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["conf.low"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["conf.high"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["method"],"name":[8],"type":["chr"],"align":["left"]}],"data":[{"1":"logCell","2":"cfu_log","3":"0.37","4":"12","5":"9.2e-29","6":"0.31","7":"0.43","8":"Pearson"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>



```r
linear_cell <-  lm(logCell ~ synID + dpi + taxa, data_cell)

# Shapiro-Wilk test for normality
shapiro.test(rstandard(linear_cell))
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  rstandard(linear_cell)
## W = 0.9, p-value <2e-16
```

```r
# Breusch-Pagan test for homogeneity of variances
ncvTest(linear_cell)
```

```
## Non-constant Variance Score Test 
## Variance formula: ~ fitted.values 
## Chisquare = 102, Df = 1, p = <2e-16
```


```r
# Kruskal-Wallis test and effect size for community complexity (synID)
kw_synID <- data_cell %>%
    group_by(taxa) %>% 
    kruskal_test(logCell ~ synID) %>% 
    mutate(p_label = case_when(p < 0.05 ~ "< 0.05", TRUE ~ as.character(p)))
keff_synID <- data_cell %>% 
    group_by(taxa) %>% 
    kruskal_effsize(logCell ~ synID, ci=TRUE, nboot=100)

# Dunn's Test
dunn_synID <- data_cell %>% 
    group_by(taxa) %>% 
    dunn_test(logCell ~ synID, p.adjust.method = "holm") %>% 
    mutate(p_label = case_when(p.adj < 0.05 ~ "< 0.05", TRUE ~ as.character(p.adj)))
```


```r
# Wilcoxon test and effect size for sampling time (dpi)
w_dpi = data_cell %>% 
    wilcox_test(formula = logCell ~ dpi)

# Fold change of population density by time of sampling (dpi)
fc_cell_dpi = data_cell %>% 
    group_by(dpi) %>% 
    summarise(median_cell = median(cell_density)) %>% 
    mutate(FC = median_cell/median_cell[1],
           logFC = log2(FC))
```







```r
areas <- c(patchwork::area(1,1,3), patchwork::area(1,2,1), patchwork::area(2,2,3))
wrap_elements(full = plt3.a) + plt3.b + plt3.c  + 
    plot_annotation(tag_levels = "A") + 
    plot_layout(guides = "collect", design = areas) & theme(legend.box.just = "center")
```

<div class="figure" style="text-align: center">
<img src="results2_celldensity_communitycomplexity_files/figure-html/figure_main_5-1.png" alt="Bacterial cell density in the arabidopsis phyllosphere" id="fig-figure-main-5-1" />
<p class="caption">Bacterial cell density in the arabidopsis phyllosphere</p>
</div>
