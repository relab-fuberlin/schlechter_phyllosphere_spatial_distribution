read.csv("results/cfu_data_processed.csv") %>% 
    na.omit %>% 
    ggplot(aes(synID, cfu_log, fill=synID))+
    
    stat_eye(
        side="right",
        adjust = 1,
        justification = -0.3,
        .width = 0,
        scale = 0.5,
        point_colour = NA
    )+
    
    geom_jitter(width = 0.1,
                alpha = 0.1)+
    
    geom_boxplot(
        fill="white",
        width=0.2,
        outlier.alpha = 0
    )+
    
    coord_cartesian(xlim=c(1,3.2))+
    
    stat_compare_means(
        aes(label = paste0(
            "Anova ~italic(p)",
            scales::label_pvalue(accuracy = 0.05)(..p..)
        )),
        parse = TRUE,
        method = "anova", 
        label.y = 11.8, 
        label.x = 2.3,
        size = 5)+
    
    stat_compare_means(
        aes(label = after_stat(p.signif)), 
        method = "t.test", 
        comparisons = list(c("C", "S2"), c("S2", "S3"), c("C", "S3")),
        size = 8,
        bracket.size = 0.5,
        symnum.args = list(cutpoints = c(0,0.05)))+
    
    scale_y_continuous(limits = c(2,12), breaks=c(2,4,6,8,10))+
    scale_fill_manual(name = "SynCom", values = syn.pal, labels = syn.lab, guide="none")+
    labs(y = bquote('Population density ['*log[10]~ "CFU g" ~ FW^-1 *"]"),
         x = "SynCom")+
    theme_rs()


read.csv("results/cfu_data_processed.csv") %>% 
    na.omit %>% 
    ggplot(aes(dpi, cfu_log, fill=synID))+
    facet_wrap(~taxa, ncol=2)+
    
    stat_halfeye(
        adjust = 1,
        justification = -0.3,
        .width = 0,
        scale = 0.6,
        point_colour = NA,
        position = position_dodge(0.9)
    )+

    geom_point(shape=21, 
        position=position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1),
        alpha = 0.1,
        size = 1)+
    
    geom_boxplot(
        aes(group=interaction(dpi,synID)),
        fill="white",
        width=0.2,
        position=position_dodge(0.9),
        outlier.alpha = 0
    )+
    scale_y_continuous(name = "Population density\n[log10 CFU g FW-1]", limits = c(2,11), expand=c(0,0), breaks=c(2,4,6,8,10))+
    scale_x_discrete(name = "Days post-inoculation [dpi]", labels = dpi.lab)+
    theme_rs()+
    scale_fill_manual(name = "SynCom", values = syn.pal, labels = syn.lab, guide="none")


read.csv("results/cfu_data_processed.csv") %>% 
    na.omit %>% 
    ggplot(aes(dpi, cfu_log))+
    
    stat_eye(
        side="right",
        adjust = 1,
        justification = -0.3,
        .width = 0,
        scale = 0.5,
        point_colour = NA
    )+
    
    geom_jitter(width = 0.1,
                alpha = 0.5)+
    
    geom_boxplot(
        fill="white",
        width=0.08,
        outlier.alpha = 0
    )+
    
    coord_cartesian(xlim=c(1.2,2.2))+
    
    stat_compare_means(aes(label = after_stat(p.signif)), 
                       method = "t.test", 
                       comparisons = list(c("07dpi", "14dpi")),
                       size = 8,
                       bracket.size = 0.5,
                       symnum.args = list(cutpoints = c(0,0.05)))+
    
    scale_y_continuous(name = "Population density\n[log10 CFU g FW-1]", limits = c(2,11), breaks=c(2,4,6,8,10))+
    
    theme_rs()


