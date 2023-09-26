read.csv(here('results', 'cfu_data_processed.csv')) %>% 
    na.omit %>% 
    ggplot(aes(synID, cfu_log, fill=synID))+
    stat_eye(
        side="right",
        adjust = 1,
        justification = -0.3,
        .width = 0,
        scale = 0.5,
        point_colour = NA)+
    geom_jitter(width = 0.1, alpha = 0.1)+
    geom_boxplot(fill="white", width=0.2, outlier.alpha = 0)+
    coord_cartesian(xlim=c(1,3.2))+
    stat_compare_means(
        aes(label = paste0(
            "Wilcoxon ~italic(p)",
            scales::label_pvalue(accuracy = 0.05)(..p..)
        )),
        parse = TRUE,
        method = "wilcox.test", 
        label.y = 10.8, 
        label.x = 2.8,
        size = 5)+
    
    stat_compare_means(
        aes(label = after_stat(p.signif)), 
        method = "t.test",
        p.adjust.method = "bonferroni",
        comparisons = list(c("C", "S2"), c("S2", "S3"), c("C", "S3")),
        vjust = 0.5,
        size = 8,
        bracket.size = 0.5,
        label.y = c(9,9.5,10),
        symnum.args = list(cutpoints = c(0,0.05)))+
    
    scale_y_continuous(limits = c(3,11), breaks=c(4,6,8,10))+
    scale_fill_manual(name = "SynCom", values = syn.pal, labels = syn.lab, guide="none")+
    labs(y = plt_bac_density_lab,
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
    scale_y_continuous(name = plt_bac_density_lab, limits = c(3,11), expand=c(0,0), breaks=c(4,6,8,10))+
    scale_x_discrete(name = plt_days_lab, labels = dpi.lab)+
    theme_rs()+
    scale_fill_manual(name = "SynCom", values = syn.pal, labels = syn.lab, guide="none")+
    
    stat_compare_means(
        aes(label = paste0(
            "Anova ~italic(p)",
            scales::label_pvalue(accuracy = 0.05)(..p..)
        )),
        parse = TRUE,
        method = "anova", 
        label.y = 8.8, 
        label.x = 2.8,
        size = 5) +
    
    geom_signif(
        inherit.aes = FALSE,
        data=anno_df,
        aes(group=interaction(synID,dpi), xmin=group1, xmax=group2, annotations=p.adj, y_position = 10),
        manual=TRUE
    )

anno_df = read.csv("results/cfu_data_processed.csv") %>% 
    na.omit %>%
    compare_means(data = ., formula = cfu_log ~ synID, group.by = c("dpi", "taxa"), ref.group = "C")

read.csv("results/cfu_data_processed.csv") %>% 
    na.omit %>% 
    filter(taxa == "Methylobacterium" & dpi=="07dpi") %>% 
    ggplot(aes(synID, cfu_log, fill=synID))+

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
    
    ggsignif::geom_signif(
        comparisons = dunnet$contrast,
        annotations = dunnet$adj.p.value,
        y_position = c(8, 9)
    )


data = read.csv("results/cfu_data_processed.csv") %>% 
    mutate(synID = factor(synID), dpi = factor(dpi), taxa = factor(taxa)) 
comp = aov(cfu_log ~ synID:taxa + synID:dpi, data = data)
comp %>% summary
comp %>% resid %>% qqnorm

tmp <- expand.grid(synID = unique(data$synID),
                   taxa = unique(data$taxa))
X <- model.matrix(~ synID * taxa, data = tmp)
glht(comp, linfct = X)

Tukey <- contrMat(table(data$synID), "Tukey")
K1 <- cbind(Tukey, matrix(0, nrow=nrow(Tukey), ncol = ncol(Tukey)))
rownames(K1) <- paste(levels(data$taxa)[1], rownames(K1), sep = ":")
K2 <- cbind(matrix(0, nrow=nrow(Tukey), ncol = ncol(Tukey)), Tukey)
rownames(K2) <- paste(levels(data$taxa)[2], rownames(K2), sep = ":")
K <- rbind(K1, K2)
colnames(K) <- c(colnames(Tukey), colnames(Tukey))

summary(glht(comp, linfct = K %*% X), test=adjusted("bonferroni"))

data$tw <- with(data, interaction(synID, taxa))
cell <- lm(cfu_log ~ tw - 1, data)
summary(glht(cell, linfct = K), test=adjusted("bonferroni"))

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
                alpha = 0.1)+
    
    geom_boxplot(
        fill="white",
        width=0.2,
        outlier.alpha = 0
    )+
    
    coord_cartesian(xlim=c(1.2,2.2))+
    
    stat_compare_means(aes(label = after_stat(p.signif)), 
                       method = "t.test", 
                       comparisons = list(c("07dpi", "14dpi")),
                       size = 8,
                       bracket.size = 0.5,
                       symnum.args = list(cutpoints = c(0,0.05)))+
    
    scale_y_continuous(name = plt_bac_density_lab, limits = c(2,11), breaks=c(2,4,6,8,10))+
    scale_x_discrete(name = plt_days_lab, labels = dpi.lab)+
    theme_rs()

read.csv("results/cfu_data_summary.csv") %>% 
    na.omit %>% 
    ggplot(aes(x=cv, y=mean, fill=synID))+
    facet_wrap(~taxa, ncol=2)+
    
    geom_point(pch=21, 
               alpha = 0.75, 
               size = 2)+
    
    scale_y_continuous(name = bquote('Mean population density ['*log[10]~ "CFU g" ~ FW^-1 *"]"), 
                       limits = c(4,9), expand=c(0,0), breaks=c(4,6,8))+
    scale_x_continuous(name = "Coefficient of Variation [%]", 
                       limits = c(0,30), expand=c(0,0), breaks=seq(0,30,10))+
    
    scale_fill_manual(name = "SynCom", values = syn.pal, labels = syn.lab)+
    guides(fill=guide_legend(override.aes = list( shape = 22, size = 5, alpha = 1)))+
    
    theme_rs()
    

