## Plots

##  Plot S5 - Effect of tritagonist on interactions (microscopy data)

### p-values
p_value_pltS5_taxa <- fold_change_S3 %>% 
    group_by(dpi) %>% 
    wilcox_test(log2FC ~ taxa) %>% 
    mutate(p.adj.signif = ifelse(p < 0.05, "*", "ns")) %>% 
    add_xy_position(x = "taxa", 
                    step.increase = 0.1,
                    dodge = 0.6) %>% 
    mutate(y.position = y.position + 4) %>% 
    filter(p.adj.signif != "ns")

p_value_pltS5a <- fold_change_S3 %>% 
    group_by(dpi, taxa) %>% 
    dunn_test(log2FC ~ tritagonist, p.adjust.method = "holm") %>% 
    mutate(p.adj.signif = ifelse(p < 0.05, "*", "ns")) %>% 
    add_xy_position(x = "taxa", 
                    step.increase = 0.1,
                    dodge = 0.8) %>%
    filter(p.adj.signif != "ns")

p_value_pltS5b <- fold_change_S3 %>% 
    group_by(dpi, taxa_pair) %>% 
    dunn_test(log2FC ~ tritagonist, p.adjust.method = "holm") %>% 
    mutate(p.adj.signif = ifelse(p < 0.05, "*", "ns")) %>% 
    add_xy_position(x = "taxa_pair", 
                    step.increase = 0.1,
                    dodge = 0.8) %>%
    filter(p.adj.signif != "ns")

p_value_pltS5c <- fold_change_S3 %>% 
    group_by(dpi, taxa) %>% 
    wilcox_test(log2FC ~ taxa_tritagonist) %>% 
    mutate(p.adj.signif = ifelse(p < 0.05, "*", "ns")) %>% 
    add_xy_position(x = "taxa", 
                    step.increase = 0.1,
                    dodge = 0.6) %>%
    filter(p.adj.signif != "ns")

p_value_pltS5d <- fold_change_S3 %>% 
    group_by(dpi) %>% 
    dunn_test(log2FC ~ taxa_pair) %>% 
    mutate(p.adj.signif = ifelse(p < 0.05, "*", "ns")) %>% 
    add_xy_position(x = "taxa_pair", 
                    step.increase = 0.5,
                    dodge = 0.8) %>% 
    mutate(y.position = y.position) %>% 
    filter(p.adj.signif != "ns")

p_value_pltS5d2 <- fold_change_S3 %>% 
    filter(taxa_pair != "SS") %>% 
    group_by(dpi, taxa_pair) %>% 
    wilcox_test(log2FC ~ taxa_tritagonist) %>% 
    mutate(p.adj.signif = ifelse(p < 0.05, "*", "ns")) %>% 
    add_xy_position(x = "taxa_pair", 
                    step.increase = 0.1,
                    dodge = 0.8)



### plots
pltS5a <- fold_change_S3 %>% 
    ggplot(aes(x = taxa, y = log2FC))+
    facet_grid(cols = vars(dpi), labeller = labeller(dpi = dpi.lab2))+
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3)+
    geom_jitter(aes(fill = tritagonist), 
                alpha = 0.5, size = 0.5, 
                position = position_jitterdodge(jitter.width = 0.05, 
                                                dodge.width = 0.8))+
    geom_boxplot(aes(fill = tritagonist),
                 outlier.alpha = 0, alpha = 0.9, size = 0.2, width = 0.8)+
    add_pvalue(data = p_value_pltS5_taxa, 
               label = "p.adj.signif", xmin = "xmax", xmax = "xmin", 
               tip.length = 0.02, size = 2,
               lineend = "round", bracket.size = 0.4, coord.flip = TRUE)+
    add_pvalue(data = p_value_pltS5a, 
               label = "p.adj.signif", xmin = "xmax", xmax = "xmin", 
               tip.length = 0.005, size = 2,
               lineend = "round", bracket.size = 0.2, coord.flip = TRUE)+
    coord_flip()+
    scale_x_discrete(name = "", labels = taxa.lab)+
    scale_y_continuous(name = bquote(Log[2]~"FC Bacterial cell density"), 
                       limits = c(-13.5, 10))+
    scale_fill_manual(name = "Tritagonist", values = sp.pal, labels = sp.lab)+
    theme_rs()+
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 3),
          axis.text.y = element_text(face = "italic"),
          strip.text = element_text(face = "plain"),
          legend.position = "bottom")

pltS5b <- fold_change_S3 %>% 
    ggplot(aes(x = taxa_pair, y = log2FC))+
    facet_grid(cols = vars(dpi), labeller = labeller(dpi = dpi.lab2))+
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3)+
    geom_jitter(aes(fill = tritagonist),
                alpha = 0.5, size = 0.5, 
                position = position_jitterdodge(jitter.width = 0.05, 
                                                dodge.width = 0.8))+
    geom_boxplot(aes(fill = tritagonist),
                 outlier.alpha = 0, alpha = 0.9, size = 0.2, 
                 position = position_dodge2(width = 0.75, 
                                            preserve = "single"))+
    add_pvalue(data = p_value_pltS5b, 
               label = "p.adj.signif", xmin = "xmax", xmax = "xmin", 
               tip.length = 0.005, size = 2,
               lineend = "round", bracket.size = 0.2, coord.flip = TRUE)+
    coord_flip()+
    scale_x_discrete(name = "", labels = taxa.lab)+
    scale_y_continuous(name = bquote(Log[2]~"FC Bacterial cell density"), 
                       limits = c(-13.5, 10))+
    scale_fill_manual(name = "Tritagonist", values = sp.pal, labels = sp.lab)+
    theme_rs()+
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 3),
          axis.text.y = element_text(face = "plain"),
          strip.text = element_text(face = "plain"),
          legend.position = "bottom")

pltS5c <- fold_change_S3 %>% 
    ggplot(aes(x = taxa, y = log2FC))+
    facet_grid(cols = vars(dpi), labeller = labeller(dpi = dpi.lab2))+
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3)+
    geom_jitter(aes(fill = taxa_tritagonist), 
                alpha = 0.5, size = 0.5, 
                position = position_jitterdodge(jitter.width = 0.3, 
                                                dodge.width = 0.6))+
    geom_boxplot(aes(fill = taxa_tritagonist),
                 outlier.alpha = 0, alpha = 0.9, size = 0.2, width = 0.6)+
    add_pvalue(data = p_value_pltS5_taxa, 
               label = "p.adj.signif", xmin = "xmax", xmax = "xmin", 
               tip.length = 0.02, size = 2,
               lineend = "round", bracket.size = 0.4, coord.flip = TRUE)+
    add_pvalue(data = p_value_pltS5c, 
               label = "p.adj.signif", xmin = "xmax", xmax = "xmin", 
               tip.length = 0.005, size = 2,
               lineend = "round", bracket.size = 0.2, coord.flip = TRUE)+
    coord_flip()+
    scale_x_discrete(name = "", labels = taxa.lab)+
    scale_y_continuous(name = bquote(Log[2]~"FC Bacterial cell density"), 
                       limits = c(-13.5, 10))+
    scale_fill_manual(name = "Tritagonist", values = taxa.pal, labels = taxa.lab)+
    theme_rs()+
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 3),
          axis.text.y = element_text(face = "italic"),
          strip.text = element_text(face = "plain"),
          legend.position = "bottom")

pltS5d <- fold_change_S3 %>% 
    ggplot(aes(x = taxa_pair, y = log2FC))+
    facet_grid(cols = vars(dpi), labeller = labeller(dpi = dpi.lab2))+
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3)+
    geom_jitter(aes(fill = taxa_tritagonist),
                alpha = 0.5, size = 0.5, 
                position = position_jitterdodge(jitter.width = 0.3, 
                                                dodge.width = 0.8))+
    geom_boxplot(aes(fill = taxa_tritagonist),
                 outlier.alpha = 0, alpha = 0.9, size = 0.2, 
                 position = position_dodge2(width = 0.8, 
                                            preserve = "single"))+
    add_pvalue(data = p_value_pltS5d, 
               label = "p.adj.signif", xmin = "xmax", xmax = "xmin", 
               tip.length = 0.005, size = 2,
               lineend = "round", bracket.size = 0.4, coord.flip = TRUE)+
    coord_flip()+
    scale_x_discrete(name = "", labels = taxa.lab)+
    scale_y_continuous(name = bquote(Log[2]~"FC Bacterial cell density"), 
                       limits = c(-13.5, 10))+
    scale_fill_manual(name = "Tritagonist", values = taxa.pal, labels = taxa.lab)+
    theme_rs()+
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 3),
          axis.text.y = element_text(face = "plain"),
          strip.text = element_text(face = "plain"),
          legend.position = "bottom")

wrap_plots(pltS5c + pltS5a + pltS5d + pltS5b,
           ncol = 1)+
    plot_annotation(tag_levels = "A") + 
    plot_layout(guides = "collect",
                widths = c(1,1)) &
    theme(legend.box.just = "left",
          legend.position = "bottom",
          panel.spacing.x = unit(0.5, "lines"),
          plot.margin = margin(1,1,0,1),
          plot.tag = element_text(size = 7)) &
    guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, barheight = unit(0.1, 'in')))

ggsave(here("results", "figS5.pdf"), width = 8, height = 5)


