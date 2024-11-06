## Plots

##  Plot 2 - Bacterial population density in the arabidopsis phyllosphere

### p-values 
p_value_plt2a = data_cfu %>% 
    mutate(synID = fct_rev(synID)) %>% 
    group_by(taxa, dpi) %>% 
    dunn_test(cfu_log ~ synID, p.adjust.method = "holm") %>% 
    mutate(p.adj.signif = ifelse(p.adj < 0.05, "*", "ns")) %>% 
    add_xy_position(x = "taxa", 
                    step.increase = 0.1,
                    dodge = 0.6) %>% 
    filter(p.adj.signif != "ns")

p_value_plt2b <- fc_cfu_taxa %>% 
    mutate(synID = fct_rev(synID)) %>% 
    group_by(dpi, synID) %>% 
    dunn_test(log2FC ~ taxa, p.adjust.method = "holm") %>% 
    mutate(p.adj.signif = ifelse(p.adj < 0.05, "*", "ns")) %>% 
    add_xy_position(x = "taxa", 
                    step.increase = 0.1,
                    dodge = 0.6) %>% 
    filter(p.adj.signif != "ns")

p_value_plt2c <- fc_cfu_strain %>% 
    mutate(synID = fct_rev(synID),
           strain = fct_rev(strain)) %>% 
    group_by(dpi, strain) %>% 
    dunn_test(log2FC ~ synID, p.adjust.method = "holm") %>% 
    add_xy_position(x = "strain",
                    step.increase = 0.075) %>% 
    mutate(p.adj.signif = ifelse(p.adj < 0.05, "*", "ns")) %>%
    filter(p.adj.signif != "ns")

### plot
plt2a <- data_cfu %>% 
    ggplot(aes(taxa, cfu_log))+
    facet_grid(cols = vars(dpi), labeller = labeller(dpi = dpi.lab2))+
    geom_jitter(aes(fill = fct_rev(synID)), alpha = 0.5, size = 0.5, 
                position = position_jitterdodge(jitter.width = 0.3, 
                                                dodge.width = 0.6))+
    geom_boxplot(aes(fill = fct_rev(synID)), 
                 outlier.alpha = 0, alpha = 0.9, size = 0.2, width = 0.6)+
    add_pvalue(data = p_value_plt2a, 
               label = "p.adj.signif", xmin = "xmax", xmax = "xmin", 
               tip.length = 0.005, size = 2,
               lineend = "round", bracket.size = 0.2, coord.flip = TRUE)+
    scale_y_continuous(name = plt_bac_density_lab, 
                       limits = c(3.5, 10.5), breaks = c(4, 6, 8, 10))+
    scale_x_discrete(name = "", labels = taxa.lab)+
    scale_fill_manual(name = "SynCom", values = syn.pal, labels = syn.lab, 
                      limits = c("C", "S2", "S3"))+
    theme_rs()+
    coord_flip()+
    theme(axis.text.y = element_text(hjust = 0.5, vjust = 1.5, face = "italic"),
          strip.text.x = element_text(face = "plain"))+
    guides(fill = guide_legend(title.position="top", title.hjust = 0.5, barheight = unit(0.1, 'in')))

plt2b <- fc_cfu_taxa %>% 
    ggplot(aes(x = taxa, y = log2FC))+
    facet_grid(cols = vars(dpi), rows = vars(synID),
               labeller = labeller(dpi = dpi.lab2))+
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3)+
    geom_jitter(aes(fill = fct_rev(synID)), 
                alpha = 0.5, size = 0.5, 
                position = position_jitterdodge(jitter.width = 0.3, 
                                                dodge.width = 0.6))+
    geom_boxplot(aes(fill = fct_rev(synID)), 
                 outlier.alpha = 0, alpha = 0.9, size = 0.2, width = 0.6)+
    add_pvalue(p_value_plt2b, 
               label = "p.adj.signif", xmin = "xmin", xmax = "xmax", 
               tip.length = 0.005, size = 2, lineend = "round", 
               bracket.size = 0.2, coord.flip = TRUE)+
    coord_flip()+
    scale_x_discrete(name = "", labels = taxa.lab)+
    scale_y_continuous(name = bquote(Log[2]~"FC Bacterial density"), 
                       limits = c(-10, 9))+
    scale_fill_manual(name = "SynCom", values = syn.pal, labels = syn.lab,
                      limits = c("C", "S2", "S3"))+
    theme_rs()+
    guides(fill = "none")+
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 3),
          axis.text.y = element_text(face = "italic"),
          strip.text.x = element_text(face = "plain"),
          strip.text.y = element_blank(),
          legend.position = "bottom")

plt2c <- fc_cfu_strain %>% 
    mutate(strain = fct_relevel(strain, rev)) %>% 
    ggplot(aes(x = fct_rev(synID), y = log2FC))+
    facet_grid(cols = vars(dpi), rows = vars(strain),
               labeller = labeller(dpi = dpi.lab2,
                                   strain = sp.lab),
               switch = "y")+
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3)+
    geom_jitter(aes(fill = fct_rev(synID)), 
                alpha = 0.5, size = 0.5, 
                position = position_jitterdodge(jitter.width = 0.3, 
                                                dodge.width = 0.6))+
    geom_boxplot(aes(fill = fct_rev(synID)), 
                 outlier.alpha = 0, alpha = 0.9, size = 0.2, width = 0.6)+
    add_pvalue(p_value_plt2c, 
               label = "p.adj.signif", xmin = "group1", xmax = "group2", 
               tip.length = 0.005, size = 2, lineend = "round", 
               bracket.size = 0.2, coord.flip = TRUE)+
    coord_flip()+
    scale_x_discrete(name = "")+
    scale_y_continuous(name = bquote(Log[2]~"FC Bacterial density"), 
                       limits = c(-10, 9))+
    scale_fill_manual(name = "SynCom", values = syn.pal, labels = syn.lab,
                      limits = c("C", "S2", "S3"))+
    theme_rs()+
    guides(fill = "none")+    
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 3),
          axis.text.y = element_blank(),
          strip.text.x = element_text(face = "plain"),
          strip.text.y.left = element_text(face = "plain", angle = 0, vjust = 0.5, hjust = 0.5),
          legend.position = "bottom")

## wrap plots
wrap_plots(plt2a, plt2b + plt2c, ncol = 1)+
    plot_annotation(tag_levels = "A") + 
    plot_layout(guides = "collect",
                heights = c(1,1)) &
    theme(legend.box.just = "left",
          legend.position = "bottom",
          panel.spacing = unit(0.5, "lines"),
          plot.margin = margin(0,1,0,1),
          plot.tag = element_text(size = 7))

# Save plot 
ggsave(here("results", "fig2.pdf"), width = 6, height = 5)
