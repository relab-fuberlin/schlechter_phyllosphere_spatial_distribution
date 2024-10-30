## Plots

##  Plot 5 - Bacterial cell density in the arabidopsis phyllosphere

### p-values 
p_value_plt5a = data_cell %>% 
    mutate(synID = fct_rev(synID)) %>% 
    group_by(taxa, dpi) %>% 
    dunn_test(logCell ~ synID, p.adjust.method = "holm") %>% 
    mutate(p.adj.signif = ifelse(p.adj < 0.05, "*", "ns")) %>% 
    add_xy_position(x = "taxa", 
                    step.increase = 0.5,
                    dodge = 0.6) %>% 
    filter(p.adj.signif != "ns")

p_value_plt5b <- fc_cell_taxa %>% 
    mutate(synID = fct_rev(synID)) %>% 
    group_by(dpi, taxa) %>% 
    dunn_test(log2FC ~ synID, p.adjust.method = "holm") %>% 
    mutate(p.adj.signif = ifelse(p.adj < 0.05, "*", "ns")) %>% 
    add_xy_position(x = "taxa", 
                    step.increase = 0.3,
                    dodge = 0.6) %>% 
    filter(p.adj.signif != "ns")

p_value_plt5c <- fc_cell_strain %>% 
    mutate(synID = fct_rev(synID)) %>% 
    group_by(dpi, strain) %>% 
    dunn_test(log2FC ~ synID, p.adjust.method = "holm") %>% 
    mutate(p.adj.signif = ifelse(p.adj < 0.05, "*", "ns")) %>% 
    add_xy_position(x = "strain", 
                    step.increase = 0.15,
                    dodge = 0.6) %>% 
    filter(p.adj.signif != "ns")

### plot
plt5a <- data_cell %>% 
    ggplot(aes(taxa, logCell))+
    facet_grid(cols = vars(dpi), labeller = labeller(dpi = dpi.lab2))+
    geom_jitter(aes(fill = fct_rev(synID)), alpha = 0.5, size = 0.5, 
                position = position_jitterdodge(jitter.width = 0.3, 
                                                dodge.width = 0.6))+
    geom_boxplot(aes(fill = fct_rev(synID)), 
                 outlier.alpha = 0, alpha = 0.9, size = 0.2, width = 0.6)+
    add_pvalue(data = p_value_plt5a, 
               label = "p.adj.signif", xmin = "xmax", xmax = "xmin", 
               tip.length = 0.005, size = 2,
               lineend = "round", bracket.size = 0.2, coord.flip = TRUE)+
    scale_y_continuous(name = plt_bac_cell_density_lab, 
                       limits = c(2, 8), breaks = c(2, 4, 6, 8))+
    scale_x_discrete(name = "", labels = taxa.lab)+
    scale_fill_manual(name = "SynCom", values = syn.pal, labels = syn.lab, 
                      limits = c("C", "S2", "S3"))+
    theme_rs()+
    coord_flip()+
    theme(axis.text.y = element_text(hjust = 0.5, vjust = 1.5, face = "italic"),
          strip.text.x = element_text(face = "plain"))+
    guides(fill = guide_legend(title.position="top", title.hjust = 0.5, barheight = unit(0.1, 'in')))

plt5b <- fc_cell_taxa %>% 
    ggplot(aes(x = taxa, y = log2FC))+
    facet_grid(cols = vars(dpi), 
               labeller = labeller(dpi = dpi.lab2))+
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3)+
    geom_jitter(aes(fill = fct_rev(synID)), 
                alpha = 0.5, size = 0.5, 
                position = position_jitterdodge(jitter.width = 0.3, 
                                                dodge.width = 0.6))+
    geom_boxplot(aes(fill = fct_rev(synID)), 
                 outlier.alpha = 0, alpha = 0.9, size = 0.2, width = 0.6)+
    add_pvalue(p_value_plt5b, 
               label = "p.adj.signif", xmin = "xmin", xmax = "xmax", 
               tip.length = 0.005, size = 2, lineend = "round", 
               bracket.size = 0.2, coord.flip = TRUE)+
    coord_flip()+
    scale_x_discrete(name = "", labels = taxa.lab)+
    scale_y_continuous(name = bquote(Log[2]~"FC Bacterial cell density"), 
                       limits = c(-10, 5))+
    scale_fill_manual(name = "SynCom", values = syn.pal, labels = syn.lab,
                      limits = c("C", "S2", "S3"))+
    theme_rs()+
    guides(fill = "none")+
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 3),
          axis.text.y = element_text(face = "italic"),
          strip.text = element_text(face = "plain"),
          legend.position = "bottom")

plt5c <- fc_cell_strain %>% 
    ggplot(aes(x = strain, y = log2FC))+
    facet_grid(cols = vars(dpi), 
               labeller = labeller(dpi = dpi.lab2))+
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3)+
    geom_jitter(aes(fill = fct_rev(synID)), 
                alpha = 0.5, size = 0.5, 
                position = position_jitterdodge(jitter.width = 0.3, 
                                                dodge.width = 0.6))+
    geom_boxplot(aes(fill = fct_rev(synID)), 
                 outlier.alpha = 0, alpha = 0.9, size = 0.2, width = 0.6)+
    add_pvalue(p_value_plt5c, 
               label = "p.adj.signif", xmin = "xmin", xmax = "xmax", 
               tip.length = 0.005, size = 2, lineend = "round", 
               bracket.size = 0.2, coord.flip = TRUE)+
    coord_flip()+
    scale_x_discrete(name = "", labels = sp.lab)+
    scale_y_continuous(name = bquote(Log[2]~"FC Bacterial cell density"), 
                       limits = c(-10, 5))+
    scale_fill_manual(name = "SynCom", values = syn.pal, labels = syn.lab,
                      limits = c("C", "S2", "S3"))+
    theme_rs()+
    guides(fill = "none")+
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 3),
          axis.text.y = element_text(face = "plain"),
          strip.text = element_text(face = "plain"),
          legend.position = "bottom")

## wrap plots
wrap_plots(plt5a, plt5b + plt5c, ncol = 1)+
    plot_annotation(tag_levels = "A") + 
    plot_layout(guides = "collect",
                heights = c(1,1)) &
    theme(legend.box.just = "left",
          legend.position = "bottom",
          panel.spacing.x = unit(0.5, "lines"),
          plot.margin = margin(0,1,0,1),
          plot.tag = element_text(size = 7))

# Save plot 
ggsave(here("results", "revision", "fig5.pdf"), width = 6, height = 5)
