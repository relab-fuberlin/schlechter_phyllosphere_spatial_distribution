## Plots

##  Plot 6 - Intraspecific spatial patterns

### p-values
w2_taxa <-  K_auc_aggregation %>% 
    mutate(synID = fct_rev(synID)) %>% 
    group_by(synID, dpi) %>% 
    wilcox_test(fractional_change ~ taxa, detailed = TRUE) %>% 
    add_y_position() %>% 
    mutate(p.adj.signif = ifelse(p < 0.05, "*", "ns"),
           y.position = y.position * 100) %>% 
    filter(p.adj.signif != "ns")

w2_strain  <-  K_auc_aggregation %>% 
    group_by(synID, dpi) %>% 
    dunn_test(fractional_change ~ strain, p.adjust.method = "holm") %>% 
    add_y_position(step.increase = 0.03) %>% 
    mutate(p.adj.signif = ifelse(p.adj < 0.05, "*", "ns"),
           y.position = y.position * 100 + 10) %>% 
    filter(p.adj.signif != "ns")

### Plots
plt6a <- K_fractions %>% 
    group_by(taxa, synID, dpi, type, r) %>% 
    summarise(mean_fraction = mean(fraction), .groups="drop") %>% 
    mutate(index = case_when(
        type == "regular_fraction" ~ 1,
        type == "random_fraction" ~ 2,
        type == "aggregate_fraction" ~ 3)) %>% 
    ggplot(aes(x=r, y=mean_fraction))+
    facet_nested(fct_relevel(taxa, c("Sphingomonas", "Methylobacterium")) ~ dpi + synID, 
                 nest_line = element_line(color="black"),
                 labeller = labeller(dpi = dpi.lab2))+
    geom_area(aes(fill = fct_reorder(type, index)), color = "black", linewidth = 0.25)+
    scale_x_continuous(name = expression(paste("Distance, ", italic(r), " (", mu,"m)")), 
                       expand = c(0,0), 
                       limit=c(0.2,30), 
                       breaks = seq(5, 25, 10))+
    scale_y_continuous(name = bquote("Relative frequency of"~hat(K)(italic(r))), 
                       expand = c(0,0), 
                       breaks = seq(0,1,0.5))+
    scale_fill_manual(name = "Spatial Pattern", labels = pattern.lab, values = pattern.pal)+
    theme_rs()+
    theme(panel.spacing.x = grid::unit(0.25, "line"),
          panel.spacing.y = grid::unit(0.5, "line"),
          strip.text.x = element_text(margin = margin(b = 1, t = 5), face = "plain"),
          strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0),
          legend.position = "bottom")

plt6b <- K_auc_aggregation %>% 
    ggplot(aes(x = taxa, y = fractional_change*100))+
    facet_grid(cols = vars(dpi), rows = vars(synID),
               labeller = labeller(dpi = dpi.lab2))+
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3)+
    geom_jitter(aes(fill = fct_rev(synID)), 
                alpha = 0.5, size = 0.5, 
                position = position_jitterdodge(jitter.width = 0.2, 
                                                dodge.width = 0.8))+
    geom_boxplot(aes(fill = fct_rev(synID)),
                 outlier.alpha = 0, alpha = 0.9, size = 0.2, width = 0.4)+
    add_pvalue(w2_taxa, 
               xmin = "group1", xmax = "group2", label = "p.adj.signif",
               tip.length = 0.005, size = 2, lineend = "round", 
               bracket.size = 0.2, coord.flip = TRUE)+
    coord_flip()+
    scale_x_discrete(name = "", labels = taxa.lab)+
    scale_y_continuous(name = "Change in spatial aggregation [%]", 
                       limits = c(-100, 40))+
    scale_fill_manual(name = "SynCom", values = syn.pal[2:3], breaks = c("S2", "S3"))+
    theme_rs()+
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 3),
          axis.text.y = element_text(face = "italic"),
          strip.text.x = element_text(face = "plain"),
          strip.text.y = element_blank(),
          legend.position = "bottom")

plt6c <- K_auc_aggregation %>% 
    ggplot(aes(x = strain, y = fractional_change * 100, fill = synID))+
    facet_grid(cols = vars(dpi), rows = vars(synID),
               labeller = labeller(dpi = dpi.lab2))+
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3)+
    geom_jitter(aes(fill = fct_rev(synID)), 
                alpha = 0.5, size = 0.5, 
                position = position_jitterdodge(jitter.width = 0.2, 
                                                dodge.width = 0.8))+
    geom_boxplot(aes(fill = fct_rev(synID)),
                 outlier.alpha = 0, alpha = 0.9, size = 0.2, width = 0.8)+
    add_pvalue(w2_strain, 
               xmin = "group1", xmax = "group2", label = "p.adj.signif",
               tip.length = 0.005, size = 2, lineend = "round", 
               bracket.size = 0.2, coord.flip = TRUE)+
    coord_flip()+
    scale_x_discrete(name = "", labels = sp.lab)+
    scale_y_continuous(name = "Change in spatial aggregation [%]", 
                       limits = c(-100, 40))+
    scale_fill_manual(name = "SynCom", values = syn.pal[2:3], breaks = c("S2", "S3"))+
    theme_rs()+
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 3),
          axis.text.y = element_text(face = "plain"),
          strip.text.x = element_text(face = "plain"),
          strip.text.y = element_text(face = "plain", angle = 0, vjust = 0.5, hjust = 0))

### Wrap plots
wrap_plots(plt6a, 
           wrap_plots(plt6b, plt6c,
                      ncol = 2,
                      guides = "collect"), 
           ncol = 1) +
    plot_annotation(tag_levels = "A") + 
    plot_layout(heights = c(1,1.2)) &
    theme(panel.spacing.x = unit(0.5, "lines"),
          panel.spacing.y = unit(0.5, "lines"),
          plot.margin = margin(1,2,0,2),
          plot.tag = element_text(size = 7),
          legend.position = "bottom")

### Save plot
ggsave(here("results", "fig6.pdf"), width = 6, height = 5)
