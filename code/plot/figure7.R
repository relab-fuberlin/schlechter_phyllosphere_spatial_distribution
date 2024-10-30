## Plots

##  Plot 7 - Interspecific spatial patterns

### p-values
pcf_pair_full <- G_auc_fold_change %>% 
    mutate(strain_pair = str_replace(strain_pair, "\\.", "\\_"),
           strain_pair = case_when(strain_pair == "mr01_smfr1" ~ "smfr1_mr01", TRUE ~ strain_pair)) %>% 
    inner_join(., unique(cfu_pair[,1:7]), 
               by = c("syncom", "synID", "dpi", "taxa_pair", "strain_pair" = "pair"))

# One-sample Wilcoxon test
w1_pair <- pcf_pair_full %>% 
    group_by(dpi, taxa_pair, taxa_excluded, type) %>% 
    wilcox_test(fractional_change ~ 1, mu = 0, detailed = TRUE) %>% 
    add_xy_position(x = "taxa_pair", 
                    dodge = 0.75) %>% 
    mutate(p.adj.signif = ifelse(p < 0.05, "*", "ns"),
           y.position = y.position * 100) %>% 
    filter(p.adj.signif != "ns")

### Plots
plt7a <- G_fractions %>% 
    group_by(taxa_pair, synID, dpi, type, r) %>% 
    summarise(mean_fraction = mean(fraction), .groups="drop") %>% 
    mutate(index = case_when(
        type == "regular_fraction" ~ 1,
        type == "random_fraction" ~ 2,
        type == "aggregate_fraction" ~ 3)) %>% 
    ggplot(aes(x = r, y = mean_fraction))+
    facet_nested(taxa_pair ~ dpi + synID, 
                 nest_line = element_line(color="black"),
                 labeller = labeller(dpi = dpi.lab2))+
    geom_area(aes(fill = fct_reorder(type, index)), color = "black", linewidth = 0.25)+
    scale_x_continuous(name = expression(paste("Distance, ", italic(r), " (", mu,"m)")), 
                       expand = c(0,0), limit=c(0.2,30), breaks = seq(5, 25, 10))+
    scale_y_continuous(name = bquote("Relative frequency of"~hat(g)(italic(r))),
                       expand = c(0,0), breaks = seq(0,1,0.5))+
    scale_fill_manual(name = "Spatial Pattern", labels = pattern.lab, values = pattern.pal)+
    theme_rs()+
    theme(panel.spacing.x = grid::unit(0.25, "line"),
          panel.spacing.y = grid::unit(0.5, "line"),
          strip.text.x = element_text(margin = margin(b = 1, t = 5), face = "plain"),
          strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0, face = "plain"),
          legend.position = "bottom")

plt7b <- pcf_pair_full %>% 
    ggplot(aes(x = taxa_pair, y = 100 * fractional_change, fill = taxa_excluded))+
    facet_grid(rows = vars(dpi), cols = vars(fct_rev(type)),
               labeller = labeller(.rows = dpi.lab2, .cols = pattern.lab))+
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3)+
    geom_jitter(aes(fill = taxa_excluded), 
                alpha = 0.5, size = 0.5, 
                position = position_jitterdodge(jitter.width = 0.2, 
                                                dodge.width = 0.8))+
    geom_boxplot(aes(fill = taxa_excluded),
                 outlier.alpha = 0, alpha = 0.9, size = 0.2, width = 0.8,
                 position = position_dodge2(width = 0.45, preserve = "single"))+
    geom_text(data = w1_pair, aes(label = p.adj.signif, y = y.position + 50, group = taxa_excluded), size = 4,
              position = position_dodge2(width = 0.9))+
    theme_rs()+
    scale_y_continuous(name = "Change in spatial pattern relative to S2 [%]", limits = c(-120, 320))+
    scale_x_discrete(name = "Taxa pair")+
    scale_fill_manual(name = "Tritagonist taxon", values = taxa.pal, labels = taxa.lab)+
    coord_flip()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(face = "italic"),
          strip.text = element_text(face = "plain"),
          panel.spacing = unit(0.5, "lines"),
          legend.position = "bottom")

wrap_plots(plt7a, plt7b, ncol = 1)+ 
    plot_annotation(tag_levels = "A") &
    theme(legend.box.just = "center",
          plot.margin = margin(0,1,0,1),
          plot.tag = element_text(size = 7))

# Save plot
ggsave(here("results", "revision", "fig7.pdf"), width = 4, height = 6)

