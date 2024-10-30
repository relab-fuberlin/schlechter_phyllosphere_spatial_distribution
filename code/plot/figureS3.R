## Plots

##  Plot S3 - Relative abundances of strain within communities
#   create a label
label <- interaction(dpi.lab, relative_fractions$syncom, sep = "\n")

#   plot
pltS3 <- relative_fractions %>% 
    ggplot(., aes(x = interaction(dpi, syncom), y = relative_fraction, fill = strain))+
    geom_col()+
    scale_fill_manual(name = "Strain", values = sp.pal, labels = sp.lab)+
    scale_x_discrete(name = "", labels = label)

wrap_plots(plt_relfraction_strain, ncol = 1)+
    labs(x = "SynCom", y = "Relative fraction") &
    theme_rs() &
    theme(
        legend.box.just = "left",
        legend.position = "bottom",
        panel.spacing.x = unit(0.1, "lines"),
        strip.text = element_text(face = "plain"),
        plot.margin = margin(1,1,0,1),
        plot.tag = element_text(size = 7)) &
    guides(fill = guide_legend(title.position="top", title.hjust = 0.5, barheight = unit(0.1, 'in')))

# Save plot  
ggsave(here("results", "figS3.pdf"), width = 8, height = 3)
