## Plots

##  Plot S7 - Change in interspecific spatial patterns

plts7 <- G_percentage %>% 
    ggplot(aes(x = type, y = percentage, fill = synID))+
    facet_grid(cols = vars(dpi),
               labeller = labeller(dpi = dpi.lab2))+
    geom_hline(yintercept = 1, linetype = 2, linewidth = 0.3)+
    geom_jitter(alpha = 0.5, size = 0.5,
                position = position_jitterdodge(jitter.width = 0.2,
                                                dodge.width = 0.6))+
    geom_boxplot(outlier.alpha = 0, alpha = 0.9, size = 0.2, width = 0.6)+
    scale_x_discrete(name = "Spatial pattern", labels = pattern.lab)+
    scale_y_continuous(name = "Relative contribution [%]",
                       limits = c(0, 100))+
    scale_fill_manual(name = "SynCom", values = syn.pal, labels = syn.lab,
                      limits = c("C", "S2", "S3"))+
    theme_rs()+
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 3),
          axis.text.y = element_text(face = "plain"),
          strip.text.x = element_text(face = "plain"),
          strip.text.y = element_blank(),
          legend.position = "bottom")

# Save plot
ggsave(here("results", "figS9.pdf"), width = 4, height = 2)    
