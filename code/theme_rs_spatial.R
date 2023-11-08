#     Ggplot theme
theme_rs <- function() {
  theme_bw() %+replace%
    theme(
      # STRIP
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(size = 6, vjust = 1, face="italic"),
      # LEGEND
      legend.box.spacing = unit(0.1, "line"),
      legend.key.size = unit(0.5, "cm"),
      legend.box.just = "center",
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 5),
      # PANEL
      panel.background = element_blank(),
      panel.border = element_rect(linewidth = 0.5, fill = NA, color = "black"),
      panel.spacing = grid::unit(1.5, "line"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # TEXT GENERAL 
      text = element_text(color = "black"),
      # AXIS
      axis.ticks = element_blank(),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6)
    )
}
