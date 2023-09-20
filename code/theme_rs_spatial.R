#     Ggplot theme
theme_rs <- function() {
  theme_bw() %+replace%
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(size=12, vjust = 1, face="italic"),
      legend.box.spacing = unit(0.1, "line"),
      legend.key.size = unit(0.1, "cm"),
      panel.background = element_blank(),
      panel.border = element_rect(size = 1, fill = NA, color = "black"),
      panel.spacing = grid::unit(1.5, "line"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      text = element_text(color = "black"),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
}
