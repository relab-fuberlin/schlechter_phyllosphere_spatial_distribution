#     Ggplot theme
theme_rs <- function() {
  theme_bw(base_size = 12) %+replace%
    theme(
      strip.background = element_rect(fill = "#f2f2f2", color = "gray", size = 0.25),
      axis.ticks = element_line(size = 0.2),
      axis.ticks.length=unit(.25, "cm"),
      legend.box.spacing = unit(0.1, "line"),
      legend.key.size = unit(0.1, "cm"),
      panel.background = element_blank(),
      panel.border = element_rect(size = 0.25, fill = NA, color = "black"),
      panel.spacing = grid::unit(0.5, "line"),
      panel.grid.major = element_line(size=0.25, color = "#f2f2f2"),
      panel.grid.minor = element_blank()
    )
}
