##  Palette SynCom spatial distribution paper

##    labels
sp.lab = c("MeL85", "MeL92", "Mr0-1", "SmFR1", "SpFA2")
names(sp.lab) <- c("meL85", "meL92", "mr01", "smfr1", "spfa2")

taxa.lab = c("Methylobacterium", "Sphingomonas")
names(taxa.lab) <- c("Methylobacterium", "Sphingomonas")

syn.lab = c("C", "S2", "S3")
names(syn.lab) <- c("C", "S2", "S3")

dpi.lab = c("7", "14")
names(dpi.lab) <- c("07dpi", "14dpi")
dpi.lab2 = c("7 dpi", "14 dpi")
names(dpi.lab2) <- c("07dpi", "14dpi")

pattern.lab = c("Aggregated", "Random", "Regular")
names(pattern.lab) <- c("aggregate_fraction", "random_fraction", "regular_fraction")

unique_pair <- c("meL85.smfr1", "meL85.spfa2", "meL85.meL92", "meL85.mr01", "meL92.smfr1", "meL92.spfa2", "meL92.mr01", "mr01.smfr1", "mr01.spfa2", "smfr1.spfa2")

pair.lab = c("MeL85-SmFR1", "MeL85-SpFA2", "MeL85-MeL92", "MeL85-Mr0-1", "MeL92-SmFR1", "MeL92-SpFA2", "MeL92-Mr0-1", "Mr0-1-SmFR1", "Mr0-1-SpFA2", "SmFR1-SpFA2")
names(pair.lab) = unique_pair

plt_bac_density_lab = bquote('Bacterial density ['*log[10]~ "CFU g" ~ FW^-1 *"]")
plt_bac_cell_density_lab = bquote('Bacterial cell density ['*log[10]~ "cell" ~ cm^-2 *"]")
plt_days_lab = "Days post-inoculation [dpi]"

##    PALETTES
sp.pal = brewer.pal(5, "RdYlBu")
taxa.pal = brewer.pal(11, "RdYlBu")[c(1,11)]
syn.pal = brewer.pal(5, "YlGnBu")[c(1,2,4)]
syn.pal = wes_palette("Cavalcanti1",3)
pattern.pal = brewer.pal(3, "PRGn")

##    DATA FRAME
palette = data.frame(
  type = c(rep("sp", 5), rep("taxa", 2), rep("syn", 3), rep("spatial", 3)),
  name = c(sp.lab, taxa.lab, syn.lab, pattern.lab),
  color = c(sp.pal, taxa.pal, syn.pal, pattern.pal))
