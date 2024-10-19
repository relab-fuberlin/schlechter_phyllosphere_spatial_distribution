##  Palette SynCom spatial distribution paper
library(RColorBrewer)
library(wesanderson)

##    labels
sp.lab = c("MeL85", "MeL92", "Mr0-1", "SmFR1", "SpFA2")
names(sp.lab) <- c("meL85", "meL92", "mr01", "smfr1", "spfa2")

taxa.lab = c("Methylobacterium", "Sphingomonas")
names(taxa.lab) <- taxa.lab

syn.lab = c("C", "S2", "S3")
names(syn.lab) <- syn.lab

dpi.lab = c("7", "14")
names(dpi.lab) <- c("07dpi", "14dpi")

dpi.lab2 <- c("7 dpi", "14 dpi")
names(dpi.lab2) <- c("07dpi", "14dpi")

pattern.lab <- c("Aggregated", "Random", "Regular")
names(pattern.lab) <- c("aggregate_fraction", "random_fraction", "regular_fraction")

unique_pair <- c("meL85.smfr1", "meL85.spfa2", "meL85.meL92", "meL85.mr01", "meL92.smfr1", "meL92.spfa2", "meL92.mr01", "mr01.smfr1", "mr01.spfa2", "smfr1.spfa2")

pair.lab <- c("MeL85-SmFR1", "MeL85-SpFA2", "MeL85-MeL92", "MeL85-Mr0-1", "MeL92-SmFR1", "MeL92-SpFA2", "MeL92-Mr0-1", "Mr0-1-SmFR1", "Mr0-1-SpFA2", "SmFR1-SpFA2")
names(pair.lab) <- unique_pair

unique_pair2 <- c("meL85_smfr1", "meL85_spfa2", "meL85_meL92", "meL85_mr01", "meL92_smfr1", "meL92_spfa2", "meL92_mr01", "mr01_smfr1", "mr01_spfa2", "smfr1_spfa2")

pair.lab2 <- c("MeL85-SmFR1", "MeL85-SpFA2", "MeL85-MeL92", "MeL85-Mr0-1", "MeL92-SmFR1", "MeL92-SpFA2", "MeL92-Mr0-1", "Mr0-1-SmFR1", "Mr0-1-SpFA2", "SmFR1-SpFA2")
names(pair.lab2) <- unique_pair2

taxapair.lab <- c("MM", "MS", "SS")
names(taxapair.lab) <- taxapair.lab

plt_bac_density_lab = bquote('Bacterial density ['*log[10]~ "CFU g" ~ FW^-1 *"]")
plt_bac_cell_density_lab = bquote('Bacterial cell density ['*log[10]~ "cell" ~ cm^-2 *"]")
plt_days_lab = "Days post-inoculation [dpi]"

##    PALETTES
sp.pal = brewer.pal(5, "RdYlBu")
taxa.pal = brewer.pal(11, "RdYlBu")[c(1,11)]
syn.pal = wes_palette("Cavalcanti1",3)
pattern.pal = brewer.pal(3, "PRGn")
taxapair.pal = brewer.pal(11, "Spectral")[c(2,7,10)]

##    DATA FRAME
palette = data.frame(
  type = c(rep("sp", 5), rep("taxa", 2), rep("syn", 3), rep("spatial", 3)),
  name = c(sp.lab, taxa.lab, syn.lab, pattern.lab),
  color = c(sp.pal, taxa.pal, syn.pal, pattern.pal))
