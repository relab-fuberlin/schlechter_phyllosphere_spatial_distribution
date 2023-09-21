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

pattern.lab = c("Aggregated", "CSR", "Regular")
names(pattern.lab) <- c("f_agg", "f_csr", "f_seg")

pair.lab = c("MeL85-MeL92", "MeL85-Mr0-1", "MeL92-Mr0-1", "SmFR1-SpFA2",
              "MeL85-SmFR1", "MeL85-SpFA2", "MeL92-SmFR1", "MeL92-SpFA2", "Mr0-1-SmFR1", "Mr0-1-SpFA2")
names(pair.lab) = c('meth85_meth92', 'meth85_mr01', 'meth92_mr01','smfr1_spfa2',
                     'meth85_smfr1', 'meth85_spfa2', 'meth92_smfr1', 'meth92_spfa2','mr01_smfr1','mr01_spfa2')

plt_bac_density_lab = bquote('Bacterial density ['*log[10]~ "CFU g" ~ FW^-1 *"]")
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
