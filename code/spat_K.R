

# Set seed
set.seed(19900725)

# Load data
H <- readRDS(here('results', 'hyperframe_strain.rds'))
H_summary <- readRDS(here('results', 'hyperframe_strain_summary.rds'))

##Kest
H_summary = H_summary[H_summary$n > 10,]
H = H[H$n > 10,]

hist(with(H, npoints(coord)))

K_inhom <- with(H, envelope(coord, Kinhom, correction ="Ripley", r=seq(0,30,0.2)))
saveRDS(K_inhom, here('results', 'stat_K_inhom.rds'))

# Create data frame
datalist = list()
for (i in 1:nrow(H)){
    datalist[[i]] <- data.frame(
        syncom = H_summary$syncom[i],
        dpi = H_summary$dpi[i],
        exp = H_summary$exp[i],
        rep = H_summary$img[i],
        strain = H_summary$strain[i],
        r = K_inhom[[i]]$r,
        obs = K_inhom[[i]]$obs,
        lo = K_inhom[[i]]$lo,
        hi = K_inhom[[i]]$hi)
}
K_df = data.table::rbindlist(datalist) %>% na.omit
saveRDS(K_df, here('results', 'stat_K_inhom_table.rds'))

