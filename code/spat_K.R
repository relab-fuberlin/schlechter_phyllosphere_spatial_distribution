

# Set seed
set.seed(19900725)

# Load data
H <- readRDS(here('results', 'hyperframe.rds'))
H_summary <- readRDS(here('results', 'hyperframe_summary.rds'))

##Kest
H_summary = H_summary[H_summary$C0 > 5 | H_summary$C1 > 5 | H_summary$C2 > 5,]
H = H[H$C0 > 5 | H$C1 > 5 | H$C2 > 5,]

H_summary %>% 
    group_by(syncom, dpi) %>% 
    tally()

hist(with(H, npoints(coord)))

K_inhom <- with(H, envelope(coord, Kinhom, correction ="Ripley", r=seq(0,30,0.2)))
saveRDS(K_inhom, here('results', 'stat_K_inhom.rds'))

# Create data frame
datalist = list()
for (i in 1:nrow(H)){
    datalist[[i]] <- data.frame(
        syncom = summary_data$syncom[i],
        dpi = summary_data$dpi[i],
        rep = summary_data$img[i],
        r = K_inhom[[i]]$r,
        obs = K_inhom[[i]]$obs,
        lo = K_inhom[[i]]$lo,
        hi = K_inhom[[i]]$hi)
}
K_df = data.table::rbindlist(datalist) %>% na.omit

