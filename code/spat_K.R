
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
