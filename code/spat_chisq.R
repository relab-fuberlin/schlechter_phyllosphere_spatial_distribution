
set.seed(12345)

# Chi-square test
# Set variables for creating the hyperframe in a two-dimensional plane
xrange <- c(0,124.94)
yrange <- c(0,100.24)
unit <- "micron"
W <- owin(xrange, yrange, unitname=unit)
formula_syncom <- "~syncom + dpi + exp + img"

# Read coordinates file
coordinates <- readRDS(here('results', 'coordinates.rds'))
sub_coord <- coordinates %>% 
    group_by(exp, dpi, comID, strain, channel, img) %>% 
    slice_sample(n=100) %>% 
    ungroup
sub_coord$channel <- factor(sub_coord$channel)

H <- hyperframe_function(coordinates, xrange, yrange, unit, formula_syncom)
H <- H[H$n > 10]

fit0 <- mppm(coord~1, H)
fit1 <- mppm(coord~id, H)

saveRDS(fit0, "chi-sq-model0.rds")
saveRDS(fit1, "chi-sq-test.rds")

q <- quadrat.test(fit1, nx=2)

##  Read-ready
chisq = readRDS("chi-sq-test.rds")

quadrat.test(chisq, nx=2)