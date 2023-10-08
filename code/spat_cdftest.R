
library(spatstat)
library(tidyverse)
library(here)

set.seed(19900725)

# Chi-square test
# Set variables for creating the hyperframe in a two-dimensional plane
xrange <- c(0,124.94)
yrange <- c(0,100.24)
unit <- "micron"
W <- owin(xrange, yrange, unitname=unit)
formula_syncom <- "~syncom + dpi + exp + img"

# Read coordinates file
coordinates <- readRDS(here('results', 'coordinates.rds'))
clist <- split(coordinates, as.formula(formula_syncom), sep="_")
clist <- lapply(clist, na.omit)
clist <- clist[sapply(clist, function(x) dim(x)[1]) > 0]
clist <- lapply(clist, dplyr::select, c(x, y))
plist <- lapply(X = clist, FUN = as.ppp, W)
H <- hyperframe(coord = plist, rep = names(plist))
H_syncom_summary <- coordinates %>% 
    group_by(syncom, dpi, exp, img, channel) %>% 
    tally() %>% 
    pivot_wider(id_cols=c(syncom, exp, dpi, img), names_from = channel, values_from = n, values_fill = 0) %>% 
    mutate(rep = paste(syncom,dpi,exp,img, sep = "_"))
reord_index <- match(H$rep, H_syncom_summary$rep)
H_syncom_summary <- H_syncom_summary[reord_index,]

# Add attributes to the hyperframe
H$C0 = H_syncom_summary$C0
H$C1 = H_syncom_summary$C1
H$C2 = H_syncom_summary$C2
H$sum = H$C0 + H$C1 + H$C2
H2 <- H[H$sum > 10]
subsample_H <- H2[sample(nrow(H2), size=500, replace = FALSE),]

# test wrong model
fit0 <- mppm(coord ~ 1, subsample_H)
cdf_fit0 <- cdf.test(fit0, "x")
qt_fit0 <- quadrat.test(fit0, method = "MonteCarlo", nsim=999)
write_rds(qt_fit0, here('results', 'quadrat_test_fit0.rds'))

# test right model
fit1 <- mppm(coord ~ id, subsample_H)
cdf.test(fit1, "x")
qt_fit0 <- quadrat.test(fit0, method = "MonteCarlo", nsim=999)
write_rds(qt_fit1, here('results', 'quadrat_test_fit1.rds'))



H <- hyperframe(X=waterstriders)
# Poisson with constant intensity for all patterns
fit1 <- mppm(X~1, H)
quadrat.test(fit1, nx=2)

# uniform Poisson with different intensity for each pattern
fit2 <- mppm(X ~ id, H)
quadrat.test(fit2, nx=2)
