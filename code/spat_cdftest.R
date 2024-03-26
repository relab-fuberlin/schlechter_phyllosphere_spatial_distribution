#!/usr/bin/env Rscript

# Test for spatial homogeneity
# Load the required libraries
library(spatstat)
library(tidyverse)
library(here)

#  Set seed
set.seed(19900725)

# Chi-square test
# Set variables for creating the hyperframe in a two-dimensional plane
xrange <- c(0, 124.94)
yrange <- c(0, 100.24)
unit <- "micron"
W <- owin(xrange, yrange, unitname=unit)
formula_syncom <- "~syncom + dpi + exp + img"

# Read coordinates file and create a hyperframe
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

# Statistical test
# Fit spatial data into a point procees model
fit1 <- mppm(coord ~ id, subsample_H)

# Run goodness-of-fit of a point process model
cdf_fit1 <- cdf.test(fit1, "x", test = "ks")
write_rds(cdf_fit1, here('results', 'cdf_fit1.rds'))

# Run dispersion test for spatial point pattern based on quadrat counts
qt_fit1 <- quadrat.test(fit1, method = "MonteCarlo", nsim = 999)
write_rds(qt_fit1, here('results', 'quadrat_test_fit1.rds'))
