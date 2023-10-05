
set.seed(19900725)

# Chi-square test
# Read coordinates file
H <- readRDS(here('results', 'hyperframe_syncom.rds'))

subsample_H <- H[H$sum > 10]
subsample_H <- subsample_H[sample(nrow(subsample_H), size=100, replace = FALSE),]

# test wrong model
fit0 <- mppm(coord ~ 1, subsample_H)
cdf.test(fit0, "x")

# test right model
fit1 <- mppm(coord ~ marks, subsample_H)
cdf.test(fit1, "x")
