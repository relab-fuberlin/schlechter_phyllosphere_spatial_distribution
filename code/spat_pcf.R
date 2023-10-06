##  Pair correlation analysis
library(spatstat)

#   Set seed
set.seed(19900725)

#   Load data
H <- readRDS(here('results', 'hyperframe_syncom.rds'))

#   PCF for S2
H_S2 <- H[grep("S2", H$syncom),]
H_S2 <- H_S2[H_S2$C0 > 10 & H_S2$C1 > 10, ]
pcf_S2 <- with(H_S2, envelope(coord, pcfcross.inhom, i = "C0", j = "C1", divisor = "d", correction = "isotropic", r=seq(0,30,0.2)))

#   PCF for S3
H_S3 <- H[H$C0 > 10 & H$C1 >10 & H$C2>10,]
pcf_S3_C0C1 <- with(H_S3, envelope(coord, pcfcross.inhom, i = "C0", j = "C1", divisor = "d", correction = "isotropic", r=seq(0,30,0.2)))
pcf_S3_C0C2 <- with(H_S3, envelope(coord, pcfcross.inhom, i = "C0", j = "C2", divisor = "d", correction = "isotropic", r=seq(0,30,0.2)))
pcf_S3_C1C2 <- with(H_S3, envelope(coord, pcfcross.inhom, i = "C1", j = "C2", divisor = "d", correction = "isotropic", r=seq(0,30,0.2)))
