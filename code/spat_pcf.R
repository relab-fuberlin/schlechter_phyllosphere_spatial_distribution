##  Pair correlation analysis

#   Set seed
set.seed(19900725)

#   Dependencies
library(spatstat)
library(tidyverse)
source(here('code', 'function_pcf.R')) # Custom PCF function

#   Data
H <- readRDS(here('results', 'hyperframe_syncom.rds')) # Read hyperframe

#   PCF for S2
H_S2 <- H[grep("S2", H$syncom)] # Filter by S2 communities. C and S3 are excluded
H_S2 <- H_S2[H_S2$C0 > 10 & H_S2$C1 > 10] # Filter by FOV with at least 10 cells per channel

pcf_S2_C0C1 <- pcf_function(H_S2, i = "C0", j = "C1")

pcf_S2 = list(NULL)
for(i in 1:nrow(H_S2)){
    tryCatch({
        print(paste("Analysing: ", i, H_S2$syncom[[i]]))
        pcf_S2[[i]] <- envelope(H_S2$coord[[i]], pcfcross.inhom, i = "C0", j = "C1", divisor = "d", correction = "isotropic", r=seq(0,30,0.2))
        pcf_S2[[i]] <- pcf_S2[[i]] %>% 
            tibble %>% 
            mutate(pair = "c0c1") %>% 
            cbind(., H_S2[i,-1], row.names=NULL)
    }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}



#   PCF for S3
H_S3 <- H[H$C0 > 10 & H$C1 >10 & H$C2>10,] # Filter by S3 and FOV with at least 10 cells per channel
pcf_S3_C0C1 <- pcf_function(H_S3, i = "C0", j = "C1")
pcf_S3_C0C2 <- pcf_function(H_S3, i = "C0", j = "C2")
pcf_S3_C1C2 <- pcf_function(H_S3, i = "C1", j = "C2")


pcf_S3_C0C1 <- with(H_S3, envelope(coord, pcfcross.inhom, i = "C0", j = "C1", divisor = "d", correction = "isotropic", r=seq(0.2,30,0.2)))
pcf_S3_C0C2 <- with(H_S3, envelope(coord, pcfcross.inhom, i = "C0", j = "C2", divisor = "d", correction = "isotropic", r=seq(0,30,0.2)))
pcf_S3_C1C2 <- with(H_S3, envelope(coord, pcfcross.inhom, i = "C1", j = "C2", divisor = "d", correction = "isotropic", r=seq(0,30,0.2)))


#   Combine results of PCF for S2 and S3 communities
pcf <- rbind(pcf_S2_C0C1, pcf_S3_C0C1, pcf_S3_C0C2, pcf_S3_C1C2)

#   Save data as RDS
write_rds(pcf, here('results', 'pcf.rds'))
