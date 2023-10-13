##  Pair correlation analysis
library(spatstat)

set.seed(19900725)

#   Load data
H <- readRDS(here('results', 'hyperframe_syncom.rds'))

pcf_function <- function(hf, i, j){
    safe_envelope = possibly(.f = envelope, otherwise=NULL)
    pcf <- map(hf$coord, safe_envelope, pcfcross.inhom, i=i, j=j, divisor = "d", correction = "isotropic", r=seq(0,30,0.2))
    pcf_tibble <- map(pcf, tibble) %>% 
        enframe(name = "rep") %>%
        unnest(cols = "value")
    
    return(pcf_tibble)
}

#   PCF for S2
H_S2 <- H[grep("S2", H$syncom)]
H_S2 <- H_S2[H_S2$C0 > 10 & H_S2$C1 > 10]
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


pcf_function(H_S2, i = "C0", j = "C1")


#   PCF for S3
H_S3 <- H[H$C0 > 10 & H$C1 >10 & H$C2>10,]
pcf_S3_C0C1 <- with(H_S3, envelope(coord, pcfcross.inhom, i = "C0", j = "C1", divisor = "d", correction = "isotropic", r=seq(0.2,30,0.2)))
pcf_S3_C0C2 <- with(H_S3, envelope(coord, pcfcross.inhom, i = "C0", j = "C2", divisor = "d", correction = "isotropic", r=seq(0,30,0.2)))
pcf_S3_C1C2 <- with(H_S3, envelope(coord, pcfcross.inhom, i = "C1", j = "C2", divisor = "d", correction = "isotropic", r=seq(0,30,0.2)))

write_rds(pcf_S2, here('results', 'pcf_S2.rds'))
write_rds(pcf_S3_C0C1, here('results', 'pcf_S3_C0C1.rds'))
write_rds(pcf_S3_C0C2, here('results', 'pcf_S3_C0C2.rds'))
write_rds(pcf_S3_C1C2, here('results', 'pcf_S3_C1C2.rds'))
