#!/usr/bin/env Rscript 

library(spatstat)

# Set variables for creating the hyperframe in a two-dimensional plane
xrange <- c(0,124.94)
yrange <- c(0,100.24)
unit <- "micron"
W <- owin(xrange, yrange, unitname=unit)
formula_syncom <- "~syncom + dpi + exp + img"
formula_strain <- "~syncom + dpi + exp + img + strain"
 
# Create hyperframe with marked points for channels
coordinates <- readRDS(here('results', 'coordinates.rds'))

hyperframe_function <- function(df, xrange, yrange, unit, formula){
    #
    # df: dataframe containing coordinates and metadata
    # xrange, yrange: each a vector of two values c(x0,x1) and c(y0,y1)
    # unit: str character indicating the unit of the two-dimensional plane
    # formula: str chr of the form "~var1+var2+varn" indicating the grouping factors
    #
    clist <- split(df, as.formula(formula), sep="_")
    clist <- lapply(clist, na.omit)
    clist <- clist[sapply(clist, function(x) dim(x)[1]) > 0]
    clist <- lapply(clist, dplyr::select, c(x, y, channel))
    list_coord <- lapply(clist, "[", c("x", "y"))
    wlist <- rep(list(W), length(list_coord))
    plist <- mapply(as.ppp, clist, W = wlist, SIMPLIFY = FALSE, multitype=TRUE)
    H <- hyperframe(coord = plist, rep = names(plist))
    
    return(H)
}

H_syncom <- hyperframe_function(coordinates, xrange, yrange, unit, formula_syncom)
H_strain <- hyperframe_function(coordinates, xrange, yrange, unit, formula_strain)

summary_data <- coordinates %>% 
    group_by(syncom, dpi, exp, img, channel) %>% 
    tally() %>% 
    pivot_wider(id_cols=c(syncom, exp, dpi, img), names_from = channel, values_from = n, values_fill = 0) %>% 
    mutate(rep = paste(syncom,dpi,exp,img, sep = "_"))

reord_index <- match(H$rep, summary_data$rep)
summary_data <- summary_data[reord_index,]

#double check with:
match(H$rep, summary_data$rep)

# Add attributes to the hyperframe
H$C0 = summary_data$C0
H$C1 = summary_data$C1
H$C2 = summary_data$C2
H$sum = H$C0 + H$C1 + H$C2
H$syncom = summary_data$syncom
H$dpi = summary_data$dpi

# Export data
saveRDS(H, here('results', 'hyperframe.rds'))
saveRDS(summary_data, here('results', 'hyperframe_summary.rds'))
