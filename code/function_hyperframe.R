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
