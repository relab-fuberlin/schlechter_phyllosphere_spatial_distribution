hyperframe_function <- function(df, xrange, yrange, unit, formula){
    #
    # df: dataframe containing coordinates and metadata
    # xrange, yrange: each a vector of two values c(x0,x1) and c(y0,y1)
    # unit: str character indicating the unit of the two-dimensional plane
    # formula: str chr of the form "~var1+var2+varn" indicating the grouping factors
    #
    W <- owin(xrange, yrange, unitname=unit)
    
    clist <- split(df, as.formula(formula), sep="_")
    clist <- map(clist, na.omit)
    clist <- clist[sapply(clist, function(x) dim(x)[1]) > 0]
    clist <- map(clist, dplyr::select, c(x, y, channel))
    plist <- map(clist, as.ppp, W)
    plist_unique <- map(plist, unique.ppp)
    H <- hyperframe(coord = plist_unique)
    row.names(H) <- names(plist_unique)
    
    return(H)
}
