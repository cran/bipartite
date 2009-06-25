robustness<- function (object)
{
    # a function to calculate robustness to extinction in the other trophic level
    # based on output by "second.extinct" (see there)
    # by Mariano Devoto, April 2009
    
    if (class(object) != "bipartite")
        stop("This function cannot be meaningfully applied to objects of this class.")
    N <- colSums(object)
    if (all(object[-nrow(object), 2] == 1)) #this IF selects the appropriate column from the data file
        y <- -object[, 3]
    else y <- -object[, 2]
    y <- (sum(y) - cumsum(y))/sum(y) #calculates the proportional cumulative sum of secondary extinctions
    x <- (object[, "no"]/max(object[, "no"])) #calculates the proportional primary extinctions
    ext.curve <- splinefun(x,y) #interpolates a function for the secondary extinctions
    ext.area <- integrate(ext.curve,0,1) # calculates the area below the curve
    return(as.numeric(ext.area[[1]]))
}