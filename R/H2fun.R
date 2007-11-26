`H2fun` <-
function(web){
    # function to calculate a measure of specialisation in a bipartite web: H2'
    #
    # web   a bipartite interaction web, with e.g. pollinators as columns and plants as rows
    #
    # returns the normalised H2' and its subcomponents H2max, H2min and H2uncorrected
    #
    # by Jochen Fründ, Uni Würzburg (some streamlining by Carsten Dormann, UFZ Leipzig); April 2007
    # This function is based on the paper by Blüthgen, Menzel & Blüthgen 2006 (BMC Ecology).
    # for web implementation see also: http://itb1.biologie.hu-berlin.de/~nils/stat/
    #
    # Example:
    # set.seed(8793)
    # web <- matrix(rpois(81, 0.5), 9)
    # H2prime(web)

    tot <- sum(web)       #sämtliche Interaktionen im Netz
    rs <- rowSums(web)    #Interaktion der Pflanze mit sämtlichen Bestäubern
    cs <- colSums(web)    #sämtliche Interaktionen des jeweiligen Bestäbers

   #--------------- H2 uncorrected------------
    H2uncorr = -sum(web/tot*log(web/tot), na.rm=TRUE)

   #--------------- H2 max ------------
    # Key idea here is to allocate interactions one-by-one into the places where
    # they fit best according to the non-integer optimal web. Easy! (Not.)
    exexpec <- outer(rs, cs/tot) # non-integer optimal web
    expec <- matrix(0, nrow(web), ncol(web)) # empty web
    difexp <- exexpec-expec # where are differences between non-integer and 0-web greatest?
    newweb <- expec  # start new web
    webfull <- matrix("no", nrow(web), ncol(web)) # makes boolean web, set to 0
    while (sum(newweb)<tot) {
       webfull[which(rowSums(newweb)==rs),]="yo" # sets columns/rows with correct cs/rs to 1
       webfull[,which(colSums(newweb)==cs)]="yo"
       OK <- webfull=="no" # matrix of potential cells
       smallestpossible <- newweb==min(newweb[OK])  # find cell with lowest number of interactions (e.g. 0)
       greatestdif <- max(difexp[smallestpossible & OK]) # find cell value with largest different between "is" and "should"
       bestone <- which(OK & smallestpossible & difexp==greatestdif ) # find cell for all three conditions
       if (length(bestone)>1) bestone <- sample(bestone,1) # select randomly a cell, if different are possible
       newweb[bestone] <- newweb[bestone]+1 # put an interaction into that cell
    }
    H2_max <- -sum(newweb/tot*log(newweb/tot), na.rm=TRUE)

   #--------------- H2 min ------------
    # The key idea here is that allocating the col- & row-sums into the
    # web-matrix will yield the minimum H2 for the matrix (since it is the
    # maximum of aggregation possible.
    newweb <- matrix(0,length(rs),length(cs))
    rsrest=rs; csrest=cs
    while (sum(rsrest)!=0) {
        newweb[which(rsrest==max(rsrest))[1],which(csrest==max(csrest))[1]]=min(c(max(rsrest),max(csrest)))
        rsrest=rs-rowSums(newweb)
        csrest=cs-colSums(newweb)
    }
    # Now we have our new web, with maximum aggregation of observations,
    # and hence minimum H2
    Pnew <- newweb/sum(newweb)
    H2_min <- -sum(Pnew*log(Pnew), na.rm=TRUE)

   #--------------- H2'  ------------
    #ranging (between 0 and 1):
    H_2prime <- (H2_max - H2uncorr) / (H2_max - H2_min)

    #output:
    c("H2"=H_2prime, "H2min"=round(H2_min,3), "H2max"=round(H2_max,3), "H2uncorr"=round(H2uncorr,3))
    
}

