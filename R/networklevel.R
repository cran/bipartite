`networklevel` <-
function(web, index="ALL", ISAmethod="Bluethgen", SAmethod="Bluethgen", extinctmethod="r",
    nrep=100, plot.it.extinction=FALSE, plot.it.dd=FALSE, CCfun=median, dist="horn",
    normalise=TRUE){
    ##
    ## web         interaction matrix, with lower trophic level in rows, higher in columns
    ##
    web <- empty(web)
    if (nrow(web) < 2 | ncol(web) <2) stop("Web is too small to calculate any reasonable index.")

    if (any(index %in% "ALL")) index <- c("number of species", "links per species",
          "connectance", "linkage density", "web asymmetry",
          "number of compartments", "generality", "vulnerability", "interaction evenness",
          "compartment diversity", "cluster coefficient", "H2", "ISA", "SA",
          "extinction slope", "degreedistribution", "niche overlap", "mean number of shared hosts", 
          "C-score", "togetherness", "V-ratio", "nestedness")
    out <- list()

    # set up enough panals for plotting:
    if (plot.it.extinction) {m=2; n=1; par(mfrow=c(m,n), mar=c(5,5,4,1))} else m=1
    if (plot.it.dd) {n=2; par(mfrow=c(m,n), mar=c(5,5,4,1))} else n=1


    #-------------------
    if ("number of species" %in% index) {
      out$"number of higher trophic species"=NCOL(web)
      out$"number of lower trophic species"=NROW(web)
    }

    #-------------------
    # mean number of links (= links/species)
    if (any(c("links per species", "interaction evenness", "linkage density", "vulnerability", "generality") %in% index)){

        L <- sum(web>0)/sum(dim(web))
        if ("links per species" %in% index) out$"links per species"=L
        
        
                #################
        # calculates linkage density for quantitative web
        # for formula see Tylianakis et al. (2006), supplement.
        # N refers to prey, P to predators
        #
        # web is a bipartite web with higher trophic level as columns and lower as rows
        preytot.mat <- matrix(rep(colSums(web), NROW(web)), NROW(web), byrow=TRUE)
        preyprop.mat <- web/preytot.mat  # = b_ik/b_.k in the first formula
        #H_Nk is the diversity index of inflow (diversity of flower visits for each pollinator)
        H_Nk <- apply(preyprop.mat, 2, function(x) -sum(x*log2(x), na.rm=TRUE))

        predtot.mat <- matrix(rep(rowSums(web), NCOL(web)), NROW(web), byrow=FALSE)
        predprop.mat <- web/predtot.mat  # = b_kj/b_.k in the second formula
        #H_Pk is the diversity index of pollinators for each plant species
        H_Pk <- apply(predprop.mat, 1, function(x) -sum(x*log2(x), na.rm=TRUE))

        # next, we need the reciprocals of this
        # note that the ifelse is only needed if the web contains prey that is
        # not eaten or predators that don't eat ...
        n_Nk <- ifelse(colSums(web)!=0, 2^H_Nk, 0)
        n_Pk <- ifelse(rowSums(web)!=0, 2^H_Pk, 0)

        # mean number of prey per predator
        G <- sum(colSums(web)/sum(web)*n_Nk)
        if ("generality" %in% index) out$"generality"=G
        # mean number of predators per prey
        V <- sum(rowSums(web)/sum(web)*n_Pk)
        if ("vulnerability" %in% index) out$"vulnerability"=V
        # linkage density
        LD_q <- 0.5*(V+G)
        if ("linkage density" %in% index) out$"linkage density"=LD_q
        #LD_qs <- LD_q/(NROW(web)+NCOL(web)) # according to Jason's text
        # interaction evenness
        p_i.mat <- web/sum(web)
        IE <- -sum(p_i.mat*log2(p_i.mat), na.rm=TRUE)/log2(sum(web!=0))

        evenness <- function(web){
            # calculates evenness of the numbers of individuals of different species in
            # a community, NOT according to formula in Müller et al. 1999, J. Anim. Ecol,
            # but according to the original formula in Alatalo 1981, Oikos
            # can be extended at some point to more indices ...
            pk <- web/sum(web)
            (Alatalo <- (1/sum(pk^2) -1) / (exp(-sum(pk*log(pk), na.rm=TRUE))-1 ))
        }

        E <- evenness(web)
        if ("interaction evenness" %in% index){
            out$"interaction evenness"=IE
            out$"Alatalo interaction evenness"=E
        }
    }
    if ("connectance" %in% index){
        # connectance: "the fraction of all possible links that are realized in a network", 
        # p. 12917 in Dunne et al. 2002
        out$connectance <- sum(web>0)/prod(dim(web))
    }

    if (any(c("number of compartments", "compartment diversity") %in% index)){
        CD <- function(co){
          if (co$n.compart>1){
            no <- NA
            for (i in 1:co$n.compart){
              comp <- which(co$cweb==i, arr.ind=TRUE) # who is in a compartment?
              no[i] <- length(unique(comp[,1])) + length(unique(comp[,2])) # how many species
            }
            no <- no/sum(dim(web)) # standardise for number of species in the web
            CD <- exp(-sum(no*log(no)))
          }  else {CD <- NA; warning("only one compartment")}
          CD
        }
       
        comps <- compart(web)
        if (class(comps)=="try-error") {
            ncompart <- compdiv <- NA
        } else  {
            ncompart <- comps$n.compart
            compdiv <- CD(comps)
        }
        out$"number of compartments"=ncompart
        out$"compartment diversity"=compdiv
    }
    
    #----------------------------------------------------------------------------
    if ("cluster coefficient" %in% index){
        cluster.coef <- function(web, full=FALSE, FUN=mean){
        # calculate cluster coefficient
        #
        # web   a bipartite web
        # full  logical; return cluster coefficients for each species in the network?
        # FUN   give a function to summarise individual cluster coefficients, defaults
        #       to 'mean'.
        #
        # The concept was introduced by Watts & Strogatz (1998) and is described under
        # in Wikipedia under http://en.wikipedia.org/w/index.php?title=Clustering_coefficient
        # Its main idea was to help identifying small-world networks, which should
        # have high cluster coefficients but low path lengths.
        #
        # Ci of species i is simply the number of realised links devided by the number of
        # possible links. At the network level, the CC for a network is the average
        # CC of its members. This is a little bit fishy, since the CCs are log-normal
        # distributed (in pollinations networks at least). Therefore with FUN also
        # other summary measures can be employed.
        # Because within 'bipartite' we look at 2-mode networks, mean C is:
        # C_lowertrophiclevel = C_highertrophiclevel = C_entirenetwork.
        #
        # Literature:
        # Watts DJ, Strogatz S (1998) Collective dynmics of 'small-world' networks. Nature 393:440-442

        web <- as.matrix(web)
        Ci.high <- colSums((web>0))/nrow(web)
        Ci.low <- rowSums((web>0))/ncol(web)
        CC <- FUN(Ci.high)
        if (full) out <- list("cluster coefficient"=CC, "CC values higher"=Ci.high,
            "CC values lower"=Ci.low) else out <- c("cluster coefficient"=CC)
        out
        }
        out$"cluster coefficient"=as.numeric(cluster.coef(web, FUN=CCfun, full=FALSE))
    }

    #---------------------------------------------------------------------------
    # Blüthgen's specialisation/non-conformity index
    if ("H2" %in% index) H2 <- as.numeric(H2fun(web)[1]) #1.element is the standardised H2 prime
        out$"H2"= ifelse(H2<0, 0, H2)
    #---------------------------------------------------------------------------
    # web asymmetry (Blüthgen et al. 2007, Fig. S2)
    if (any(c("SA", "ISA", "web asymmetry") %in% index)){

        if ("web asymmetry" %in% index) out$"web asymmetry" <- (NCOL(web)-NROW(web))/sum(dim(web))

        #----------------------------------------------------------------------------
        # Dependence asymmetry (Bascompte et al. 2006; Blüthgen et al. 2007, Fig. S2)
        depL <- web/matrix(rowSums(web), nrow=NROW(web), ncol=NCOL(web), byrow=FALSE)
        depH <- web/matrix(colSums(web), nrow=NROW(web), ncol=NCOL(web), byrow=TRUE)

        if (ISAmethod=="Bascompte" & "ISA" %in% index)  out$"dependence asymmetry"=mean(abs(depL-depH)/max(depL, depH))
        if (ISAmethod=="Bluethgen" & "ISA" %in% index) {
            web2 <- web
            # delete cells for species encountered only once:
            web2[, which(colSums(web)==1)] <- 0
            web2[which(rowSums(web)==1), ] <- 0
            rowsummat <- matrix(rowSums(web2), nrow=NROW(web2), ncol=NCOL(web2), byrow=FALSE)
            colsummat <- matrix(colSums(web2), nrow=NROW(web2), ncol=NCOL(web2), byrow=TRUE)
            depL <- web2/rowsummat
            depH <- web2/colsummat

            depL[depL<=0] <- NA
            depH[depH<=0] <- NA
            # now we need a correction to account for the fact that links with few (e.g. 2)
            # observations will have a minimum depL of 1/2: all on one species: depL=1,
            # one on each of two: depL=0.5
            depLprime <- (depL - 1/rowsummat)/(1 - 1/rowsummat) # assumes depLmin = 1/web2 and depLmax=1
            depHprime <- (depH - 1/colsummat)/(1 - 1/colsummat)
            out$"interaction strength asymmetry"=mean(as.matrix(depHprime-depLprime), na.rm=TRUE) #ranges from -1 to 1 /sum(depLprime, depHprime, na.rm=TRUE)
        }

        #----------------------------------------------------------------------------
        # Specialisation asymmetry (Blüthgen et al. 2007, Fig. S2)
        ## USE CODE FROM IndicesSpeciesLevel.r
        ## 2 options for calculating the "mean" SA:
        # either as Blüthgen et al: average weighted by number of interactions in the cell
        # or as mean of logarithms (since the dependencies follow a lognormal distribution)

        di <- dfun(web)$dprime
        dj <- dfun(t(web))$dprime
        if (SAmethod=="log"){
            lgmeani <- mean(log(di)); lgmeanj <- mean(log(dj))
            SA <- (lgmeanj-lgmeani)/sum(lgmeani, lgmeanj)  # ij-sequence changed because log changes sequence, too
        }
        if (SAmethod=="Bluethgen"){
            wmeani <- sum(di*rowSums(web))/sum(web)
            wmeanj <- sum(dj*colSums(web))/sum(web)
            SA <- (wmeanj-wmeani)/sum(wmeani, wmeanj)
        }
        out$"specialisation asymmetry"=SA

    }
    #----------------------------------------------------------------------------
    # species extinction curve:
    if ("extinction slope" %in% index){
        extL <- try(second.extinct(web=web, method=extinctmethod, nrep=nrep, participant="lower"))
        slopeL <- try(slope.bipartite(extL, col="green", pch=16, type="b", plot.it=plot.it.extinction))
        out$"extinction slope lower trophic level"=as.numeric(slopeL)
        
        extH <- try(second.extinct(web=web, method=extinctmethod, nrep=nrep, participant="higher"))
        slopeH <- try(slope.bipartite(extH, col="green", pch=16, type="b", plot.it=plot.it.extinction))
        out$"extinction slope higher trophic level"=as.numeric(slopeH)
    }

    #----------------------------------------------------------------------------
    # degree distribution fits:
    if ("degreedistribution" %in% index){
        dd <- try(degreedistr(web, plot.it=plot.it.dd, pure.call=FALSE))

        out$"degree distribution lower trophic level"=dd$"lower trophic level dd fits"
        out$"degree distribution higher trophic level"=dd$"higher trophic level dd fits"

    }

    #----------------------------------------------------------------------------
    # mean similarity of niches (niche overlap, sensu Krebs, Ecological Methodology)
    # vegdist requires "sites" to be in rows, therefore the web has to be transposed
    # to calculate dissimilarity between higher level species; similarity is simply
    # 1-dissimilarity:
    if ("niche overlap" %in% index) {
      NOhigher <- mean(1-vegdist(t(web), method=dist))
      NOlower <- mean(1-vegdist(web, method=dist))
      out$"higher trophic level niche overlap" <- NOhigher
      out$"lower trophic level niche overlap" <- NOlower
    }

    #-------------------
    # The Stone & Roberts's indices: S, T and C score:
    #-------------------
    if ("mean number of shared hosts" %in% index) {
      out$"mean number of shared hosts" <- 
            mean(designdist(web>0, method="J", terms="minimum"))
#      out$"mean number of shared mutualists lower trophic level" <- 
#            mean(designdist(t(web)>0, method="J", terms="minimum"))
    } 
    #-------------------
    if ("togetherness" %in% index){
      out$togetherness <- togetherness(web, normalise=normalise, na.rm=TRUE)
    }
    #-------------------
    if ("C-score" %in% index){
      out$"C-score" <- C.score(web, normalise=normalise, na.rm=TRUE)
    }

    #-------------------
    if ("V-ratio" %in% index){
      out$"V-ratio" <- V.ratio(web)
    }

    #-------------------
    if ("nestedness" %in% index){
      out$nestedness <- nestedness(web, null.models=FALSE)$temperature
    }

    return(out)
}
#networklevel(Safariland, index="H2")
#networklevel(Safariland, plot.it.dd=TRUE, plot.it.extinction=TRUE)