`specieslevel` <-
function(web, index="ALL", logbase="e", low.abun=NULL, high.abun=NULL) {
    # function to calculate bipartite web indices at the species level
    #
    # web    interaction matrix, with higher trophic level as columns
    # index  vector of indices to be calculated for each trophic level of the web;
    #        options are: "specs" for number of species, "species degree", "dependence",
    #        "d" for Blüthgen's d, "species strength" as sum of dependencies,
    #        "interaction" for interaction push/pull (our version of dependence
    #        asymmetry: see details), "PSI" for pollination service index (or pollinator
    #        support index, depending on the trophic level), "niche overlap" or "ALL"
    #        for all the aforementioned (default).
    # dist   distance metric to be used to calculate niche overlap; defaults to
    #        Horn's index, which is the recommendation of Krebs (Ecological Methodology);
    #        for other options see vegdist {vegan}
    # indices to be calculated:
    # "specs" for number of species,
    # "d" for Blüthgen's species-level diversity,
    # "species degree" for the sum of interactions per species,
    # "species strength" for Bascompte's summed dependence values per species (i.e. quantitative species degree),
    # "PSI" for pollination service index (only for plants in pollination webs).
    #
    # Carsten Dormann, Jochen Fründ & Denis Lippok, April/May 2007


    # m <- matrix(c(4,7,0,0,9,1,2,0,5), 3, byrow=TRUE)

    web <- empty(web) # delete unobserved species

    allindex <- c("species number", "degree", "ND", "dependence", "strength", "interaction", "PSI", "NS", "BC", "CC", "Fisher", "diversity", "effective partners", "d")

    if ("ALL" %in% index) index <- allindex
    if ("ALLBUTD" %in% index) index <- allindex[-c(1,4)]
#    out <- list("higher trophic level"=1, "lower trophic level"=1)
    out <- list()

   #----------------------------------------------------------------------------
   # species number
    if ("species number" %in% index) {
        spH <- c("number of species"=NCOL(web))
        spL <- c("number of species"=NROW(web))
        out$"higher trophic level"$"number of species" <- spH
        out$"lower trophic level"$"number of species" <- spL
    }

    #----------------------------------------------------------------------------
    # species degrees:
    if ("degree" %in% index){
      sdL <- rowSums(web>0)
      sdH <- colSums(web>0)
      out$"higher trophic level"$"species degree" <- sdH
      out$"lower trophic level"$"species degree" <- sdL
    }
    
    if ("ND" %in% index){
      nds <- ND(web)
      out$"higher trophic level"$"normalised degree" <- nds[[2]]
      out$"lower trophic level"$"normalised degree" <- nds[[1]]
    }

    #----------------------------------------------------------------------------
    # dependence values, following the lead by Bascompte et al. 2006 (Science) and
    # modifications suggested by Blüthgen et al. 2007 (Current Biology)

    if (any(c("strength", "dependence", "interaction") %in% index)){
      depL <- web/matrix(rowSums(web), nrow=NROW(web), ncol=NCOL(web), byrow=FALSE)
      depH <- web/matrix(colSums(web), nrow=NROW(web), ncol=NCOL(web), byrow=TRUE)
      if ("dependence" %in% index){
          out$"higher trophic level"$"dependence" <- depH
          out$"lower trophic level"$"dependence" <- depL
      }

      if ("strength" %in% index){
          # strength = sum of dependences for a species (referenced in Bascompte et al. 2006)
          SL <- rowSums(depH) # a plant's strength are the sum of the dependencies of all its pollinators
          SH <- colSums(depL) # accordingly ...
          out$"higher trophic level"$strength <- SH
          out$"lower trophic level"$strength <- SL
      } 

      #----------------------------------------------------------------------------
      # Interaction asymmetry (Vazquez et al. 2007, Oikos); rather similar to dependence above, really
      if ("interaction" %in% index) {
        Dij <- depH-depL  # positive values indicate a stronger effect of i (=plants) on j (bees) than vice versa
        Ailow <- rowSums(Dij)/rowSums(web>0)
        Aihigh <- colSums(-Dij)/colSums(web>0)
        out$"higher trophic level"$"interaction push/pull" <- Aihigh
        out$"lower trophic level"$"interaction push/pull" <- Ailow
      }
    }

    #----------------------------------------------------------------------------
    # Pollination webs only: pollination service index for each pollinator species
    if ("PSI" %in% index){
        PSI <- function(web, beta=1){
            # calculates the average contribution per visit for each pollinator species
            # (which in itself depends on the specialisation and abundance of the bees,
            # as well as the abundance of the plant species)
            #
            # web   a pollination web, with plants in rows
            # beta  a parameter accounting for the fact that two repeated landings are
            #       required to transfer pollen: one for the source, one for the sink; a
            #       value of 2 (default) assumes independent spacing of the plants; a value
            #       of 1 implies infinitely storing pollen from source to sink plant
            #
            # developed by Dormann, Blüthgen & Gruber, 3 May 2007
            #
            # example:
            # m <- matrix(c(4,4,0,4,1,7), nrow=3, byrow=TRUE)
            # PSI(m, beta=1)

            Wi. <- matrix(rep(colSums(web), NROW(web)), nrow=NROW(web), byrow=TRUE)
            W.j <- matrix(rep(rowSums(web), NCOL(web)), ncol=NCOL(web), byrow=FALSE)

            PSImat <- web/W.j * (web/Wi.)^beta
            (PSI <- colSums(PSImat))
        }
        psi <- PSI(web)
        out$"higher trophic level"$"Pollination Service Index PSI" <- psi
        out$"lower trophic level"$"Pollinator Support Index PSI" <- PSI(t(web))

    }
    
    #----------------------------------------------------------------------------
    if ("NS" %in% index){
#      require(sna) # which brings the function geodist used in nodespec
      NS <- nodespec(web)
      out$"higher trophic level"$"node specialisation index" <- NS$higher
      out$"lower trophic level"$"node specialisation index" <- NS$lower
    }
    if ("BC" %in% index){
      bcs <- BC(web)
      out$"higher trophic level"$"betweenness" <- bcs[[2]]
      out$"lower trophic level"$"betweenness" <- bcs[[1]]
    }
    if ("CC" %in% index){
      ccs <- CC(web)
      out$"higher trophic level"$"closeness" <- ccs[[2]]
      out$"lower trophic level"$"closeness" <- ccs[[1]]
    }


    #----------------------------------------------------------------------------
    if ("Fisher" %in% index){
      ff.low <- try(suppressWarnings(fisher.alpha(web, MARGIN=1)), silent=TRUE)
      ff.high <- try(suppressWarnings(fisher.alpha(web, MARGIN=2)), silent=TRUE)
      out$"lower trophic level"$"Fisher alpha" <- if (!inherits(ff.low, "try-error")) ff.low else rep(NA, times=NROW(web))
      out$"higher trophic level"$"Fisher alpha" <- if (!inherits(ff.high, "try-error")) ff.high else rep(NA, times=NCOL(web))

    }


    #----------------------------------------------------------------------------    
    if (any(c("diversity", "effective partners") %in% index)){
        preytot.mat <- matrix(rep(colSums(web), NROW(web)), NROW(web), byrow=TRUE)
        preyprop.mat <- web/preytot.mat  # = b_ik/b_.k in the first formula
        #H_Nk is the diversity index of inflow (diversity of flower visits for each pollinator)
        predtot.mat <- matrix(rep(rowSums(web), NCOL(web)), NROW(web), byrow=FALSE)
        predprop.mat <- web/predtot.mat  # = b_kj/b_.k in the second formula

        if (logbase==2 | logbase=="2"){
           H_Nk <- apply(preyprop.mat, 2, function(x) -sum(x*log2(x), na.rm=TRUE))
           #H_Pk is the diversity index of pollinators for each plant species
           H_Pk <- apply(predprop.mat, 1, function(x) -sum(x*log2(x), na.rm=TRUE))
           # next, we need the reciprocals of this
           # note that the ifelse is only needed if the web contains prey that is
           # not eaten or predators that don't eat ...
           n_Nk <- ifelse(colSums(web)!=0, 2^H_Nk, 0)
           n_Pk <- ifelse(rowSums(web)!=0, 2^H_Pk, 0)
        }
        if (logbase=="e"){ # same code as above, just with "e"
             H_Nk <- apply(preyprop.mat, 2, function(x) -sum(x*log(x), na.rm=TRUE))
             H_Pk <- apply(predprop.mat, 1, function(x) -sum(x*log(x), na.rm=TRUE))
             n_Nk <- ifelse(colSums(web)!= 0, exp(H_Nk), 0)
             n_Pk <- ifelse(rowSums(web)!= 0, exp(H_Pk), 0)
        }
        if ("diversity" %in% index){
              out$"higher trophic level"$"partner diversity" <- H_Nk
              out$"lower trophic level"$"partner diversity" <- H_Pk
        }
        
        if ("effective partners" %in% index){
              out$"higher trophic level"$"effective partners" <- n_Nk
              out$"lower trophic level"$"effective partners" <- n_Pk
        }
    }

   #----------------------------------------------------------------------------
   # species-level standardised diversity index d
   # di represents the lower level's partner diversity, correcting for their relative
   # abundances; a common pollinator, e.g., will thus have to be more overrepresented than
   # a rare pollinator to have the same contribution to the index
    if ("d" %in% index){
        dsL <- dfun(web, abuns=high.abun)[[1]]
        dsH <- dfun(t(web), abuns=low.abun)[[1]]
        out$"higher trophic level"$"d" <- dsH
        out$"lower trophic level"$d <- dsL
    }

    #---------------------------------------------------------------------------
    if (!("dependence" %in% index)) {
        out[[1]] <- as.data.frame(out[[1]])
        out[[2]] <- as.data.frame(out[[2]])
    }
    
    out
}

