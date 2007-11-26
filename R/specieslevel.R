`specieslevel` <-
function(web, index="ALL") {
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

    if (index=="ALL") index <- c("specs", "species degree", "dependence", "d", "species strength",
              "interaction", "PSI", "niche overlap")
#    out <- list("higher trophic level"=1, "lower trophic level"=1)
    out <- list()

   #----------------------------------------------------------------------------
   # species number
    if ("specs" %in% index) {
        spH <- c("number of species"=NCOL(web))
        spL <- c("number of species"=NROW(web))
        out$"higher trophic level"$"number of species" <- spH
        out$"lower trophic level"$"number of species" <- spL
    }

   #----------------------------------------------------------------------------
   # species-level standardised diversity index d
   # di represents the lower level's partner diversity, correcting for their relative
   # abundances; a common pollinator, e.g., will thus have to be more overrepresented than
   # a rare pollinator to have the same contribution to the index
    if ("d" %in% index){
        dsL <- dfun(web)[[1]]
        dsH <- dfun(t(web))[[1]]
        out$"higher trophic level"$d <- dsH
        out$"lower trophic level"$d <- dsL
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
    # species degrees:
    if ("species degree" %in% index){
      sdL <- rowSums(web>0)
      sdH <- colSums(web>0)
      out$"higher trophic level"$"species degree" <- sdH
      out$"lower trophic level"$"species degree" <- sdL
    }

    #----------------------------------------------------------------------------
    # dependence values, following the lead by Bascompte et al. 2006 (Science) and
    # modifications suggested by Blüthgen et al. 2007 (Current Biology)

    if ("dependence" %in% index){
      depL <- web/matrix(rowSums(web), nrow=NROW(web), ncol=NCOL(web), byrow=FALSE)
      depH <- web/matrix(colSums(web), nrow=NROW(web), ncol=NCOL(web), byrow=TRUE)
      out$"higher trophic level"$"dependence" <- depH
      out$"lower trophic level"$"dependence" <- depL

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
    out
}

