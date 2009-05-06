`second.extinct` <-
function(web, participant="higher", method="abun", nrep=10, details=FALSE){
    # returns a matrix quantifying secondary extinctions following the deletion of
    # a species from the web
    # web           an interaction web
    # participant   "high" or "low" or "both"
    # method        deletion following "random" or "abundance"
    # nrep          only for "random": number of repititions of random extinction sequence

    if (details==FALSE & pmatch(participant, c("both", "lower", "higher"))==1)
    {
          warning("\nFor random extinctions of both participants extinction sequences
            will differ in length. Simply averaging sequences can hence not be used. Thus,
            option 'details' will be set to FALSE internally.\n")
          details <- TRUE
    }

    one.second.extinct <- function(web=web, participant=participant, method=method)
    {
#        dead <- data.frame("no"=1:NCOL(web), "ext.plant"=rep(NA, NCOL(web)), "ext.poll"=rep(NA, NCOL(web)))
        dead <- matrix(nrow=0, ncol=3)
        colnames(dead) <- c("no", "ext.lower", "ext.higher")
        m2 <- web
        #for (i in 1:(NCOL(web)-1))
        i=1
        while (min(dim(m2))>1)
        {
            # extinct a species and count secondary extinctions:
            n <- extinction(m2, participant=participant, method=method)
            dead <- rbind(dead, c(i, attributes(m2 <- empty(n, count=TRUE))$empty))
            i <- i+1
        }
        dead2 <- rbind(dead, c(NROW(dead)+1, attributes(empty(m2, count=TRUE))$empty)) # counts extinction knock-on for the last species

      # If I use a random sequence, sometimes it takes only 23 steps to let all pollinators go extinct, sometimes 24. So, the matrix "dead" will have different dimensions. Correct this:
        if (nrow(dead)!=nrow(dead2)) { # Is dead of the right length?
          # ci <- which(sapply(1:3, function(x) all(dead[-nrow(dead),x]==1))==TRUE)  # a clumsy way to find out the participant subject to primary extinction: for this column, there is a 1 in every row;
          # BETTER: 
            ci <- pmatch(participant, c("both", "lower", "higher")) 
            dead2 <- matrix(0, nrow=ifelse(ci==2, nrow(web), ncol(web)), ncol=3)  # new dead-matrix
            colnames(dead2) <- colnames(dead) #use names from "dead"
            dead2[,1] <- 1:nrow(dead2) # step numbers
            dead2[1:nrow(dead),2:3] <- dead[1:nrow(dead),2:3] # plug matrix "dead" into "dead2", leaving the potentially longer matrix dead2 filled with zeros below nrow(dead).

        } else dead2=dead
       dead2

#        dead
    }

    if (is.vector(method)) sequence = method
    if (pmatch(method, c("abundance", "random"))==1)
    {
        out <- one.second.extinct(web=web, participant=participant, method=method)
    } else
    {
        o <- replicate(nrep, one.second.extinct(web=web, participant=participant, method=method), simplify=FALSE)

        if (details) out <- o else
        {
            z <- o[[1]]
            for (k in 2:length(o)) z <- z + o[[k]]
            out <- z/length(o)
        }
    }

    class(out) <- "bipartite"
    attr(out, "exterminated") <- c("both", "lower", "higher")[pmatch(participant, c("both", "lower", "higher"))]
    out


}

