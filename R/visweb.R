                        `visweb` <-
function (web, type="nested",  prednames=TRUE, preynames=TRUE, labsize=1, plotsize=12,
                   square="interaction", text="compartment", frame=NULL,
                   textsize=1, textcol="red", pred.lablength=NULL, prey.lablength=NULL,
                   clear=TRUE)
### image funktion für foodwebs
### Version 1.0d
### written by Bernd Gruber, 2006-05-07
###



{

     if (!is.matrix(web)) {web <- as.matrix(web); warning("Object converted to matrix.")}
     ### if frame was not defined, preset is to be sorted by diagonal
     if (type!="diagonal" && is.null(frame)) frame <- FALSE

     ## delete species without interactions
     if (clear==TRUE) web<-empty(web)

     ## order web by least numbers of crossings (along the diagonal) by a cca
     if (type=="diagonal")
        {
           require(vegan)
           web<-empty(web)
           ca <- cca(web)
           web <- web[order(summary(ca)$sites[,1], decreasing=TRUE), order(summary(ca)$species[,1], decreasing=TRUE)]
        }
     if (type=="nested")
        {
           web<-empty(web)
           web<- web[order(rowSums(web), decreasing=TRUE), order(colSums(web), decreasing=TRUE)]

        }

      n.pred <- ncol(web)
      n.prey <- nrow(web)
      plotsize = plotsize/2.54 # convert to cm
      mcol=max(web)
      #x bigger than y
      if (n.pred>n.prey)
        {
        wx <- plotsize
        wy <- (plotsize)/n.pred*n.prey
        }  else
        {
        wy <- plotsize
        wx <- (plotsize)/n.prey*n.pred
        }

   #calculate max length of prednames in inch
   m.predsize =max(strwidth(colnames(web),units="inches"))
   #calculate max length of preynames in inch
   m.preysize =max(strwidth(rownames(web),units="inches"))


   # calculate text size:
      cellsize = wx/n.pred #same as wy/n.prey
      if (substr(text,1,1)=="i") s <- as.character(max(web)) else s="A"
      lettersize = strwidth(s, units="inches")
      clratio = cellsize/lettersize
   # set up plotting region:
       mm<-max(m.predsize,m.preysize)
      par(pin=c(wx,wy),  omi=c(0,0,0,0),  mai=c(mm,mm,0.0,0.0) )

      plot(1, type="n", axes=FALSE, xlim=c(0,n.pred), ylim=c(0,n.prey),asp=1,xlab="",ylab="" )


        ## labels on x-axis (prednames)
        pnl<-0
    if (prednames && !is.null(colnames(web))) {
        for (ii in 1:n.pred) {
            s <- colnames(web)[ii]
            if (!is.null(s)) {
                pop <- regexpr("_", s)
                if (pop > 0) {
                  cn <- paste(c(substr(s, 1, min(pop - 1, 9)),
                    "\n", substr(s, pop + 1, pop + 9)), sep = "",
                    collapse = "")
                }
                else {
                  if (!is.null(pred.lablength)) {
                    cn <- substr(s, 1, pred.lablength)
                  }
                  else {
                    cn <- s
                  }
                }
                pnl[ii] <- cn
            }
        }

    axis(1, (1:n.pred) - 0.5, labels = pnl, tick = FALSE, mgp = c(0,
        0, 0), las = 2, cex.axis = 0.4 * labsize * clratio)
    }
    ynl <- 0
    if (preynames  && !is.null(rownames(web))) {
        for (i in 1:n.prey) {
            s <- rownames(web)[n.prey - i + 1]
            if (!is.null(s)) {
                pop <- regexpr("_", s)
                if (pop > 0) {
                  cr <- paste(c(substr(s, 1, min(pop - 1, 9)),
                    "\n", substr(s, pop + 1, pop + 9)), sep = "",
                    collapse = "")
                }
                else {
                  if (!is.null(prey.lablength)) {
                    cr <- substr(s, 1, prey.lablength)
                  }
                  else {
                    cr <- s
                  }
                }
                ynl[i] <- cr
            }
        }

    axis(2, (1:n.prey) - 0.5, labels = ynl, tick = FALSE, mgp = c(0,
        0, 0), las = 2, cex.axis = 0.4 * labsize * clratio)
    }
            ### Crossfunktion für calculation of compartments ....
            cross = function(web,start,comp) {
                n.r=nrow(web)
                n.c=ncol(web)
                r=start[1]
                c=start[2]
                web[r,c]=-comp     #assign a negative compartment number to interaction

                for (i in 1:n.r) #schaue senkrecht
                  {
                  if (web[i,c]>0)  web<-cross(web,c(i,c),comp)
                  }

                for (i in 1:n.c)    #schaue waagrecht
                  {
                   if (web[r,i]>0) web<-cross(web,c(r,i),comp)
                  }
                return (web)
            }

            comp = 1     #start with the first compartment
            w <- web
            ### assign compartment for each interaction
            while (max(w)>0) #any interactions left?
                {
                    start=which(w==max(w),arr.ind=TRUE)[1,]  #start at the highest number of interactions (arbitrary)
                    w<- cross(w,start,comp) #start recursion until no more neighbours in one compartment are found
                    comp = comp+1   #go to the next compartment
                }
            w <- abs(w) #make compartments positive numbers

        #draw coloured squares and labels
        nl <- length(unique(rank(web)))-1   # number of levels
        lev <- as.numeric(names(table(web)))


        for (i in 1:n.prey)
        {
        for (ii in 1:n.pred)
          {
          if (substr(square,1,1)=="c") c= 1-(w[n.prey-i+1,ii])*(1/(max(w)))  else
          if (substr(square,1,1)=="i") c= (1 - (which(lev==web[n.prey-i+1,ii])-1)/nl) else
          if (substr(square,1,1)=="b") c = floor((1-web[n.prey-i+1,ii]/mcol)) else
              c=1

          c=gray(c)
          rect(ii-1,i-1,ii,i, col=c)
          tc <- ""
          if (substr(text,1,1)=="i")
              {
              if (web[n.prey-i+1,ii]>0) tc <- web[n.prey-i+1,ii]
              } else if (substr(text,1,1)=="c")
              {
              if (w[n.prey-i+1,ii]>0) tc <- LETTERS[w[n.prey-i+1,ii]]
              }
          text(ii-0.5,i-0.5,tc, col =textcol, cex=1*textsize*clratio*0.5, adj=c(0.5,0.5))
          }
        }


        if (is.null(frame)) frame=TRUE
        if (frame)
        {
        # one square for each compartment
        for (i in 1:max(abs(w)))
          {
          squares <- which(w==i,arr.ind=TRUE)
          minsqx <- min(squares[,2])-1
          maxsqx <- max(squares[,2])
          minsqy <- min(squares[,1])-1
          maxsqy <- max(squares[,1])
          rect(minsqx,n.prey-maxsqy,maxsqx,n.prey-minsqy, lwd=3)
          }
        }
}

