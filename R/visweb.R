`visweb` <-
function (web, type="nested",  prednames=TRUE, preynames=TRUE, labsize=1, plotsize=12,
                   square="interaction", text="compartment", frame=NULL,
                   textsize=1, textcol="red", pred.lablength=NULL, prey.lablength=NULL,
                   clear=TRUE )
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
           web<- web[order(rowSums(web),decreasing=TRUE),order(colSums(web),decreasing=TRUE)]

        }

      n.pred <- length(web[1,])
      n.prey <- length(web[,1])
      plotsize = plotsize/2.54 # convert to cm
      mcol=max(web)
      #x bigger than y
      if (n.pred>n.prey) 
        {
        wx <- plotsize-1
        wy <- (plotsize-1)/n.pred*n.prey
        }  else
        {
        wy <- plotsize-1
        wx <- (plotsize-1)/n.prey*n.pred
        }
        
      
      
        
   #   windows(width=wx+2/2.54, height=wy+2/2.54)
      par(pin=c(wx,wy),mai=c(0.394,0.394,0.394,0.394) ,  omi=c(0,0,0,0) )
      plot(1,type="n",axes=FALSE, xlim=c(0,n.pred), ylim=c(0,n.prey), asp=1     )
      
      cellsize = wx/n.pred #same as wy/n.prey
      if (substr(text,1,1)=="i") s <- as.character(max(web)) else s="A"
      lettersize = strwidth(s,units="inches")
      clratio = cellsize/lettersize
   
        ## labels on x-axis (prednames)

        if (prednames==TRUE)
        {     
        for (ii in 1:n.pred)
          {
          s <- colnames(web)[ii]
          if (!is.null(s)) 
            {
            pop <- regexpr("_",s)
            if (pop>0) cn<-paste(c(substr(s,1,min(pop-1,9)),"\n", 
                substr(s,pop+1,pop+9)),sep="",collapse="") else 
                if (!is.null(pred.lablength))  cn <- substr(s,1,pred.lablength) else
                                          cn <- s
                
            mtext(cn, side=1,at=c(ii-0.5), cex=0.3*labsize*clratio )
            }
          }      
        }
        ## labels on y-axis (preynames)
        if (preynames==TRUE)
        {
        for (i in 1:n.prey)
          {
          s <- rownames(web)[n.prey-i+1]
          if (!is.null(s)) 
            {
            pop <- regexpr("_",s)
            if (pop>0) cr<-paste(c(substr(s,1,min(pop-1,9)),"\n", 
                substr(s,pop+1,pop+9)),sep="",collapse="") else 
                if (!is.null(prey.lablength))  cr <- substr(s,1,prey.lablength) else
                                          cr <- s

            mtext(cr, side=2, at=c(i-0.5), cex=0.3*labsize*clratio )
            } 
          }     
        }
     
            ### Crossfunktion für calculation of compartments ....
            cross = function(web,start,comp)    
            {   
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
        
            comp=1     #start with the first compartment
            w<-web
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
          if (substr(square,1,1)=="i") c= (1- (which(lev==web[n.prey-i+1,ii])-1)/nl) else
          if (substr(square,1,1)=="b") c =floor((1-web[n.prey-i+1,ii]/mcol)) else
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

