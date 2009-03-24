
`plotweb` <-
function(web, method = "cca", empty = TRUE, labsize = 1, ybig = 1,
    y_width = 0.1, spacing = 0.05, arrow="no", col.interaction="grey80",
    col.pred = "grey10", col.prey="grey10", lab.space=1,
    bor.col.interaction ="black", bor.col.pred="black", bor.col.prey="black",
    lablength = NULL, sequence=NULL,low.abun=NULL,low.abun.col="green",
    bor.low.abun.col ="black",
    high.abun=NULL, high.abun.col="red", bor.high.abun.col="black",
    text.rot=0, text.high.col="black", text.low.col="black")
{
  op <- par(no.readonly = TRUE)
  if (empty) web <- empty(web) else method <- "normal"
  web <- as.matrix(web) # to convert data.frames into matrix: needed for cumsum

  prey.order <- 1:dim(web)[1]
  pred.order <- 1:dim(web)[2]



  meths <- c("normal", "cca")
  meths.match <- pmatch(method, meths)
  if (is.na(meths.match)) stop("Choose plot-method: normal/cca.\n")
  if (meths.match==2)
  { #for the option "cca" the web is first re-arranged and then treated as "normal"; thus: there is no "else" to this "if".

  ## Problem: cca sometimes doesn't get the compartments right!!!!!
  ## Function "compart" returns a matrix with links assigned to compartments
  ## So, we need to extract the compartments there and put them in sequence:
    co <- compart(web)
    if (co$n.compart>1){ #do the arrangement for each compartment separately
      require(vegan)
      row.seq <- NULL
      col.seq <- NULL
      for (m in 1:co$n.compart){
        comp.member <- which(abs(co$cweb)==m, arr.ind=TRUE)
        rs <- unique(comp.member[,1])
        cs <- unique(comp.member[,2])
        if (length(rs) < 3 | length(cs) < 3){
          row.seq <- c(row.seq, rs)
          col.seq <- c(col.seq, cs)
        } else { #works fine for webs with only one compartment
          ca <- cca(web[rs, cs])
          row.seq <- c(row.seq, rs[order(summary(ca)$sites[,1], decreasing=TRUE)])
          col.seq <- c(col.seq, cs[order(summary(ca)$species[,2], decreasing=TRUE)])

        }

      }
      web <- web[row.seq, col.seq]
       prey.order <- row.seq
       pred.order <- col.seq
    } else {
      ca <- cca(web)
      web <- web[order(summary(ca)$sites[,1], decreasing=TRUE), order(summary(ca)$species[,1], decreasing=TRUE)]
      prey.order <- order(summary(ca)$sites[,1], decreasing=TRUE)
      pred.order <- order(summary(ca)$species[,1], decreasing=TRUE)
    }
  } # end for meths.match==2 condition; start of the "normal" plotting

   if (!is.null(sequence)) {
       cs <- which(sequence$seq.pred %in% colnames(web))
       rs <- which(sequence$seq.prey %in% rownames(web))

       for (i in 1:dim(web)[2]) pred.order[i] <- which(sequence$seq.pred[i]==colnames(web))
       for (i in 1:dim(web)[1]) prey.order[i] <- which(sequence$seq.prey[i]==rownames(web))
       web <- web[sequence$seq.prey[rs], sequence$seq.pred[cs]]
   }

   websum <- sum(web)
   difff <- diffh <-0 # if no abundances are set leave plotsize as before

        ###rearrange if lowfreq is set  # lowfreq is a named vector!!!
        if (!is.null(low.abun)) {

        lowfreq = rowSums(web)

        dummy <- lowfreq

        for (i in 1:length(low.abun) )
        {
        ind <- which(names(low.abun)[i] == names(dummy))

        lowfreq[ind] <- lowfreq[ind]+low.abun[i]

        }
        #websum <- sum(lowfreq)
        difff = (lowfreq-rowSums(web))/websum
        }

        ###rearrange if highfreq is set  # lowfreq is a named vector!!!
        if (!is.null(high.abun)) {

        highfreq = colSums(web)

        dummy <- highfreq

        for (i in 1:length(high.abun) )
        {
        ind <- which(names(high.abun)[i] == names(dummy))

        highfreq[ind] <- highfreq[ind]+high.abun[i]

        }
        #websum <- sum(highfreq)
        diffh = (highfreq-colSums(web))/websum
        }






        if (is.null(high.abun)) pred_prop <- colSums(web)/websum else pred_prop <- highfreq/websum
        if (is.null(low.abun)) prey_prop <- rowSums(web)/websum else prey_prop <- lowfreq/websum
        n.pred <- length(pred_prop)
        n.prey <- length(prey_prop)
        pred_x <- 0
        pred_xold <- -1
        pred_versatz <- 0
        pred_y <- 1.5
        prey_x <- 0
        prey_xold <- -1
        prey_versatz <- 0
        prey_y <- 0.5
        if (length(colnames(web)) == 0)
            colnames(web) <- colnames(web, do.NULL = FALSE)
        if (length(rownames(web)) == 0)
            rownames(web) <- rownames(web, do.NULL = FALSE)
        if (!is.null(lablength))
            colnames(web) <- substr(colnames(web), 1, lablength)
        if (!is.null(lablength))
            rownames(web) <- substr(rownames(web), 1, lablength)

        par(mai = c(0.2, 0.2, 0.2, 0.2))
        pred_spacing =  (n.prey - 1)/(n.pred - 1)
        prey_spacing =   (n.pred - 1)/(n.prey -1)
        pred_spacing <- pred_spacing*spacing
        prey_spacing <- prey_spacing*spacing

        if (n.pred>n.prey) prey_spacing <- pred_spacing*(n.pred-1)/(n.prey-1) else pred_spacing<- prey_spacing*(n.prey-1)/(n.pred-1)

        if (!is.null(low.abun)) pred_spacing <- pred_spacing+sum(difff)/n.pred

        if (!is.null(high.abun)) prey_spacing <- prey_spacing+sum(diffh)/n.prey


        wleft = 0
        wright = (max(n.pred, n.prey)) * min(prey_spacing,pred_spacing) +1+max(sum(diffh),sum(difff))

        wup <- 1.6 + y_width + lab.space * 0.05 #we need some space for labels
        wdown <- 0.4 - y_width - lab.space * 0.05 #we need some space for labels

        plot(0, type = "n", xlim = range(wleft, wright), ylim = range(wdown/ybig,
            wup * ybig), axes = FALSE, xlab = "", ylab = "")
        pred_x = 0
        hoffset <- 0
        links <- 0
        rechts <- 0
        hoehe <- strheight(colnames(web)[1], cex = 0.6)
        for (i in 1:n.pred) {
            rect(pred_x, pred_y - y_width, pred_x + pred_prop[i],
                pred_y + y_width, col = col.pred[(pred.order[i]-1) %% (length(col.pred))+1], border=bor.col.pred)
            #### coloured boxes at the end if highfreq is given
            if (!is.null(high.abun))
              {
              rect(pred_x + pred_prop[i]-diffh[i], pred_y - y_width, pred_x + pred_prop[i],
                pred_y + y_width, col = high.abun.col[(pred.order[i]-1) %% (length(high.abun.col))+1],                           border=bor.high.abun.col[(pred.order[i]-1) %% (length(bor.high.abun.col))+1])
              }

            breite <- strwidth(colnames(web)[i], cex = 0.6 *
                labsize)
            links <- pred_x + pred_prop[i]/2 - breite/2
            if (links < rechts && i > 1)
                hoffset = hoffset + hoehe
            else {
                rechts <- pred_x + pred_prop[i]/2 + breite/2
                hoffset <- 0
            }
            if (text.rot==90) {hoffset=0; ad =c(0,0.3)} else ad=c(0.5,0.4)
            text(pred_x + pred_prop[i]/2, pred_y + y_width +
                hoehe + hoffset, colnames(web)[i], cex = 0.6 *
                labsize, offset = 0, srt=text.rot, adj=ad,
                col=text.high.col[(pred.order[i]-1) %% (length(text.high.col))+1])
            pred_x <- pred_x + pred_prop[i] + pred_spacing
        }
        prey_x <- 0
        links <- 0
        rechts <- 0
        hoehe <- strheight(rownames(web)[1], cex = 0.6)
        hoffset <- hoehe
        for (i in 1:n.prey) {
            rect(prey_x, prey_y - y_width, prey_x + prey_prop[i],
                prey_y + y_width, col = col.prey[(prey.order[i]-1) %% (length(col.prey))+1], border=bor.col.prey[(prey.order[i]-1) %% (length(bor.col.prey))+1])
            #### coloured boxes at the end if lowfreq is given
            if (!is.null(low.abun))
              {
              rect(prey_x + prey_prop[i]-difff[i], prey_y - y_width, prey_x + prey_prop[i],
                prey_y + y_width, col = low.abun.col[(prey.order[i]-1) %% (length(low.abun.col))+1], border=bor.low.abun.col[(prey.order[i]-1) %% (length(bor.low.abun.col))+1])
              }
            breite <- strwidth(rownames(web)[i], cex = 0.6 *
                labsize)
            links <- prey_x + prey_prop[i]/2 - breite/2
            if (links < rechts && i > 1)
                hoffset = hoffset + hoehe
            else {
                rechts <- prey_x + prey_prop[i]/2 + breite/2
                hoffset <- hoehe
            }
            if (text.rot==90) {hoffset=hoehe; ad =c(1,0.3)} else ad=c(0.5,0.4)
            text(prey_x + prey_prop[i]/2, prey_y - y_width -
                hoffset, rownames(web)[i], cex = 0.6 * labsize,
                offset = 0, srt=text.rot, adj=ad,
                col=text.low.col[(prey.order[i]-1) %% (length(text.low.col))+1])
            prey_x <- prey_x + prey_prop[i] + prey_spacing

        }
        px <- c(0, 0, 0, 0)
        py <- c(0, 0, 0, 0)
        pred_x <- 0
        zwischenweb <- web
        XYcoords <- matrix(ncol = 2, nrow = length(zwischenweb))
        for (i in 1:length(zwischenweb)) {
            XYcoords[i,1:2 ] <- which(zwischenweb == max(zwischenweb),
                arr.ind = TRUE)[1, ]

            zwischenweb[XYcoords[i, 1], XYcoords[i, 2]] <- -1
        }
        y1 <- pred_y - y_width
        y2 <- y1
        y3 <- prey_y + y_width
        y4 <- y3
        for (p in 1:sum(web > 0)) {
            i <- XYcoords[p, 1]
            j <- XYcoords[p, 2]

            if (j == 1 & i == 1)
                x1 <- 0   else x1 <- (j - 1) * pred_spacing + cumsum(web)[(j -
                1) * nrow(web) + (i - 1)]/websum
            if (!is.null(high.abun) && j>1) x1 <- x1 +cumsum(diffh)[j-1]
            x2 <- x1 + web[i, j]/websum
            if (arrow=="up" || arrow=="both") {x2<-(x1+x2)/2; x1<-x2}
            if (arrow=="up.center"|| arrow=="both.center")
               {
               if (j!=1)  {x2 <-  (j - 1) * pred_spacing + cumsum(web)[(j -
                1) * nrow(web) ]/websum +colSums(web)[j]/websum/2
               if (!is.null(high.abun)) x2 <- x2 +cumsum(diffh)[j-1]
               x1<-x2
                }  else
                   {
                   x2=colSums(web)[j]/websum/2; x1<-x2
                   }
               }
                tweb <- t(web)
            if (j == 1 & i == 1)
                x3 <- 0   else x3 <- (i - 1) * prey_spacing + cumsum(tweb)[(i -
                1) * nrow(tweb) + (j - 1)]/websum
            if (!is.null(low.abun) && i>1) x3 <- x3 +cumsum(difff)[i-1]
            x4 <- x3 + tweb[j, i]/websum
            if (arrow=="down" || arrow=="both") {x4<-(x3+x4)/2; x3<-x4}
            if (arrow=="down.center" || arrow=="both.center")
               {
               if (i!=1)  {x3 <-  (i - 1) * prey_spacing + cumsum(tweb)[(i -
                1) * nrow(tweb) ]/websum +colSums(tweb)[i]/websum/2
                if (!is.null(low.abun)) x3<- x3 +cumsum(difff)[i-1]
                x4<-x3

                }  else
                   {
                   x3=colSums(tweb)[i]/websum/2; x4=x3;
                   }
               }
            # calculate color of interaction based on web order
            icol <- col.interaction[((prey.order[XYcoords[p,1]]-1)* (length(pred.order))+(pred.order[XYcoords[p,2]]-1)) %% (length(col.interaction))+1]
            bicol <- bor.col.interaction[((prey.order[XYcoords[p,1]]-1)* (length(pred.order))+(pred.order[XYcoords[p,2]]-1)) %% (length(bor.col.interaction))+1]
            polygon(c(x1, x2, x4, x3), c(y1, y2, y4, y3), col = icol, border=bicol) #bor.col.interaction[(p-1) %% (length(bor.col.interaction))+1])
        }
par(op)
}

