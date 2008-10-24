library(bipartite)
#------- helper functions: enter only once
extract <- function(x){
    #extracts degrees, d', interaction push/pull, strength and PSI for each species
    o <- x$"higher trophic level"
    list(degree=o$"species degree",
         d=o$"Blüthgen's d"$dprime,
         ipp=o$"interaction push/pull",
         strength=o$strength,
         PSI=o$"Pollination Service Index PSI" )
}
M4 <- function(x, ...) c(min=min(x, ...), mean=mean(x, ...), median=median(x, ...), max=max(x, ...))
#--------

data()
data(Safariland) # or other data set
plotweb(Safariland, method="cca") #method="normal"
visweb(Safariland, type="diagnoal", text="none")
?visweb #for details and options


set.seed(101) #fixes a random number for reproducibility
# run extinction sequence randomly:
(SEran <- slope.bipartite(second.extinct(Safariland, participant="higher", method="random", nrep=100))) # write down slope: 2.755149
# run extinction sequence by increasing abundance:
(SEabun <- slope.bipartite(second.extinct(Safariland, participant="higher", method="abun"))) # 6.419259

#calculate some indices at the species level:
s <- specieslevel(Safariland)
names(s)
str(s$"higher trophic level") # look at the structure of the output


# calculate the 4 Ms for each of the above measures:
(res <- sapply(extract(s), M4))
# re-organise table into a vector, and retain names:
result <- as.vector(res)
names(result) <- as.vector(outer(rownames(res), colnames(res), paste))

# calculate some network-level indices:
n <- unlist(networklevel(Safariland, index=c("web asymmetry", "cluster coefficient", "SA", "H2'", "ISA")))
# append to species-level results
result <- c("sec.ext.random"=SEran, "sec.ext.abun"=SEabun, result, n)
#write as text file (i.e. give name such as "safariland_indices.txt")
write.table(t(result), file=file.choose(), row.names=F, sep="\t") # can now be imported into spreadsheet

# repeat for each data set!

datasets <- c("Safariland", "barrett1987", "elberling1999", "kato1990", "memmott1999",
    "mosquin1967", "motten1982", "olesen2002aigrettes", "olesen2002flores", "schemske1978",
    "small1976", "vazarr", "vazcer", "vazllao", "vazmasc", "vazmasnc", "vazquec", "vazquenc")
results.m <- matrix(nrow=length(datasets), ncol)
for (i in datasets){
    idata <- data(i) # doesn't work, because data requires i as a real name!
    set.seed(101) #fixes a random number for reproducibility
    # run extinction sequence randomly:
    (SEran <- slope.bipartite(second.extinct(idata, participant="higher", method="random", nrep=100)))
    # run extinction sequence by increasing abundance:
    (SEabun <- slope.bipartite(second.extinct(idata, participant="higher", method="abun")))
    s <- specieslevel(Safariland)
    (res <- sapply(extract(s), M4))
    # re-organise table into a vector, and retain names:
    result <- as.vector(res)
    names(result) <- as.vector(outer(rownames(res), colnames(res), paste))
    # calculate some network-level indices:
    n <- unlist(networklevel(Safariland, index=c("web asymmetry", "cluster coefficient", "SA", "H2'", "ISA")))
    # append to species-level results
    result <- c("sec.ext.random"=SEran, "sec.ext.abun"=SEabun, result, n)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# analysis of all the above
# save spreadsheet as text-file, tabulator-delimited
extdata <- read.delim(file.choose())# if decimal "point" is a comma, use: dec=","

# loop through correlations for sec.ext.random.exponent:
(cor.ran <- sapply(3:24, function(x) cor(x=extdata[,1], y=extdata[,x], method="kendall")))
which.max(abs(cor.ran)) # best correlation
colnames(extdata)[which.max(cor.ran)+2]  #find out the name: +2 because first two columns were responses
cor(x=extdata[,1], y=extdata[,which.max(cor.ran)+2], method="ken") #check that value is correct

# repeat for sec.ext.abun.exponent:
(cor.abun <- sapply(3:24, function(x) cor(x=extdata[,2], y=extdata[,x], method="kendall")))
which.max(abs(cor.abun)) # best correlation
colnames(extdata)[which.max(cor.abun)+2]  #find out the name: +2 because first two columns were responses
cor(x=extdata[,2], y=extdata[,which.max(cor.abun)+2], method="ken") #check that value is correct


# plot findings:

par(mfrow=c(1,2)) #make two panels in plot region

# random extinction pattern:
plot(extdata[,1] ~ extdata[,"max strenght"], pch=16, cex=2, cex.lab=1.5,
    xlab="max strength", ylab="sec. extinction slope", main="random extinctions")
lines(lowess(extdata[,1] ~ extdata[,"max strenght"], f=1), lwd=3, col="grey")

# abundance extinction pattern:
plot(extdata[,2] ~ extdata[,"mean degree"], pch=16, cex=2, cex.lab=1.5,
    xlab="mean degree", ylab="sec. extinction slope", main="abundance-related extinctions")
lines(lowess(extdata[,2] ~ extdata[,"mean degree"], f=1), lwd=3, col="grey")



# look at correlation amongst specialisation measures:
library(Hmisc)
plot(varclus(extdata, similarity="pearson")) #  or: spearm,  hoeffding
round(cor(extdata, method="kendall"),2)
pairs(extdata[,2:7], pch="+") # doesn't fit on screen for all


#-------------------------------------------------------------------------------
null.test <- function(web, N=30, ...){
     # S L O W !! More repetitions will take longer ...
     #
     # A little null-model function to check, if the observed values actually are
     # much different to what one would expect under random numbers given the
     # observed row and column totals (i.e. information in the structure of the
     # web, not only in its species' abundances)
     #
     # web  a bipartite matrix with interactions
     # N    number of random webs to be created
     # ...  options passed on to other functions called internally:
     #             can go to t.test (for setting the conf.level to other than 0.95)
     #             or to index of networklevel or specieslevel (see there)
     #
     # returns a table with selected indices in rows and columns giving observed
     # values, mean null model values with CI-values, and t.test statistics (two-tailed)
     #
     # Does only work for indices directly returned by the function called (e.g. NOT for
     # degree distributions).
     
     obs <- networklevel(web, ...)
     null.list <- r2dtable(n=N, r=rowSums(web), c=colSums(web))
     res <- sapply(null.list, networklevel, ...)
     # t.test the difference between observed and expected from null model:
     t.test.mat <- matrix(nrow=nrow(res), ncol=6)
     rownames(t.test.mat) <- rownames(res)
     colnames(t.test.mat) <- c("obs", "null mean", "lower CI", "upper CI", "t", "P")
     for (i in 1:nrow(res)){
         t.res <- try(t.test(unlist(res[i,]), mu=obs[[i]], na.action=na.omit, ...))
         t.test.mat[i,] <- c(obs[[i]], mean(unlist(res[i,]), na.rm=T), t.res$conf.int,
                            t.res$statistic, t.res$p.value)
     }
     t.test.mat
}
#example:
null.test(web, index=c("number of species", "number of links", "linkage density", "web asymmetry",
          "number of compartments", "generality", "vulnerability", "interaction evenness",
          "compartment diversity", "cluster coefficient", "H2", "ISA", "SA",
          "extinction slope"), nrep=2, N=10)

