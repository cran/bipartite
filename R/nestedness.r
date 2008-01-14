
nestedness <-function(m, null.models=TRUE,
                      n.nulls=100, popsize=30, n.ind=7,n.gen=2000
                       )
{

# calculates matrix temperature using the binmatnest programm from Miguel Rodriguez-Girones
# Rodríguez-Gironès & Santamaría (2006). A new algorithm to calculate the nestedness
# temperature of presence-absence matrices. Journal of Biogeography 33:924-935.
#
# make sure matrix is a valid one as error proofing in the C++ function does not fully work
# and R may crash...
m<- ifelse(m>0,1,0)   # create a binary matrix
if (popsize<n.ind) n.ind <-popsize- 1 # you cannot pick more individuals then there are in the population...
bmn <- .C("bmn",
          mat=as.integer(t(m)),   #column dominated....
          n.rows = as.integer(ncol(m)),      # notice swapping of cols
          n.cols = as.integer(nrow(m)),      # and rows
          temperature = as.double(-1.0),
          n.nullmodels = as.integer(n.nulls),
          population.size = as.integer(popsize),
          n.individuals = as.integer(n.ind),
          n.generations = as.integer(n.gen),
          nullmodels = as.integer(null.models),
          p.null1 = as.double(-1.0),
          mean.temp.null1 = as.double(-1.0) ,
          var.temp.null1 = as.double(-1.0),
          p.null2 = as.double(-1.0),
          mean.temp.null2 = as.double(-1.0) ,
          var.temp.null2 = as.double(-1.0),
          p.null3 = as.double(-1.0),
          mean.temp.null3 = as.double(-1.0) ,
          var.temp.null3 = as.double(-1.0),
          pack.order.row = as.integer(rep(-1,nrow(m))),  # notice swapping of cols
          pack.order.col = as.integer(rep(-1,ncol(m))),  # and rows
          PACKAGE="bipartite")
# swap nrows and ncols again (due to different priorities of rows and cols)
s <-bmn$n.rows
bmn$n.rows <- bmn$n.cols
bmn$n.cols <- s

bmn$packed.matrix <- m[bmn$pack.order.row,bmn$pack.order.col]

bmn
}
