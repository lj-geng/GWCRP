###########################################################################
## Dahl's method to summarize the samples from the MCMC
###########################################################################

getDahl <- function(out, burn)
{
  #===========================================================#
  ## Input: out = the result from GWCRP ##
  ##        burn = the number of burn-in iterations ##
  
  ## Output:
  ##         zout = estimated clustering configuration, n * 1 vector##
  ##         betaout = estimated matrix of regression coefficients (dimension p1)
  ##                   and baselince hazards (dimension p2), k * (p1+p2) matrix ##
  #===========================================================#
  
  iters <- out$Iterates[-(1:burn)]
  n <- length(iters[[1]][[1]])
  niters <- length(iters)
  membershipMatrices <- lapply(iters, function(x){
    clusterAssign <- x$zout
    outer(clusterAssign, clusterAssign, FUN = "==")
  })
  membershipAverage <- Reduce("+", membershipMatrices)/niters
  SqError <- sapply(membershipMatrices, function(x, av) sum((x - av)^2),
                    av = membershipAverage)
  DahlIndex <- which.min(SqError)
  DahlAns <- iters[[DahlIndex]]
  attr(DahlAns, "iterIndex") <- burn + DahlIndex
  attr(DahlAns, "burnin") <- burn
  DahlAns
}

###########################################################################
## Reorder the zout and betaout from getDahl function, to make zout increasing
###########################################################################

reorderDahl <- function(x)
{
  #===========================================================#
  ## Input: x = the output of getDahl function ##
  
  ## Output:
  ##         zout = estimated clustering configuration, n * 1 vector##
  ##         betaout = estimated matrix of regression coefficients (dimension p1)
  ##                   and baselince hazards (dimension p2), k * (p1+p2) matrix ##
  #===========================================================#
  
  zout <- x$zout
  betaout <- x$betaout
  for (i in 1:length(unique(x$zout))){
    zout[x$zout == unique(x$zout)[i]] <- i
    betaout[i, ] <- x$betaout[unique(x$zout)[i], ]
  }
  Dahlout <- list(zout=zout, betaout=betaout)
  attr(Dahlout, "iterIndex") <- attr(x, "iterIndex")
  attr(Dahlout, "burnin") <- attr(x, "burnin")
  return(Dahlout)
}