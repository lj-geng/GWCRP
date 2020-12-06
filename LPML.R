# packages needed: "matrixStats"
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("matrixStats")
ipak(packages)

###########################################################################
## Calculate LPML given posterior samples after burnin
###########################################################################

LPML <- function(dat, out, p1=p1, cuts) 
{
  #===========================================================#
  ## Input:  dat = original survival dataset  ##
  ##         out = the result from GWCRP      ##
  ##         p1  = number of covariate coefficients ##
  ##         cuts = cutting points for pch    ##
  
  ## Output:
  ##         LPML = LPML ##
  #===========================================================#
  
  breaks = c(cuts, max(dat$time)+1)
  Xmat <- as.matrix(subset(dat, select = -c(time, status, site)))
  p2 <- length(cuts) # number of pieces
  LPML <- 0
  for (i in 1:nrow(dat)){
    x <- Xmat[i, ]
    time <- dat$time[i]
    status <- dat$status[i]
    loglik <- lapply(out$Iterates, function(y){
      cluster <- y$zout[dat$site[i]]
      beta <- y$betaout[cluster, ][1:p1]
      lambda <- y$betaout[cluster, ][-(1:p1)]
      exbeta <- exp(x %*% beta)
      ll <- status * lambda[as.numeric(
        cut(time, breaks, include.lowest = TRUE, right = FALSE))] +
        status * (x %*% beta)
      temp <- 0
      for (k in 1:p2){
        temp <- temp + exp(lambda[k]) * exbeta *
          max(0, min(time, breaks[k+1]) - breaks[k])
      }
      ll <- ll - temp
      return(ll)
    })
    LPML <- LPML - logSumExp( - unlist(loglik) - log(length(out$Iterates)))
  }
  return(LPML)
}