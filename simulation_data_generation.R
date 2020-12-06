# packages needed: "MASS", "Matrix", "survival"
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("MASS", "Matrix", "survival")
ipak(packages)

# options(digits=10) 

###########################################################################
## Generate survival data with piecewise constant baseline hazard
## (no time-dependent covariates)
###########################################################################

gendatapch <- function(n, beta, lambda, timebreaks, Crate, tmax)
{
  #===========================================================#
  ## Input: n = number of subjects
  ##        beta = the vector of regreesion coefficients ##
  ##        lambda = the vector of baseline hazards ##
  ##        Crate = the parameter for generate censoring time ##
  ##        timebreaks  = timebreaks of pch
  ##        tmax = the upper limit for observe/censoring time ##
  
  ## Output:
  ##        a data frame containing:
  ##         survival time, status and covariates ##
  #===========================================================#
  
  p1 <- length(beta) # number of covariates
  p2 <- length(lambda) # number of basline hazard pieces 

  # covariates generated from N(0, 1)
  x <- rnorm(n*p1, 0, 1) 
  x <- matrix(x, n, p1)
  colnames(x) <- paste0("x", 1:p1)
  
  t <- timebreaks
  dt <- diff(t)
  LD <- matrix(0, nrow=p2, ncol=p2-1) # create a help matrix
  LD[lower.tri(LD)] <- 1
  H <- -log(1-runif(n))/exp(x %*% beta)
  Ht <- as.vector(LD %*% (lambda[-p2]*dt))
  
  # Generate survival time
  time <- t[p2]+(H-Ht[p2])/lambda[p2] 
  for (i in 1:(p2-1)){
    time <- ifelse(Ht[i] < H & H <= Ht[i+1], 
                   t[i]+(H-Ht[i])/lambda[i], time)
  }
  # Genrate censoring time
  Ctime <- rexp(n, Crate)
  Ctime[Ctime>tmax] <- tmax
  status <- (time < Ctime) + 0
  time <- status * time + (1-status) * Ctime
  data.frame(time=time, status=status, x)
}

###########################################################################
## Compute the hessian matrix of (beta, log(lambda))
###########################################################################

computehessian <- function(mydata, beta, lambda, quantile)
{
  #===========================================================#
  ## Input: mydata = survival dataset ##
  ##        beta = regression coefficients, p1 * 1 vector ##
  ##        lambda = baseline hazards (positive), p2 * 1 vector ##
  ##        (Note that beta and lambda are the MLE)
  ##        quantile = timebreaks of pch, p2 * 1 vector with the first element 0 ##
  
  ## Output:
  ##         hessian matrix of (beta, log(lambda)), p * p matrix ##
  #===========================================================#

  Xmat <- as.matrix(mydata[, -(1:2)])
  p1 <- length(beta)
  p2 <- length(lambda)
  n <- nrow(Xmat)
  hessian <- matrix(0, (p1+p2), (p1+p2))
  exbeta <- exp(Xmat %*% beta) # exbeta is a n*1 matrix
  statusbench <- mydata$status[mydata$status==1] 
  timebench <- mydata$time[mydata$status==1]
  quantile <- c(quantile, max(mydata$time)+1)
  
  dvec <- tapply(statusbench, 
                 cut(timebench, quantile, right=FALSE, include.lowest=TRUE), 
                 sum)
  dvec[is.na(dvec)] <- 0
  hessian[(1+p1):(p1+p2), (1+p1):(p1+p2)] <- diag(dvec)
  
  for(r in 1:p1){
    for(s in 1:p1){
      temp1 <- 0
      for(j in 1:p2){
        temp2 <- 0
        for(i in 1:n){
          temp2 <- temp2 + Xmat[i, r] * Xmat[i, s] * exbeta[i, 1]* 
            max((min(mydata[i,1], quantile[j+1])-quantile[j]), 0) 
        }
        temp1 <- temp1 + temp2 * lambda[j]
      }
      hessian[r, s] <- temp1  
    }
  }
  
  for(r in 1:p1){
    for(s in (p1+1): (p1+p2)){
      temp1 <- 0
      j <- s-p1
      for(i in 1:n){
        temp1 <- temp1 + Xmat[i,r] * exbeta[i,1] * 
          max((min(mydata[i,1], quantile[j+1])-quantile[j]), 0) 
      }
      hessian[r, s] <- hessian[s, r] <- temp1 * lambda[j] 
    }
  }
  
  return(hessian)
}