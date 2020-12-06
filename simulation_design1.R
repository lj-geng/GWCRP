## Simulation Design 1, true cluster: 1, 2, 3

# packages needed: "MASS", "Matrix", "survival", "matrixStats"
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("MASS", "Matrix", "survival", "matrixStats", "eha", "parallel")
ipak(packages)


source("simulation_data_generation.R")
source("GWCRP.R")
source("getDahl.R")
source("LPML.R")

load("d_matrix.RData")

# true clustering configuration of Deigsn 1
a <- rep(0, 64)
a[c(3, 4, 26, 29, 32, 36, 38, 44, 45, 46, 47, 48, 52, 53, 55, 59)] <- 2
a[c(7, 8, 9, 11, 13, 14, 15, 16, 18, 21, 22, 25, 30, 31, 33, 34, 
    35, 37, 41, 42, 43, 54, 56, 60, 62, 64)] <- 3
a[c(1, 2, 5, 6, 10, 12, 17, 19, 20, 23, 24, 27, 28, 39, 40, 49, 
    50, 51, 57, 58, 61, 63)] <- 1

p1 <- 3 # number of regression coefficients / length of beta
p2 <- 3 # number of piecewise baseline hazards / langth of lambda
p <- p1 + p2
n <- 64 # number of counties
timebreaks <- c(0, 1.5, 6) # three pieces of pch

# true beta's and lambda's of three clusters
truebeta1 <- c(1, 0.5, 1)
truebeta2 <- c(1.5, 1, 1)
truebeta3 <- c(2, 0.5, 1.5)
truebetas <- list(truebeta1, truebeta2, truebeta3)
lambda1 <- 0.045 * c(1, 0.8, 1)
lambda2 <- 0.045 * c(1, 0.8, 0.8)
lambda3 <- 0.045 * c(0.8, 1, 1.1)
lambdas <- list(lambda1, lambda2, lambda3)

nsims <- 100 # simulation replicates

## ======= Generate simulation data ======= ##
simdata <- list()
simdata_betahat <- list()
simdata_hesshat <- list()
for (sim in 1:nsims){
  cat(sim)
  set.seed <- sim+123 # myseed same for different h, varies with sim replicates
  dat <- data.frame()
  betahat <- matrix(0, n, p)
  hesshat <- array(0, c(p, p, n))
  nsite <- 60 # subjects at each site
  for (i in 1:64){
    tempdat <- gendatapch(truebetas[[a[i]]], n = nsite,
                          lambda = lambdas[[a[i]]],
                          timebreaks = timebreaks, Crate = 0.01, tmax = 150)
    mod <- phreg(Surv(time, status) ~ x1+x2+x3, data = tempdat,
                 dist = "pch", cuts = timebreaks[-1])
    betahat[i, ] <- c(mod$coefficients, log(mod$hazards))
    hesshat[, , i] <-
      computehessian(mydata=tempdat, beta=mod$coefficients, lambda=mod$hazards,
                     quantile=timebreaks)
    tempdat$site <- i
    dat <- rbind(dat, tempdat)
  }
  simdata[[sim]] <- dat
  simdata_betahat[[sim]] <- betahat
  simdata_hesshat[[sim]] <- hesshat
}

simdata_sigmahat <- simdata_hesshat
for (i in 1:length(simdata_hesshat)){
  for (j in 1:64){
    simdata_sigmahat[[i]][,,j] <- ginv(simdata_hesshat[[i]][,,j])
  }
}

# save(simdata, file = "simdata_weak123_11142019.Rdata")
# save(simdata_betahat, simdata_hesshat, simdata_sigmahat,
#      file = "simestdata_weak123_11142019.Rdata")
save(simdata, file = "simdata_design1.Rdata")
save(simdata_betahat, simdata_hesshat, simdata_sigmahat, 
     file = "simestdata_design1.Rdata")
## ======= generate data finished ======= ##

load("simdata_design1.Rdata")
load("simestdata_design1.Rdata")

sigma0 <- 100 * diag(p)
alpha <- 1
initNClusters <- 10 # initial number of clusters for MCMC
niterations <- 2000 # iterations for MCMC
burnin <- 500 # burnin for MCMC

sim <- function(isim){
  # library(MASS)
  # library(mvtnorm)
  dat <- simdata[[isim]]
  betahat <- simdata_betahat[[isim]]
  sigmahat <- simdata_sigmahat[[isim]]

  simresult <- list()
  for (i in 1:length(h)){
    set.seed(h[i]*10000 + isim)
    out <- GWCRP(betahat, sigmahat, d_matrix, h[i], niterations, sigma0,
                 alpha, initNClusters)
    Dahlout <- reorderDahl(getDahl(out, burn=burnin))
    LPML <- LPML(dat, list(Iterates = out$Iterates[(burnin+1):niterations]),
                 p1 = p1, cuts = timebreaks)
    # simresult[[i]] <- list(h = h[i], out = out, Dahlout = Dahlout)
    simresult[[i]] <- list(h = h[i], out = out, Dahlout = Dahlout, LPML = LPML)
  }
  return(simresult)
}

h <- c(seq(0, 2, by = 0.2), seq(3, 10, by = 1))
cl <- makeCluster(Sys.getenv()["SLURM_NTASKS"], type = "MPI")
clusterExport(cl, varlist = ls())
simresult <- clusterApply(cl, 1:nsims, sim)

# save.image("weak123_11142019.RData")
save.image("simresult_design1.RData")




