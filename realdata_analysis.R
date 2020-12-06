###########################################################################
## Louisiana Respiratory Cancer Data Analysis
###########################################################################

# packages needed: "MASS", "Matrix", "matrixStats", "mvtnorm", "eha", "parallel"
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("MASS", "Matrix", "matrixStats", "mvtnorm", "eha", "parallel")
ipak(packages)

source("simulation_data_generation.R")
source("GWCRP.R")
source("getDahl.R")
source("LPML.R")

# load graph distance matrix of pairs of counties in Louisiana
load("d_matrix.RData")

# load respiratory cancer Data in Louisiana
load("LA_RESPIR64.Rdata")
# select four regression coefficients
dat <- dat64[, c("SRV_TIME_MON", "DTH_CLASS", "AGE_DX",
                 "SEX", "GRADE", "HST_STGA", "site")]
dat$SEX <- as.numeric(dat$SEX) - 1
dat$GRADE <- as.numeric(dat$GRADE) - 1
dat$HST_STGA <- as.numeric(dat$HST_STGA) - 1
colnames(dat)[1:2] <- c("time", "status")
# save(dat, file = "respdata.RData")
# 
# load("respdata.RData")

n <- 64
# set four pieces for piecewise constant baseline hazards
timebreaks <- c(0, quantile(dat$time, c(1/6, 2/6, 1/2))+0.01) # c(0, 1.01, 4.01, 9.01)
p1 <- 4
p2 <- length(timebreaks)
p <- p1 + p2

betahat <- matrix(0, n, p)
hesshat <- array(0, c(p, p, n))
sigmahat <- array(0, c(p, p, n))

for (i in 1:64){

  mod <- phreg(Surv(time, status) ~ AGE_DX + SEX + GRADE + HST_STGA,
               dat[dat$site == i, ], dist = "pch", cuts = timebreaks[-1])
  mod$var
  betahat[i, ] <- c(mod$coefficients, mod$hazards)
  hesshat[, , i] <-
    computehessian(mydata=dat[dat$site == i, -7],
                   beta=mod$coefficients, lambda=mod$hazards,
                   quantile=timebreaks)
  sigmahat[, , i] <- solve(
    computehessian(mydata=dat[dat$site == i, -7],
                   beta=mod$coefficients, lambda=mod$hazards,
                   quantile=timebreaks), diag(1, p))
}
colnames(betahat) <- c(mod$covars, paste0("lambda", 1:p2))
# save(betahat, hesshat, sigmahat, file = "resp_estdat.Rdata")
# 
# load("resp_estdat.Rdata")

sigma0 <- 100 * diag(p)
alpha <- 1
initNClusters <- 10 # initial number of clusters for MCMC
niterations <- 5000 # number of iterations for MCMC
burnin <- 2000 # number of burnin iterations for MCMC

sim <- function(i){
  set.seed(h[i]*1000+1)
  out <- GWCRP(betahat, sigmahat, d_matrix, h[i], niterations, sigma0,
               alpha, initNClusters)
  Dahlout <- reorderDahl(getDahl(out, burn=burnin))
  LPML <- LPML(dat, list(Iterates = out$Iterates[(burnin+1):niterations]),
               p1 = p1, cuts = timebreaks)
  result <- list(h = h[i], out = out, Dahlout = Dahlout, LPML = LPML)
  return(result)
}

h <- c(seq(0, 10, by = 0.1))
cl <- makeCluster(Sys.getenv()["SLURM_NTASKS"], type = "MPI")
clusterExport(cl, varlist = ls())
result <- clusterApply(cl, 1:length(h), sim)
save(result, file = "resp_out.RData")