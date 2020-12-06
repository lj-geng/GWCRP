# packages needed: "mvtnorm","MASS"
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("mvtnorm","MASS")
ipak(packages)

###########################################################################
## GWCRP MCMC Collapsed sampler  
###########################################################################

GWCRP <- function(betahat, sigmahat, d_matrix, h, niterations, 
                  sigma0, alpha, initNClusters)
{
  ## Model: beta_hat|z,beta \sim MVN(beta_{z_i},sigmahat_{z_i}) ##
  ##        z|alpha \sim GWCRP(alpha, h)
  ##        beta_i \sim MVN(0, sigma0)
  
  #===========================================================#
  ## Input: 
  ##        betahat, n*p matrix ##
  ##        sigmahat, p*p*n arrary ##
  ##        d_matrix, n*n graph distance matrix ##
  ##        h = decay coefficient for graph distance ##
  ##        niterations = the total number of iterations ##
  ##        sigma0 = hyperparameters for the prior on beta in multivariate normal ##
  ##        alpha = the parameter in Dirichlet distribution that 
  ##                controls the relative size of clusters ##
  ##        initNClusters = the initial number of clusters ##
  
  ## Output: 
  ##         zout = clustering configuration, n by 1 vector ##
  ##         betaout = estimate beta, K by p matrix ##
  #===========================================================#
  
  n <- dim(betahat)[1]
  p <- dim(betahat)[2]
  mu0 <- rep(0, p)
  W <- exp(-d_matrix*h)
  diag(W) <- 0
  W[d_matrix==1] <- 1
  W <- W * n * (n-1) / sum(W)
  
  sigmahat_inv <- sigmahat
  sigmahat_invb <- matrix(0, n, p)
  for (s in 1:n){
    sigmahat_inv[, , s] <- ginv(sigmahat[, , s])
    sigmahat_invb[s, ] <- sigmahat_inv[, , s] %*% betahat[s, ]
  }
  sigma0_inv <- ginv(sigma0)

  #===== first initialize clustering and beta
  clusterAssign <- c(sample(1:initNClusters, size=initNClusters, replace=FALSE),
                     sample(1:initNClusters, size=n-initNClusters, replace=TRUE))

  beta <- rmvnorm(initNClusters, mean=mu0, sigma=sigma0) 
  # beta is a initNClusters * p matrix
  History <- vector("list", niterations)
  
  ## start Gibb's sampling
  for (iter in 1:niterations)
  {
    ## update z ##
    clusterSizes = table(as.factor(clusterAssign))
    nClusters = length(clusterSizes)
    
    for (i in 1:n)
    { # determine whether ith component is a singleton 
      if (clusterSizes[clusterAssign[i]] > 1){
        # if not a singleton, have nClusters + 1 choices
        clusterSizes[clusterAssign[i]] = clusterSizes[clusterAssign[i]] - 1
        # the probs for choosing exsiting clusters
        clusterProbs = sapply(1:nClusters, function(x) {
          sum(W[i, clusterAssign == clusterAssign[x]]) * 
            dmvnorm(betahat[i, ], beta[x, ], sigmahat[, , i])
        })
        # the prob for choosing a new cluster
        clusterProbs[nClusters+1] <- alpha * 
          dmvnorm(betahat[i, ], rep(0, p), sigmahat[, , i] + sigma0)
        
        # choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1, prob = clusterProbs)
        clusterAssign[i] <- cluster.i
        
        if (cluster.i > nClusters) {
          beta = rbind(beta, rmvnorm(1, mean = mu0, sigma=sigma0))
        }

       } else {
        # if singleton, have nClusters choices
        clusterAssign[clusterAssign > clusterAssign[i]] <- 
          clusterAssign[clusterAssign > clusterAssign[i]] - 1
        if (dim(beta)[1] > 2){
          beta <- beta[-clusterAssign[i], ]
        } else {
          beta <- matrix(beta[-clusterAssign[i], ], nrow = nClusters-1)
        }
        # the probs for choosing exsiting clusters
        clusterProbs = sapply(1:(nClusters-1), function(x) {
          sum(W[i, clusterAssign == clusterAssign[x]]) * 
            dmvnorm(betahat[i, ], beta[x, ], sigmahat[, , i])
        })
        # the prob for choosing a new cluster
        clusterProbs[nClusters] <- alpha * 
          dmvnorm(betahat[i, ], rep(0, p), sigmahat[, , i] + sigma0)
        
        # choose the cluster number for ith observation
        cluster.i <- sample(1:nClusters, size = 1, prob = clusterProbs)
        clusterAssign[i] <- cluster.i
        
        if (cluster.i == nClusters) {
          beta = rbind(beta, rmvnorm(1, mean = mu0, sigma=sigma0))
        }
       }
      
      clusterSizes <- table(as.factor(clusterAssign))
      nClusters <- length(clusterSizes)
    }
    
    ## update beta ##
    for (r in 1:nClusters){
      sigmahat_u <- ginv(
        apply(sigmahat_inv[, , clusterAssign == r], 1:2, sum) + sigma0_inv)
      mu0_u <- sigmahat_u %*% apply(matrix(sigmahat_invb[clusterAssign == r, ], 
                                           nrow = clusterSizes[r], 
                                           ncol=p), 2, sum)
      beta[r,] = rmvnorm(1, mean = mu0_u, sigma=sigmahat_u)
    }

    History[[iter]] <- list(zout = clusterAssign, betaout = beta)
     cat(" iteration:", iter,"\n",clusterAssign,"\n")
    
  }# for loop over iterations
  
  list(Iterates = History)
}