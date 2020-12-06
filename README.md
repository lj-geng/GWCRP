# GWCRP

`GWCRP.R` is the MCMC Collapsed sampler function for GWCRP method, and returns
the posterior samples of the MCMC.

`getDahl.R` obtains the posterior estimates of cluster memberships and 
model parameters by Dahl's method.

`LPML.R` calculates the LPML based on the posterior samples.

`simulation_data_generation.R` contains the function of generating the
simulated survival data and 
the function of computing the Hessian matrix.

`d_matrix.RData` is the matrix of graph distances between pairs of counties in
Louisiana.

`simulation_design1.R` gives the code for simulation study. 

`realdata_analysis.R` gives the code for analyzing the SEER respiratory cancer
data in Louisiana.