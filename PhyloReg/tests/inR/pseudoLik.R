library(Rcpp)
library(RcppArmadillo)
library(roptim)
library(Matrix)
sourceCpp("./src//Ising.cpp") 
sourceCpp("./src/pseudolike.cpp")

n_para <- 100
n_host <- 100

p_para <- 2
p_host <- 2

set.seed(43)
beta_host <- c(-0,0.1)
beta_para <- c(-0.1,-0.1)

gamma_host <- 0.08
gamma_para <- 0.08

host_para_mat <- matrix(ifelse(runif(n_para*n_host)>0.5,1,-1),n_host,n_para)
design_para <- cbind(1,rnorm(n_para))
design_host <- cbind(rnorm(n_host),rnorm(n_host))

dist_host <- matrix(rexp(n_host^2),n_host,n_host)
dist_host <- dist_host+t(dist_host)
diag(dist_host) <- 0
dist_host <- dist_host/max(dist_host)

dist_para <- matrix(rexp(n_para^2),n_para,n_para)
dist_para <- dist_para+t(dist_para)
diag(dist_para) <- 0
dist_para <- dist_para/max(dist_para)


data_mcmc <- bootstrap_helper(beta_host, beta_para, gamma_host, gamma_para,
                              host_para_mat, dist_host, dist_para, design_host,
                              design_para, 1000, 100,10)


host_para_mat <- matrix(data_mcmc[2,],n_host,n_para)
template <- matrix(data_mcmc[1,],n_host,n_para)
image(host_para_mat)

fitting <- maxPL_cpp(host_para_mat, dist_host,dist_para,design_host,design_para,"BFGS")


bootstrap_mcmc <- bootstrap_helper(fitting$beta_host, fitting$beta_para, 
                              fitting$gamma_host, fitting$gamma_para,
                              template, dist_host, dist_para, design_host,
                              design_para, 1000, 10000,50)

all_boot <- lapply(1:200,function(i,boot_mcmc,n_host,n_para){
    matrix(boot_mcmc[i,],n_host,n_para)
},bootstrap_mcmc,n_host,n_para)

all_fit <- lapply(all_boot, maxPL_cpp, dist_host,dist_para,design_host,design_para,"BFGS")

hist(sapply(all_fit,function(w){w$gamma_para}))
