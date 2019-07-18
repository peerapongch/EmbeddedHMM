setMKLthreads(2)
load('./data/model2_poisson3.RData')
load('./data/model2_poisson3_env.RData')

source('../../samplers/model2_pgmet.R')
N <- 250
es <- c(0.3,1)
L <- 5000
N.mcmc <- 50
seed <- 1
system.time(
  pgmet_out <- pgmetModel2(ssm_poisson,N,L,es,N.mcmc=N.mcmc,seed=seed)
)
save(pgmet_out,file='./data/poisson3_pgmet.RData')

# library(coda)
# x <- as.mcmc(pgbs_out$X_sample[,100,1:3])
# plot(x)
# 
# source('../../evaluation/function_mcmc_mu.R')
# pgbs_mu <- mci_mu(abs(pgbs_out$X_sample),T,dim)
# 
# plot(abs(ssm_poisson$X[,1]),type='l')
# lines(pgbs_mu[,1],col='red')
