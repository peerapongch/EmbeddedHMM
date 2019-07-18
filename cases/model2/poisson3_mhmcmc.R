load('./data/model2_poisson3.RData')
load('./data/model2_poisson3_env.RData')

source('../../samplers/model2_mhmcmc.R')
N <- 2000
es <- c(0.2,0.8)
seed <- 1
system.time(
  mhmcmc_out <- mhmcmcModel2(ssm_poisson,N,es,seed=seed)
)
# save(ehmm_out,file='./data/ehmm_backward_model1_poisson1_L50_N6000_AR_SH.RData')

library(coda)
x <- as.mcmc(mhmcmc_out$X_sample[,100,1:3])
plot(x)

source('../../evaluation/function_mcmc_mu.R')
mhmcmc_mu <- mci_mu(mhmcmc_out$X_sample,T,dim)

plot(ssm_poisson$X[,1],type='l')
lines(mhmcmc_mu[,1],col='red')
