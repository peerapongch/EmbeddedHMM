load('./data/model2_poisson3.RData')
load('./data/model2_poisson3_env.RData')

source('../../samplers/model2_pgmet.R')
N <- 1000
es <- c(0.2,0.8)
L <- 1000
seed <- 1
system.time(
  pgbs_out <- pgbsModel2(ssm_poisson,N,L,seed=seed)
)
# save(ehmm_out,file='./data/ehmm_backward_model1_poisson1_L50_N6000_AR_SH.RData')

library(coda)
x <- as.mcmc(pgbs_out$X_sample[,100,1:3])
plot(x)

source('../../evaluation/function_mcmc_mu.R')
pgbs_mu <- mci_mu(abs(pgbs_out$X_sample),T,dim)

plot(abs(ssm_poisson$X[,1]),type='l')
lines(pgbs_mu[,1],col='red')
