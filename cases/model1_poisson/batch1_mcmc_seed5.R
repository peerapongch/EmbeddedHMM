load('./data/model1_poisson2.RData')
load('./data/model1_poisson2_env.RData')
source('../../samplers/model1_mhmcmc.R')

es <- c(0.2,0.8)
N <- 1000000
seed <- 5

s<- system.time(
  mcmc_out <- mcmcGaussianSSM(N,es,ssm_poisson,thin.factor=10,seed=seed)
)
tpi <- s[[3]]/N
tps <- tpi*10

save(tpi,tps,mcmc_out,file=paste('./data/batch1_mcmc_seed',seed,'.RData',sep=''))