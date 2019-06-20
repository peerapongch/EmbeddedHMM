load('./data/model1_poisson1.RData')
load('./data/model1_poisson1_env.RData')
source('../../samplers/model1_mhmcmc.R')

### MHMCMC with autoregressive update ###
N <- 2800000
es <- c(0.2,0.8)
seed <- 1

system.time(
  mcmc_out <- mcmcGaussianSSM(N,es,ssm_poisson,obs='Poisson',seed=seed) 
)

save(mcmc_out,file='./data/test1_poisson1_mcmc_N2800000_seed1.RData')